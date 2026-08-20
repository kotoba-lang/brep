(ns brep.feature-eval-test
  "Verification for the 2026-07-27 additions to `brep.feature` and
  `brep.assembly`:

  - `sketch->ring` / `extrude-prism` — an extrude now consumes its sketch
    profile instead of producing a unit square regardless of it.
  - `evaluate` — boolean features return an explicit error instead of the
    previous silent no-op.
  - `evaluate-mesh` — booleans applied for real via `brep.mesh-csg`.
  - `brep.assembly/solve-translations` — axis-aligned translation constraints
    actually move instances.

  Every geometric assertion compares against a closed-form value (a prism's
  volume is its footprint times its height; a box with a rectangular hole is
  the difference of two such products), so these are verification tests, not
  change detectors."
  (:require [clojure.test :refer [deftest is testing]]
            [clojure.string :as string]
            [brep.feature :as f]
            [brep.assembly :as a]
            [brep.kernel :as k]
            [brep.step :as step]
            [brep.tessellate :as tess]))

(defn- rect-sketch [id w h]
  (f/sketch-feature id (f/sketch-plane-xy)
                    [(f/sketch-line [0 0] [w 0])
                     (f/sketch-line [w 0] [w h])
                     (f/sketch-line [w h] [0 h])
                     (f/sketch-line [0 h] [0 0])]))

(defn- offset-rect-sketch [id x0 y0 w h]
  (f/sketch-feature id (f/sketch-plane-xy)
                    [(f/sketch-line [x0 y0] [(+ x0 w) y0])
                     (f/sketch-line [(+ x0 w) y0] [(+ x0 w) (+ y0 h)])
                     (f/sketch-line [(+ x0 w) (+ y0 h)] [x0 (+ y0 h)])
                     (f/sketch-line [x0 (+ y0 h)] [x0 y0])]))

(defn- mesh-volume
  "Signed volume of a closed triangle mesh via the divergence theorem."
  [{:keys [positions indices]}]
  (Math/abs
   (double
    (reduce
     + 0.0
     (for [t (range (quot (count indices) 3))]
       (let [g (fn [kk] (nth positions (nth indices (+ (* t 3) kk))))
             [ax ay az] (g 0) [bx by bz] (g 1) [cx cy cz] (g 2)]
         (/ (+ (* ax (- (* by cz) (* bz cy)))
               (- (* ay (- (* bx cz) (* bz cx))))
               (* az (- (* bx cy) (* by cx))))
            6.0)))))))

;; ---------------------------------------------------------------------------
;; sketch->ring
;; ---------------------------------------------------------------------------

(deftest sketch-ring-chains-segments-in-any-order
  (testing "a closed rectangle chains regardless of segment order/direction"
    (let [sk (f/sketch-feature "s" (f/sketch-plane-xy)
                               [(f/sketch-line [10 0] [10 10])   ; deliberately not
                                (f/sketch-line [0 10] [0 0])     ; in walk order,
                                (f/sketch-line [10 10] [0 10])   ; and some reversed
                                (f/sketch-line [0 0] [10 0])])
          ring (f/sketch->ring sk)]
      (is (= 4 (count ring)))
      (is (= #{[0 0] [10 0] [10 10] [0 10]} (set ring))))))

(deftest sketch-ring-refuses-open-and-multi-loop
  (testing "an open chain yields nil, never a guess"
    (is (nil? (f/sketch->ring
               (f/sketch-feature "s" (f/sketch-plane-xy)
                                 [(f/sketch-line [0 0] [10 0])
                                  (f/sketch-line [10 0] [10 10])
                                  (f/sketch-line [10 10] [5 20])])))))
  (testing "two disjoint closed loops yield nil"
    (is (nil? (f/sketch->ring
               (f/sketch-feature "s" (f/sketch-plane-xy)
                                 (concat (:entities (rect-sketch "a" 10 10))
                                         (:entities (offset-rect-sketch "b" 100 100 5 5)))))))))

(deftest sketch-ring-polygonises-circles
  (let [sk (f/sketch-feature "s" (f/sketch-plane-xy) [(f/sketch-circle [0 0] 10)])
        ring (f/sketch->ring sk 16)]
    (is (= 16 (count ring)) "circle becomes the requested number of segments")))

;; ---------------------------------------------------------------------------
;; extrude consumes the profile
;; ---------------------------------------------------------------------------

(deftest extrude-uses-the-sketch-not-a-unit-square
  (testing "prism volume = footprint x height, exactly"
    (let [w 305.4 h 192.4 d 133.4
          tree (-> (f/feature-tree)
                   (f/add-feature (rect-sketch "s1" w h))
                   (f/add-feature (f/extrude-feature "e1" "s1" [0 0 1] d :new)))
          [st [solid edges verts]] (f/evaluate tree)]
      (is (= :ok st))
      (is (= 6 (k/face-count solid)) "a 4-gon prism has 4 sides + 2 caps")
      (is (= 12 (k/edge-count solid)))
      (is (= 8 (k/vertex-count solid edges verts)))
      (is (= [[0.0 0.0 0.0] [w h d]]
             (mapv #(mapv double %) (k/bounding-box solid edges verts))))
      (is (< (Math/abs (- (tess/volume solid edges verts) (* w h d))) 1e-6)
          "this is the regression that matters: the old evaluator returned a 1x1 box"))))

(deftest extruded-prism-is-step-exportable
  (testing "planar faces + linear edges survive a STEP round-trip"
    (let [tree (-> (f/feature-tree)
                   (f/add-feature (rect-sketch "s1" 60 40))
                   (f/add-feature (f/extrude-feature "e1" "s1" [0 0 1] 20 :new)))
          [_ [solid edges verts]] (f/evaluate tree)
          text (step/write-step solid edges verts)
          [s2 e2 v2] (step/read-step text)]
      (is (re-find #"^ISO-10303-21;" text))
      (is (= (k/face-count solid) (k/face-count s2)))
      (is (< (Math/abs (- (tess/volume solid edges verts) (tess/volume s2 e2 v2))) 1e-6)
          "volume must survive write->read"))))

(deftest extrude-refuses-an-unusable-named-sketch
  (testing "a named-but-open sketch errors rather than silently becoming 1x1"
    (let [tree (-> (f/feature-tree)
                   (f/add-feature (f/sketch-feature "s1" (f/sketch-plane-xy)
                                                    [(f/sketch-line [0 0] [10 0])]))
                   (f/add-feature (f/extrude-feature "e1" "s1" [0 0 1] 5 :new)))
          [st msg] (f/evaluate tree)]
      (is (= :error st))
      (is (re-find #"closed loop" msg)))))

(deftest extrude-without-a-sketch-keeps-legacy-unit-square
  (testing "no :sketch-ref at all still yields the documented unit prism"
    (let [tree (-> (f/feature-tree)
                   (f/add-feature (f/extrude-feature "e1" nil [0 0 1] 3 :new)))
          [st [solid edges verts]] (f/evaluate tree)]
      (is (= :ok st))
      (is (< (Math/abs (- (tess/volume solid edges verts) 3.0)) 1e-9)
          "1 x 1 x 3"))))

;; ---------------------------------------------------------------------------
;; booleans
;; ---------------------------------------------------------------------------

(deftest boolean-in-evaluate-is-an-explicit-error-not-a-no-op
  (testing "the old behaviour silently returned the base solid unchanged"
    (let [tree (-> (f/feature-tree)
                   (f/add-feature (rect-sketch "s1" 100 100))
                   (f/add-feature (f/extrude-feature "e1" "s1" [0 0 1] 100 :new))
                   (f/add-feature (offset-rect-sketch "s2" 25 25 50 50))
                   (f/add-feature (f/extrude-feature "e2" "s2" [0 0 1] 100 :cut)))
          [st msg] (f/evaluate tree)]
      (is (= :error st))
      (is (re-find #"evaluate-mesh" msg) "the error must name the working path"))))

(deftest evaluate-mesh-applies-a-real-difference
  (testing "cutting a 50x50x100 hole from a 100x100x100 block removes 0.25 of it"
    (let [base (-> (f/feature-tree)
                   (f/add-feature (rect-sketch "s1" 100 100))
                   (f/add-feature (f/extrude-feature "e1" "s1" [0 0 1] 100 :new)))
          cut (-> base
                  (f/add-feature (offset-rect-sketch "s2" 25 25 50 50))
                  (f/add-feature (f/extrude-feature "e2" "s2" [0 0 1] 100 :cut)))
          [st1 m1] (f/evaluate-mesh base)
          [st2 m2] (f/evaluate-mesh cut)]
      (is (= :ok st1))
      (is (= :ok st2))
      (is (< (Math/abs (- (mesh-volume m1) 1.0e6)) 1.0) "base = 100^3 mm3")
      (is (< (Math/abs (- (mesh-volume m2) 0.75e6)) 1.0)
          "after the cut = 750,000 mm3 — proves the boolean is real"))))

;; ---------------------------------------------------------------------------
;; assembly translation solver
;; ---------------------------------------------------------------------------

(defn- two-part-assembly []
  (let [asm (a/assembly "t")
        [i1 asm] (a/add-instance asm (a/part-ref 1 "base") {:translate [0.0 0.0 0.0]} "base")
        [i2 asm] (a/add-instance asm (a/part-ref 2 "arm") {:translate [0.0 0.0 0.0]} "arm")]
    [i1 i2 asm]))

(deftest distance-constraint-actually-moves-the-instance
  (let [[i1 i2 asm] (two-part-assembly)
        asm (a/add-constraint asm (assoc (a/distance-constraint i1 "f" i2 "f" 25.0) :axis 2))
        r (a/solve-translations asm)]
    (is (= :ok (:status r)))
    (is (= [[0.0 0.0 0.0] [0.0 0.0 25.0]]
           (mapv (comp :translate :transform) (a/instances (:assembly r))))
        "solve must change geometry, not just validate references")))

(deftest mate-constraint-makes-instances-coincident-on-the-axis
  (let [asm (a/assembly "t")
        [i1 asm] (a/add-instance asm (a/part-ref 1 "a") {:translate [0.0 0.0 7.0]} "a")
        [i2 asm] (a/add-instance asm (a/part-ref 2 "b") {:translate [0.0 0.0 40.0]} "b")
        asm (a/add-constraint asm (assoc (a/mate-constraint i1 "f" i2 "f") :axis 2))
        r (a/solve-translations asm)]
    (is (= :ok (:status r)))
    (is (= 7.0 (nth (:translate (:transform (second (a/instances (:assembly r))))) 2)))))

(deftest solver-refuses-what-it-cannot-do
  (testing "rotation constraints are refused, not ignored"
    (let [[i1 i2 asm] (two-part-assembly)
          asm (a/add-constraint asm (assoc (a/angle-constraint i1 "f" i2 "f" 1.5) :axis 0))
          r (a/solve-translations asm)]
      (is (= :error (:status r)))
      (is (re-find #"rotation is not implemented" (:message r)))))
  (testing ":identity transforms are refused"
    (let [asm (a/assembly "t")
          [_ asm] (a/add-instance asm (a/part-ref 1 "a") (a/affine3-identity) "a")
          [_ asm] (a/add-instance asm (a/part-ref 2 "b") (a/affine3-identity) "b")]
      (is (= :error (:status (a/solve-translations asm))))))
  (testing "a missing :axis is refused"
    (let [[i1 i2 asm] (two-part-assembly)
          asm (a/add-constraint asm (a/distance-constraint i1 "f" i2 "f" 5.0))]
      (is (= :error (:status (a/solve-translations asm)))))))

(deftest solver-reports-non-convergence-instead-of-swallowing-it
  (testing "a cyclic constraint pair hits the cap and says so"
    (let [[i1 i2 asm] (two-part-assembly)
          asm (-> asm
                  (a/add-constraint (assoc (a/distance-constraint i1 "f" i2 "f" 10.0) :axis 0))
                  (a/add-constraint (assoc (a/distance-constraint i2 "f" i1 "f" 10.0) :axis 0)))
          r (a/solve-translations asm 1e-9 8)]
      (is (= :unconverged (:status r)))
      (is (pos? (:residual r)))
      (is (re-find #"cyclic or" (:message r))))))

;; ---------------------------------------------------------------------------
;; The feature registry (2026-08-20).
;;
;; `apply-feature` is an open multimethod: registering a method IS the claim
;; that this kernel evaluates that feature kind. Before this, the constructors
;; above (`fillet-feature`, `sweep-feature`, …) read as capability while
;; `evaluate-mesh` handled `:extrude` alone and answered everything else with a
;; generic "not supported as a boolean operand" — a vocabulary and an evaluator
;; kept as two lists that nothing compared. The tests below pin the property
;; that replaces the second list: what is unregistered says so, by name.
;; ---------------------------------------------------------------------------

(defn- square-sketch [id plane s]
  (f/sketch-feature
   id plane
   [(f/sketch-line [0 0] [s 0]) (f/sketch-line [s 0] [s s])
    (f/sketch-line [s s] [0 s]) (f/sketch-line [0 s] [0 0])]))

(defn- bbox [{:keys [positions]}]
  [(mapv (fn [i] (apply min (map #(nth % i) positions))) [0 1 2])
   (mapv (fn [i] (apply max (map #(nth % i) positions))) [0 1 2])])

(deftest registry-answers-for-itself
  (testing "the supported set is derived from the registry, not restated"
    (is (contains? (f/supported-feature-kinds) :extrude))
    (is (contains? (f/supported-feature-kinds) :sketch))
    (is (not (contains? (f/supported-feature-kinds) :fillet))))

  (testing "an unregistered kind refuses BY NAME and names what is registered"
    (let [tree (-> (f/feature-tree)
                   (f/add-feature (square-sketch 1 (f/sketch-plane-xy) 4))
                   (f/add-feature (f/extrude-feature 2 1 [0 0 1] 3 :new))
                   (f/add-feature (f/fillet-feature 3 [1] 1.0)))
          [status msg] (f/evaluate-mesh tree)]
      (is (= :error status))
      (is (string/includes? msg ":fillet"))
      (is (string/includes? msg "no evaluator is registered"))
      (is (string/includes? msg ":extrude"))
      (is (string/includes? msg ":loft")))))

(deftest loft-builds-a-frustum-between-two-profiles
  (let [tree (-> (f/feature-tree)
                 (f/add-feature (square-sketch 1 (f/sketch-plane-xy) 6))
                 (f/add-feature
                  (square-sketch 2 (f/sketch-plane-custom [1.0 1.0 6.0] [0 0 1]) 4))
                 (f/add-feature (f/loft-feature 3 [1 2] :new)))
        [status mesh] (f/evaluate-mesh tree)]
    (testing "it builds"
      (is (= :ok status)))
    (testing "the solid spans both profile planes"
      (let [[[_ _ zmin] [_ _ zmax]] (bbox mesh)]
        (is (= 0.0 (double zmin)))
        (is (= 6.0 (double zmax)))))
    (testing "side walls plus both caps, all triangles"
      (is (zero? (mod (count (:indices mesh)) 3)))
      (is (= 36 (count (:indices mesh)))))
    (testing "every index addresses a position that exists"
      (is (every? #(< -1 % (count (:positions mesh))) (:indices mesh)))))

  (testing "mismatched vertex counts are refused, not silently resampled"
    (let [tree (-> (f/feature-tree)
                   (f/add-feature (square-sketch 1 (f/sketch-plane-xy) 6))
                   (f/add-feature
                    (f/sketch-feature 2 (f/sketch-plane-custom [0.0 0.0 5.0] [0 0 1])
                                            [(f/sketch-line [0 0] [4 0])
                                             (f/sketch-line [4 0] [2 4])
                                             (f/sketch-line [2 4] [0 0])]))
                   (f/add-feature (f/loft-feature 3 [1 2] :new)))
          [status msg] (f/evaluate-mesh tree)]
      (is (= :error status))
      (is (string/includes? msg "same vertex count"))
      (is (string/includes? msg "resampling is not implemented")))))

(deftest pattern-repeats-the-body-along-a-direction
  (let [base-tree (-> (f/feature-tree)
                      (f/add-feature (square-sketch 1 (f/sketch-plane-xy) 4))
                      (f/add-feature (f/extrude-feature 2 1 [0 0 1] 3 :new)))
        [_ base] (f/evaluate-mesh base-tree)
        [status mesh] (f/evaluate-mesh
                       (f/add-feature base-tree
                                            (f/pattern-feature 3 [2] [1 0 0] 3 10)))]
    (testing "it builds"
      (is (= :ok status)))
    (testing "the bounding box grows by (count - 1) x spacing along the direction"
      (let [[[x0 _ _] [x1 _ _]] (bbox base)
            [[px0 _ _] [px1 _ _]] (bbox mesh)]
        (is (= (double x0) (double px0)))
        (is (= (+ (double x1) 20.0) (double px1)))))
    (testing "the direction is normalised, so its length does not scale the spacing"
      (let [[_ scaled] (f/evaluate-mesh
                        (f/add-feature base-tree
                                             (f/pattern-feature 3 [2] [7 0 0] 3 10)))]
        (is (= (bbox mesh) (bbox scaled))))))

  (testing "degenerate patterns are refused"
    (let [t (-> (f/feature-tree)
                (f/add-feature (square-sketch 1 (f/sketch-plane-xy) 4))
                (f/add-feature (f/extrude-feature 2 1 [0 0 1] 3 :new)))]
      (is (= :error (first (f/evaluate-mesh
                            (f/add-feature t (f/pattern-feature 3 [2] [1 0 0] 1 10))))))
      (is (= :error (first (f/evaluate-mesh
                            (f/add-feature t (f/pattern-feature 3 [2] [0 0 0] 3 10))))))
      (is (= :error (first (f/evaluate-mesh
                            (f/add-feature t (f/pattern-feature 3 [2] [1 0 0] 3 0))))))))

  (testing "selective patterning is refused rather than over-applied"
    (let [t (-> (f/feature-tree)
                (f/add-feature (square-sketch 1 (f/sketch-plane-xy) 4))
                (f/add-feature (f/extrude-feature 2 1 [0 0 1] 3 :new))
                (f/add-feature (square-sketch 4 (f/sketch-plane-xy) 2))
                (f/add-feature (f/extrude-feature 5 4 [0 0 1] 9 :add))
                (f/add-feature (f/pattern-feature 6 [2] [1 0 0] 3 10)))
          [status msg] (f/evaluate-mesh t)]
      (is (= :error status))
      (is (string/includes? msg "Selective patterning is not implemented")))))
