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
