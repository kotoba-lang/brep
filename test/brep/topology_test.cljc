(ns brep.topology-test
  (:require [clojure.test :refer [deftest is testing]]
            [brep.topology :as topo]
            [brep.mesh-csg :as csg]
            [brep.feature :as f]))

(defn- square [id plane x0 y0 x1 y1]
  (f/sketch-feature id plane
                    [(f/sketch-line [x0 y0] [x1 y0]) (f/sketch-line [x1 y0] [x1 y1])
                     (f/sketch-line [x1 y1] [x0 y1]) (f/sketch-line [x0 y1] [x0 y0])]))

(defn- box-mesh []
  (second (f/evaluate-mesh
           (-> (f/feature-tree)
               (f/add-feature (square 1 (f/sketch-plane-xy) 0 0 4 4))
               (f/add-feature (f/extrude-feature 2 1 [0 0 1] 3 :new))))))

(deftest weld-reports-what-it-merged
  (let [w (topo/weld-mesh {:positions [[0 0 0] [0 0 0] [1 0 0] [0 1 0]]
                           :indices [0 2 3 1 2 3]})]
    (testing "coincident positions collapse"
      (is (= 3 (count (:positions w))))
      (is (= 1 (:merged w))))
    (testing "the two triangles survive because they were distinct"
      (is (= 2 (/ (count (:indices w)) 3)))
      (is (= 0 (:degenerate w)))))

  (testing "a triangle that loses a corner to the weld is dropped AND counted"
    ;; A weld that silently deletes triangles turns a hole into a 'clean' mesh.
    (let [w (topo/weld-mesh {:positions [[0 0 0] [0 0 0] [1 0 0]]
                             :indices [0 1 2]})]
      (is (= 0 (count (:indices w))))
      (is (= 1 (:degenerate w))))))

(deftest box-topology-is-a-closed-orientable-solid
  (let [t (topo/topology (box-mesh))]
    (testing "the soup welds down to the box's own 8 corners"
      (is (= 8 (count (:vertices t))))
      (is (= 24 (+ (count (:vertices t)) (:merged t)))))
    (testing "12 triangles, 18 edges, Euler characteristic 2"
      (is (= 12 (count (:faces t))))
      (is (= 18 (count (:edges t))))
      (is (= 2 (topo/euler-characteristic t))))
    (testing "closed and orientable"
      (is (:manifold? t))
      (is (:orientable? t))
      (is (empty? (topo/boundary-edges t))))))

(deftest convexity-needs-consistently-wound-faces
  ;; `tessellate-solid` emits each face independently, so two triangles meeting
  ;; at an edge may wind the SAME way; their computed normals then disagree by a
  ;; sign that has nothing to do with the shape. Before `orient-faces`, a plain
  ;; box measured 8 convex / 4 concave / 6 flat. A box has 12 convex edges (its
  ;; own), 6 flat ones (the diagonals splitting each square face) and no concave
  ;; edge anywhere.
  (let [t (topo/topology (box-mesh))
        by (fn [c] (count (topo/edges-where t #(= c (:convexity %)))))]
    (is (= 12 (by :convex)))
    (is (= 6 (by :flat)))
    (is (= 0 (by :concave)))
    (testing "the box's own edges are right angles"
      (is (every? #(< (Math/abs (- (:dihedral (get (:edges t) %)) (/ Math/PI 2))) 1.0e-9)
                  (topo/sharp-edges t :convex 1.0))))))

(deftest open-meshes-are-reported-open
  (let [t (topo/topology {:positions [[0 0 0] [1 0 0] [0 1 0]] :indices [0 1 2]})]
    (testing "a lone triangle has three boundary edges and is not manifold"
      (is (= 3 (count (topo/boundary-edges t))))
      (is (not (:manifold? t))))
    (testing "a boundary edge carries no dihedral or convexity"
      ;; On a boundary edge the question has no single answer, and answering it
      ;; anyway is how a chamfer ends up cutting into a hole.
      (is (every? #(and (nil? (:dihedral %)) (nil? (:convexity %)))
                  (vals (:edges t)))))))

(defn- box-mesh-at [x0 y0 x1 y1 z h]
  (let [m (second (f/evaluate-mesh
                   (-> (f/feature-tree)
                       (f/add-feature (square 1 (f/sketch-plane-xy) x0 y0 x1 y1))
                       (f/add-feature (f/extrude-feature 2 1 [0 0 1] h :new)))))]
    (update m :positions (fn [ps] (mapv (fn [[x y zz]] [x y (+ zz z)]) ps)))))

(defn- closure [op a b]
  (let [t (topo/topology (csg/mesh-boolean op a b))]
    {:euler (topo/euler-characteristic t)
     :boundary (count (topo/boundary-edges t))
     :closed? (and (:manifold? t) (empty? (topo/boundary-edges t)))}))

(deftest booleans-close-when-nothing-is-cut
  ;; `mesh-boolean` now welds and re-winds its inputs. It has to: a BSP boolean
  ;; reads front/back off each polygon's own normal, and `tessellate-solid`
  ;; emits faces independently with no guarantee of consistent winding, so every
  ;; caller arriving from a feature tree was handing it exactly that.
  ;;
  ;; The sharpest case is the first one: two boxes standing 20 apart, nothing to
  ;; intersect, and before this the union came back with 4 boundary edges.
  (let [big (box-mesh-at 0 0 10 10 0 4)]
    (testing "a union of two boxes that never meet is two closed shells"
      (let [r (closure :union big (box-mesh-at 20 20 24 24 0 4))]
        (is (:closed? r))
        (is (= 4 (:euler r)))))

    (testing "a union of two boxes sharing a face is one closed solid"
      (is (:closed? (closure :union big (box-mesh-at 10 0 20 10 0 4)))))

    (testing "a difference by a body that misses entirely leaves the original closed"
      (is (:closed? (closure :difference big (box-mesh-at 20 20 24 24 0 4)))))

    (testing "an intersection of two overlapping boxes is closed"
      (is (:closed? (closure :intersection big (box-mesh-at 5 5 15 15 0 4)))))))

(deftest booleans-that-actually-cut-are-still-open
  ;; Pinned, not accepted. Measured 2026-08-21 AFTER the winding fix: the cases
  ;; where the two bodies genuinely cut each other still leak. Part of the damage
  ;; is T-junctions (a vertex sitting in the middle of a neighbouring triangle's
  ;; edge — geometrically closed, topologically cracked) and the larger part is
  ;; polygons that are simply gone:
  ;;
  ;;   union, overlapping    18 boundary edges   6 T-junctions, 12 holes
  ;;   difference, hole      28 boundary edges   8 T-junctions, 20 holes
  ;;   difference, notch     14 boundary edges   6 T-junctions,  8 holes
  ;;
  ;; The csg.js transcription has been checked against the original step by step
  ;; and matches. The leading suspect is that this code is fed TRIANGLES where
  ;; csg.js is fed quads, so each planar face arrives pre-split and its fragments
  ;; fall below the three-vertex floor in `split-polygon`.
  ;;
  ;; These assertions say "still broken" on purpose. The day the CSG is fixed
  ;; they fail, and that failure is the notification.
  (let [big (box-mesh-at 0 0 10 10 0 4)]
    (testing "overlapping union leaks"
      (is (not (:closed? (closure :union big (box-mesh-at 5 5 15 15 0 4))))))
    (testing "a through hole leaks"
      (is (not (:closed? (closure :difference big (box-mesh-at 3 3 7 7 -1 6))))))
    (testing "a corner notch leaks"
      (is (not (:closed? (closure :difference big (box-mesh-at 8 8 12 12 -1 6))))))
    (testing "and the extrusion they all start from is closed, so the pipeline is not at fault"
      (is (:manifold? (topo/topology (box-mesh-at 0 0 10 10 0 4)))))))

(deftest coplanar-merge-puts-the-faces-back-together
  ;; csg.js takes a cube as six quads. A tessellated solid arrives as twelve
  ;; triangles, so every plane that should cut one polygon cuts two, and the
  ;; extra fragments fall below the three-vertex floor in the BSP's splitter.
  ;; This is the inverse operation, applied before a boolean.
  (let [polys (topo/merge-coplanar (box-mesh))]
    (testing "12 triangles become the box's own 6 quads"
      (is (= 6 (count polys)))
      (is (= [4 4 4 4 4 4] (mapv count polys))))

    (testing "every merged loop winds with the faces it came from"
      ;; `region-loop` chains border edges without regard to winding, so half the
      ;; loops come back reversed unless they are re-oriented — which would
      ;; reintroduce at the polygon level the defect `orient-faces` removes at
      ;; the face level. Leaving it out turned three closed boolean cases open.
      (let [centre [2.0 2.0 1.5]
            outward? (fn [pts]
                       (let [n (topo/-loop-normal pts)]
                         (pos? (reduce + (map * n (mapv - (first pts) centre))))))]
        (is (every? outward? polys))))

    (testing "a region whose boundary is not one simple loop stays triangulated"
      ;; Two triangles meeting only at a point share no edge, so they are two
      ;; regions, not one polygon with a pinch.
      (let [pinched {:positions [[0 0 0] [1 0 0] [0 1 0] [-1 0 0] [0 -1 0]]
                     :indices [0 1 2 0 3 4]}]
        (is (= 2 (count (topo/merge-coplanar pinched))))))))
