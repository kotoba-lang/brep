(ns brep.topology-test
  (:require [clojure.test :refer [deftest is testing]]
            [brep.topology :as topo]
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

(deftest mesh-booleans-do-not-return-closed-solids
  ;; Measured 2026-08-21. `brep.mesh-csg/mesh-boolean` renders, and its triangle
  ;; count changes as expected, but the result is not watertight in ANY case
  ;; tried — including a union of two boxes that do not touch, where there is
  ;; nothing to intersect. Pinned here so the day it is fixed, this test says so.
  ;;
  ;; Stable across weld tolerances from 1e-9 to 1e-2, so this is real geometry,
  ;; not a welding artifact.
  (let [run (fn [tree] (topo/topology (second (f/evaluate-mesh tree))))
        base (-> (f/feature-tree)
                 (f/add-feature (square 1 (f/sketch-plane-xy) 0 0 10 10))
                 (f/add-feature (f/extrude-feature 2 1 [0 0 1] 4 :new)))]
    (testing "the extrusion it starts from IS closed, so the pipeline is not at fault"
      (is (:manifold? (run base))))

    (testing "a corner cut leaves the mesh open"
      (let [t (run (-> base
                       (f/add-feature (square 3 (f/sketch-plane-xy) 5 5 12 12))
                       (f/add-feature (f/extrude-feature 4 3 [0 0 1] 9 :cut))))]
        (is (not (:manifold? t)))
        (is (pos? (count (topo/boundary-edges t))))))

    (testing "so does a union of two boxes that never meet"
      (let [t (run (-> (f/feature-tree)
                       (f/add-feature (square 1 (f/sketch-plane-xy) 0 0 4 4))
                       (f/add-feature (f/extrude-feature 2 1 [0 0 1] 2 :new))
                       (f/add-feature (square 3 (f/sketch-plane-xy) 20 20 24 24))
                       (f/add-feature (f/extrude-feature 4 3 [0 0 1] 2 :add))))]
        (is (not (:manifold? t)))
        (is (pos? (count (topo/boundary-edges t))))))))
