(ns brep.mesh-csg
  "Portable BSP-based boolean operations for closed triangle meshes."
  (:require [brep.topology :as topo]))

(def ^:private epsilon 1.0e-7)
(def ^:private coplanar 0)
(def ^:private front 1)
(def ^:private back 2)
(def ^:private spanning 3)

(defn- v+ [a b] (mapv + a b))
(defn- v- [a b] (mapv - a b))
(defn- scale [v s] (mapv #(* s %) v))
(defn- dot [a b] (reduce + (map * a b)))
(defn- cross [[ax ay az] [bx by bz]]
  [(- (* ay bz) (* az by)) (- (* az bx) (* ax bz)) (- (* ax by) (* ay bx))])
(defn- normalize [v]
  (let [length (#?(:clj Math/sqrt :cljs js/Math.sqrt) (dot v v))]
    (if (pos? length) (scale v (/ 1.0 length)) [0.0 0.0 0.0])))

(defn- plane-from-vertices [vertices]
  (let [[a b c] (map :position (take 3 vertices))
        normal (normalize (cross (v- b a) (v- c a)))]
    {:normal normal :w (dot normal a)}))

(defn- polygon [vertices]
  {:vertices (vec vertices) :plane (plane-from-vertices vertices)})

(defn- flipped-polygon [poly]
  (let [vertices (mapv #(update % :normal (fn [normal] (mapv - normal)))
                       (reverse (:vertices poly)))
        plane (:plane poly)]
    {:vertices vertices :plane {:normal (mapv - (:normal plane)) :w (- (:w plane))}}))

(defn- interpolate-vertex [a b t]
  {:position (v+ (:position a) (scale (v- (:position b) (:position a)) t))
   :normal (normalize (v+ (:normal a) (scale (v- (:normal b) (:normal a)) t)))})

(defn- split-polygon [plane poly]
  (let [types (mapv (fn [vertex]
                      (let [distance (- (dot (:normal plane) (:position vertex)) (:w plane))]
                        (cond (< distance (- epsilon)) back
                              (> distance epsilon) front
                              :else coplanar)))
                    (:vertices poly))
        polygon-type (reduce bit-or coplanar types)]
    (case polygon-type
      0 (if (pos? (dot (:normal plane) (get-in poly [:plane :normal])))
          {:coplanar-front [poly]} {:coplanar-back [poly]})
      1 {:front [poly]}
      2 {:back [poly]}
      3 (let [vertices (:vertices poly) count-vertices (count vertices)
              [front-vertices back-vertices]
              (reduce (fn [[front-result back-result] i]
                        (let [j (mod (inc i) count-vertices)
                              ti (nth types i) tj (nth types j)
                              vi (nth vertices i) vj (nth vertices j)
                              front-result (cond-> front-result (not= ti back) (conj vi))
                              back-result (cond-> back-result (not= ti front) (conj vi))]
                          (if (= spanning (bit-or ti tj))
                            (let [direction (v- (:position vj) (:position vi))
                                  t (/ (- (:w plane) (dot (:normal plane) (:position vi)))
                                       (dot (:normal plane) direction))
                                  vertex (interpolate-vertex vi vj t)]
                              [(conj front-result vertex) (conj back-result vertex)])
                            [front-result back-result])))
                      [[] []] (range count-vertices))]
          (cond-> {}
            (>= (count front-vertices) 3) (assoc :front [(polygon front-vertices)])
            (>= (count back-vertices) 3) (assoc :back [(polygon back-vertices)]))))))

(defn- partition-polygons [plane polygons]
  (reduce (fn [result poly]
            (merge-with into result (split-polygon plane poly)))
          {:coplanar-front [] :coplanar-back [] :front [] :back []} polygons))

(declare build-node clip-to)

(defn- build-node
  ([polygons] (build-node nil polygons))
  ([node polygons]
   (if (empty? polygons)
     node
     (let [plane (or (:plane node) (:plane (first polygons)))
           parts (partition-polygons plane polygons)]
       {:plane plane
        :polygons (vec (concat (:polygons node) (:coplanar-front parts)
                               (:coplanar-back parts)))
        :front (build-node (:front node) (:front parts))
        :back (build-node (:back node) (:back parts))}))))

(defn- clip-polygons [node polygons]
  (if-not (:plane node)
    polygons
    (let [parts (partition-polygons (:plane node) polygons)
          front-polygons (vec (concat (:coplanar-front parts) (:front parts)))
          back-polygons (vec (concat (:coplanar-back parts) (:back parts)))
          front-polygons (if (:front node) (clip-polygons (:front node) front-polygons)
                             front-polygons)
          back-polygons (if (:back node) (clip-polygons (:back node) back-polygons) [])]
      (vec (concat front-polygons back-polygons)))))

(defn- clip-to [node other]
  (when node
    (assoc node
           :polygons (clip-polygons other (:polygons node))
           :front (clip-to (:front node) other)
           :back (clip-to (:back node) other))))

(defn- invert-node [node]
  (when node
    (let [plane (:plane node)]
      {:plane {:normal (mapv - (:normal plane)) :w (- (:w plane))}
       :polygons (mapv flipped-polygon (:polygons node))
       :front (invert-node (:back node))
       :back (invert-node (:front node))})))

(defn- all-polygons [node]
  (if node
    (vec (concat (:polygons node) (all-polygons (:front node)) (all-polygons (:back node))))
    []))

(defn- mesh->polygons [{:keys [positions indices normals]}]
  (mapv (fn [[a b c]]
          (let [face-normal (normalize (cross (v- (nth positions b) (nth positions a))
                                              (v- (nth positions c) (nth positions a))))]
            (polygon (mapv (fn [index]
                             {:position (nth positions index)
                              :normal (if (seq normals) (nth normals index) face-normal)})
                           [a b c]))))
        (partition 3 indices)))

(defn- polygons->mesh [polygons]
  (let [triangles (mapcat (fn [poly]
                            (let [vertices (:vertices poly)]
                              (map (fn [i] [(first vertices) (nth vertices i) (nth vertices (inc i))])
                                   (range 1 (dec (count vertices))))))
                          polygons)
        vertices (vec (mapcat identity triangles))]
    {:positions (mapv :position vertices)
     :normals (mapv :normal vertices)
     :indices (vec (range (count vertices)))}))

(defn mesh-boolean
  "Apply `:union`, `:difference`, or `:intersection` to two closed triangle
  meshes.

  **The inputs are welded and re-wound before use.** A BSP boolean classifies
  every polygon as front or back of a plane, and that classification is read off
  the polygon's own normal — so a mesh whose faces are not wound consistently
  feeds the tree normals that point in directions unrelated to the shape.
  `brep.tessellate/tessellate-solid` emits each face independently and does NOT
  guarantee consistent winding, so every caller coming from a feature tree was
  handing this function exactly that. Measured 2026-08-21: a union of two boxes
  standing 20 apart — nothing to intersect at all — came back with 4 boundary
  edges. Welding and orienting first makes that case, and a union of two boxes
  sharing a face, come back closed.

  ⚠ **The result is still not watertight in every case.** Booleans whose
  operands actually cut each other leak, measured as of 2026-08-21:

      union, disjoint / touching        closed
      union, overlapping                18 boundary edges (6 T-junctions, 12 holes)
      difference, no overlap            closed
      difference, through hole          28 boundary edges (8 T-junctions, 20 holes)
      difference, corner notch          14 boundary edges (6 T-junctions, 8 holes)
      intersection, overlapping         closed

  So part of the damage is T-junctions — vertices sitting in the middle of a
  neighbouring triangle's edge, where the surface is geometrically closed but
  topologically cracked — and the larger part is polygons that are simply gone.
  The transcription of the csg.js algorithm here has been checked against the
  original step by step and matches it; the leading suspect is that this
  function is fed TRIANGLES where csg.js is fed quads, so each planar face
  arrives already split and its fragments fall below the three-vertex floor in
  `split-polygon`. Feeding it coplanar-merged polygons is the next thing to try.
  `brep.topology` is what can both merge them and check the answer.

  Until then: check the result with `brep.topology/topology` before relying on
  it being a solid. Volume, STEP export, machining stock and printing all
  depend on closure, and this does not yet provide it in the cutting cases."
  [operation mesh-a mesh-b]
  (let [a (build-node (mesh->polygons (topo/welded-oriented mesh-a)))
        b (build-node (mesh->polygons (topo/welded-oriented mesh-b)))
        result
        (case operation
          :union (let [a1 (clip-to a b)
                       b1 (clip-to b a1)
                       b2 (-> b1 invert-node (clip-to a1) invert-node)]
                   (build-node a1 (all-polygons b2)))
          :difference (let [a1 (invert-node a)
                            a2 (clip-to a1 b)
                            b1 (clip-to b a2)
                            b2 (-> b1 invert-node (clip-to a2) invert-node)]
                        (-> (build-node a2 (all-polygons b2)) invert-node))
          :intersection (let [a1 (invert-node a)
                              b1 (clip-to b a1)
                              b2 (invert-node b1)
                              a2 (clip-to a1 b2)
                              b3 (clip-to b2 a2)]
                          (-> (build-node a2 (all-polygons b3)) invert-node)))]
    (polygons->mesh (all-polygons result))))
