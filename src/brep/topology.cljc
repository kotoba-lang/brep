(ns brep.topology
  "Topology for a triangle mesh: welded vertices, face adjacency across edges,
  and per-edge dihedral angle and convexity.

  **Why this exists.** `brep.feature/evaluate-mesh` accumulates
  `{:positions :indices}` — a triangle soup with no edge identity and no face
  adjacency. Every feature that acts on EDGES rather than on profiles (chamfer,
  fillet, shell) needs exactly what a soup does not carry, and the feature
  constructors name their targets with BREP edge ids that the mesh side does
  not have. Mapping those ids onto soup edges by guesswork would return a
  successfully-built solid with the WRONG edges treated. This namespace supplies
  the missing half so those features can select edges by a property of the mesh
  they are actually modifying.

  Nothing here modifies geometry. It answers questions about it.

  All of it is pure `.cljc`, no I/O."
  (:require [brep.kernel :as k]
            [brep.config :as config]))

;; ---------------------------------------------------------------------------
;; Welding
;; ---------------------------------------------------------------------------

(defn- quantise
  "Snap a point to an integer lattice of cell size `eps`, as the weld key.

  A tolerance-based weld is not an equivalence relation — a and b may each be
  within eps of c while being 2*eps apart — so any implementation has to choose
  where to break the chain. Snapping to a lattice chooses it deterministically
  and in O(n): the answer never depends on the order the points arrive in,
  which a nearest-neighbour walk cannot promise. The cost is the usual one for
  lattices: two points either side of a cell boundary and closer than eps stay
  apart. `weld-mesh` reports the counts so a caller can see what happened
  rather than inferring it."
  [eps [x y z]]
  [(Math/round (/ (double x) eps))
   (Math/round (/ (double y) eps))
   (Math/round (/ (double z) eps))])

(defn weld-mesh
  "Merge coincident positions and reindex the triangles.

  Returns `{:positions :indices :merged :degenerate}` where `:merged` is how
  many input positions collapsed into another and `:degenerate` how many
  triangles lost a corner in the process and were dropped. Both are reported,
  never silent: a weld that quietly deletes triangles turns a hole into a
  'clean' mesh."
  ([mesh] (weld-mesh mesh (* 1.0e3 (get-in config/config [:epsilon :point-merge]))))
  ([{:keys [positions indices]} eps]
   (let [{:keys [pts key->id remap]}
         (reduce (fn [{:keys [pts key->id remap] :as acc} p]
                   (let [key (quantise eps p)]
                     (if-let [id (get key->id key)]
                       (assoc acc :remap (conj remap id))
                       (let [id (count pts)]
                         {:pts (conj pts p)
                          :key->id (assoc key->id key id)
                          :remap (conj remap id)}))))
                 {:pts [] :key->id {} :remap []}
                 positions)
         tris (partition 3 indices)
         kept (vec (for [[a b c] tris
                         :let [t [(remap a) (remap b) (remap c)]]
                         :when (= 3 (count (set t)))]
                     t))]
     {:positions pts
      :indices (vec (mapcat identity kept))
      :merged (- (count positions) (count pts))
      :degenerate (- (count tris) (count kept))})))

(defn- edge-key [a b] (if (< a b) [a b] [b a]))

(defn- signed-volume
  "Six times the signed volume of a closed triangle mesh. Positive when the
  faces wind counter-clockwise seen from outside, i.e. when normals point out."
  [pts faces]
  (reduce (fn [acc [a b c]]
            (+ acc (k/v-dot (nth pts a) (k/v-cross (nth pts b) (nth pts c)))))
          0.0 faces))

(defn orient-faces
  "Return `[faces orientable?]` with every face wound consistently, and with
  normals pointing OUT when the mesh is closed.

  Convexity cannot be read off a triangle soup as it arrives. `tessellate-solid`
  emits each face independently, so two triangles meeting at an edge may wind
  the same way rather than opposite ways, and their computed normals then
  disagree by a sign that has nothing to do with the shape. Measured on a plain
  extruded box: 8 edges read convex, 4 read concave and 0 read as the 12 convex
  edges a box has.

  So: flood-fill the adjacency, flipping any neighbour that traverses the shared
  edge in the same direction as the face it was reached from — that is what
  'consistently wound' means — then flip everything at once if the enclosed
  signed volume came out negative. `orientable?` is false for a Möbius-like
  mesh, where no consistent choice exists; callers must not read convexity off
  one."
  [pts faces]
  (let [n (count faces)
        edge->fs (reduce (fn [m [fi [a b c]]]
                           (reduce (fn [m e] (update m e (fnil conj []) fi))
                                   m [(edge-key a b) (edge-key b c) (edge-key c a)]))
                         {} (map-indexed vector faces))
        directed? (fn [face i j]
                    ;; true when `face` traverses edge i->j in that order
                    (let [[a b c] face]
                      (or (= [i j] [a b]) (= [i j] [b c]) (= [i j] [c a]))))
        flip (fn [[a b c]] [a c b])]
    (loop [oriented (vec faces)
           seen #{0}
           queue [0]
           ok true]
      (if (empty? queue)
        (let [remaining (remove seen (range n))]
          (if (seq remaining)
            ;; another connected shell: orient it independently
            (recur oriented (conj seen (first remaining)) [(first remaining)] ok)
            (let [v (signed-volume pts oriented)]
              [(if (neg? v) (mapv flip oriented) oriented) ok])))
        (let [fi (peek queue)
              face (nth oriented fi)
              [a b c] face
              nbrs (mapcat (fn [[i j]]
                             (map (fn [fj] [fj i j])
                                  (remove #{fi} (get edge->fs (edge-key i j) []))))
                           [[a b] [b c] [c a]])
              [oriented seen queue ok]
              (reduce (fn [[oriented seen queue ok] [fj i j]]
                        (let [nf (nth oriented fj)
                              same? (directed? nf i j)]
                          (cond
                            (seen fj) [oriented seen queue (and ok (not same?))]
                            same? [(assoc oriented fj (flip nf)) (conj seen fj) (conj queue fj) ok]
                            :else [oriented (conj seen fj) (conj queue fj) ok])))
                      [oriented seen (pop queue) ok]
                      nbrs)]
          (recur oriented seen queue ok))))))

;; ---------------------------------------------------------------------------
;; Adjacency
;; ---------------------------------------------------------------------------

(defn- face-normal [pts [a b c]]
  (k/v-normalize (k/v-cross (k/v- (nth pts b) (nth pts a))
                            (k/v- (nth pts c) (nth pts a)))))

(defn- edge-kind [n]
  (case (long n) 1 :boundary 2 :manifold :non-manifold))

(defn- opposite-corner
  "The corner of `face` that is not on `edge`."
  [face [i j]]
  (first (remove #{i j} face)))

(defn- classify-edge
  "Dihedral angle and convexity for a manifold edge.

  Convexity is decided by where the far corner of the second face sits relative
  to the first face's plane: outside it (along the normal) means the solid folds
  AWAY from the edge, which is a concave edge; inside means convex. Comparing
  normals alone cannot tell the two apart — a 90° convex corner and a 90° concave
  one give the same angle between normals."
  [pts faces edge f1 f2]
  (let [n1 (face-normal pts (nth faces f1))
        n2 (face-normal pts (nth faces f2))
        cosang (max -1.0 (min 1.0 (k/v-dot n1 n2)))
        angle (Math/acos cosang)
        far (nth pts (opposite-corner (nth faces f2) edge))
        on-plane (nth pts (first edge))
        side (k/v-dot n1 (k/v- far on-plane))]
    {:dihedral angle
     :convexity (cond (< angle 1.0e-9) :flat
                      (> side 1.0e-9) :concave
                      :else :convex)}))

(defn topology
  "Build the topology of `mesh` (a `{:positions :indices}` triangle soup).

  Returns
  `{:vertices :faces :edges :merged :degenerate :epsilon :manifold?}` where
  `:edges` maps `[i j]` (sorted vertex pair — the stable identity a soup lacks)
  to `{:faces [...] :kind :manifold|:boundary|:non-manifold :dihedral :convexity}`.
  `:dihedral` and `:convexity` are present only for manifold edges: on a
  boundary or non-manifold edge the question has no single answer, and
  answering it anyway is how a chamfer ends up cutting into a hole."
  ([mesh] (topology mesh (* 1.0e3 (get-in config/config [:epsilon :point-merge]))))
  ([mesh eps]
   (let [{:keys [positions indices merged degenerate]} (weld-mesh mesh eps)
         raw-faces (vec (map vec (partition 3 indices)))
         [faces orientable?] (if (seq raw-faces) (orient-faces positions raw-faces) [raw-faces true])
         edge->faces (reduce (fn [m [fi [a b c]]]
                               (reduce (fn [m e] (update m e (fnil conj []) fi))
                                       m
                                       [(edge-key a b) (edge-key b c) (edge-key c a)]))
                             {}
                             (map-indexed vector faces))
         edges (into {}
                     (map (fn [[e fs]]
                            (let [kind (edge-kind (count fs))]
                              [e (merge {:faces fs :kind kind}
                                        (when (= :manifold kind)
                                          (classify-edge positions faces e
                                                         (first fs) (second fs))))])))
                     edge->faces)]
     {:vertices positions
      :faces faces
      :edges edges
      :merged merged
      :degenerate degenerate
      :epsilon eps
      :orientable? orientable?
      :manifold? (every? #(= :manifold (:kind %)) (vals edges))})))

;; ---------------------------------------------------------------------------
;; Queries — what an edge feature actually asks
;; ---------------------------------------------------------------------------

(defn euler-characteristic
  "V - E + F. 2 for a closed surface of genus 0; a different value is the
  cheapest signal that a mesh is not the closed solid it is being treated as."
  [{:keys [vertices faces edges]}]
  (- (+ (count vertices) (count faces)) (count edges)))

(defn edges-where
  "Edge keys whose entry satisfies `pred`, in sorted order so a selection is
  reproducible."
  [{:keys [edges]} pred]
  (vec (sort (keep (fn [[e info]] (when (pred info) e)) edges))))

(defn sharp-edges
  "Manifold edges of the given `convexity` whose dihedral angle is at least
  `min-angle` radians. This is the selector an edge feature can actually use on
  a mesh: it names edges by a property of the geometry in front of it, not by a
  BREP id the mesh never carried."
  [topo convexity min-angle]
  (edges-where topo #(and (= :manifold (:kind %))
                          (= convexity (:convexity %))
                          (>= (:dihedral %) min-angle))))

(defn boundary-edges
  "Edges with a single adjacent face — the mesh has a hole there."
  [topo]
  (edges-where topo #(= :boundary (:kind %))))
