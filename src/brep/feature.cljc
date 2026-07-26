(ns brep.feature
  "Parametric feature tree: sketch, extrude, revolve, fillet, chamfer,
  boolean, etc. Restored from kami-cad's `feature` module (deleted
  PR #82). A `FeatureId` is a plain number (the original's newtype
  `FeatureId(u64)` wrapper adds no behavior). Sketch constraint kinds use
  the same keyword vocabulary as `kotoba-lang/engineer`'s
  `engineer.constraint/kinds` (:coincident/:parallel/etc.) — the original
  duplicated `ConstraintKind` as `SketchConstraintKind` for serde
  independence rather than a hard crate dependency; this mirrors that by
  using the same keywords without an explicit inter-repo dependency."
  (:require [brep.kernel :as k]
            [brep.tessellate :as tess]
            [brep.mesh-csg :as csg]))

(def boolean-ops #{:new :add :cut :intersect})

(def sketch-constraint-kinds
  "Mirrors engineer.constraint/kinds — see namespace docstring."
  #{:coincident :parallel :perpendicular :tangent :equal :horizontal :vertical
    :fixed :symmetric :concentric :midpoint :collinear :distance :angle
    :radius :diameter})

;; SketchPlane
(defn sketch-plane-xy [] {:kind :xy})
(defn sketch-plane-xz [] {:kind :xz})
(defn sketch-plane-yz [] {:kind :yz})
(defn sketch-plane-custom [origin normal] {:kind :custom :origin origin :normal normal})

(defn sketch-plane-normal [plane]
  (case (:kind plane)
    :xy [0.0 0.0 1.0]
    :xz [0.0 1.0 0.0]
    :yz [1.0 0.0 0.0]
    :custom (k/v-normalize (:normal plane))))

(defn sketch-plane-origin [plane]
  (case (:kind plane)
    (:xy :xz :yz) [0.0 0.0 0.0]
    :custom (:origin plane)))

;; SketchEntity variants
(defn sketch-line [start end] {:kind :line :start start :end end})
(defn sketch-arc [center radius start-angle end-angle]
  {:kind :arc :center center :radius radius :start-angle start-angle :end-angle end-angle})
(defn sketch-circle [center radius] {:kind :circle :center center :radius radius})
(defn sketch-spline [control-points] {:kind :spline :control-points control-points})
(defn sketch-dimension [entity-ref value] {:kind :dimension :entity-ref entity-ref :value value})
(defn sketch-constraint [kind entity-refs] {:kind :constraint :constraint-kind kind :entity-refs entity-refs})

;; Feature variants — each a map with :kind + :id + variant-specific fields
(defn sketch-feature [id plane entities] {:kind :sketch :id id :plane plane :entities entities})
(defn extrude-feature [id sketch-ref direction distance operation]
  {:kind :extrude :id id :sketch-ref sketch-ref :direction direction :distance distance :operation operation})
(defn revolve-feature [id sketch-ref axis angle operation]
  {:kind :revolve :id id :sketch-ref sketch-ref :axis axis :angle angle :operation operation})
(defn fillet-feature [id edges radius] {:kind :fillet :id id :edges edges :radius radius})
(defn chamfer-feature [id edges distance] {:kind :chamfer :id id :edges edges :distance distance})
(defn sweep-feature [id profile-ref path-ref operation]
  {:kind :sweep :id id :profile-ref profile-ref :path-ref path-ref :operation operation})
(defn loft-feature [id profiles operation] {:kind :loft :id id :profiles profiles :operation operation})
(defn shell-feature [id removed-faces thickness] {:kind :shell :id id :removed-faces removed-faces :thickness thickness})
(defn pattern-feature [id source-features direction count spacing]
  {:kind :pattern :id id :source-features source-features :direction direction :count count :spacing spacing})
(defn boolean-feature [id operation tool-body] {:kind :boolean :id id :operation operation :tool-body tool-body})

(defn feature-id [feature] (:id feature))

;; FeatureTree
(defn feature-tree
  "A fresh, empty parametric feature tree."
  []
  {:entries []})

(defn add-feature [tree feature]
  (update tree :entries conj {:feature feature :suppressed false}))

(defn suppress [tree id]
  (update tree :entries
          (fn [entries] (mapv (fn [e] (if (= (feature-id (:feature e)) id) (assoc e :suppressed true) e)) entries))))

(defn unsuppress [tree id]
  (update tree :entries
          (fn [entries] (mapv (fn [e] (if (= (feature-id (:feature e)) id) (assoc e :suppressed false) e)) entries))))

(defn reorder
  "Move the feature with `id` to `new-index`."
  [tree id new-index]
  (let [entries (:entries tree)
        pos (first (keep-indexed (fn [i e] (when (= (feature-id (:feature e)) id) i)) entries))]
    (if-not pos
      tree
      (let [entry (nth entries pos)
            without (vec (concat (subvec entries 0 pos) (subvec entries (inc pos))))
            idx (min new-index (count without))]
        (assoc tree :entries (vec (concat (subvec without 0 idx) [entry] (subvec without idx))))))))

(defn tree-len [tree] (count (:entries tree)))
(defn tree-empty? [tree] (empty? (:entries tree)))

(defn- sketch-by-id
  "Non-suppressed :sketch features in `entries`, keyed by :id — revolve
  looks its profile up by :sketch-ref, and (unlike extrude) actually
  needs it."
  [entries]
  (into {}
        (keep (fn [{:keys [feature suppressed]}]
                (when (and (not suppressed) (= :sketch (:kind feature)))
                  [(:id feature) feature])))
        entries))

(defn- revolve-profile-cylinder
  "The only revolve profile this evaluator supports today: `sketch`'s
  :entities contains exactly one sketch-line, and both its 2D [x y]
  endpoints share the same x (an axis-parallel profile — x is radius, y
  is axial offset from the sketch plane's origin). Returns
  {:radius :axial-min :axial-max}, or nil if the sketch doesn't match
  that shape (an angled line — a cone/frustum profile — or any other
  entity count/kind is not yet implemented; no cone tessellation exists
  in brep.tessellate yet to render one correctly)."
  [sketch]
  (let [lines (filter #(= :line (:kind %)) (:entities sketch))]
    (when (= 1 (count lines))
      (let [{:keys [start end]} (first lines)
            [x0 y0] start [x1 y1] end]
        (when (< (Math/abs (- x0 x1)) 1e-9)
          {:radius x0 :axial-min (min y0 y1) :axial-max (max y0 y1)})))))


;; ---------------------------------------------------------------------------
;; Sketch profile -> prism (2026-07-27).
;;
;; Before this, `evaluate`'s `:extrude :new` branch built a **unit-square**
;; prism and ignored the sketch entirely ("sketch entities are not yet
;; consumed here, only :direction/:distance"). So a feature tree describing a
;; 305 x 192 enclosure panel evaluated to a 1 x 1 box, and nothing in the
;; return value said so. The sketch was decorative.
;;
;; These helpers make the profile real for the case that is both common and
;; exactly representable in STEP: a single closed loop of straight segments,
;; extruded along the sketch-plane normal. The resulting prism has only planar
;; faces and linear edges, which is precisely the subset `brep.step/write-step`
;; accepts — so an extruded sketch now round-trips to STEP.
;; ---------------------------------------------------------------------------

(defn- pt2= [[ax ay] [bx by]]
  (and (< (Math/abs (- (double ax) (double bx))) 1e-7)
       (< (Math/abs (- (double ay) (double by))) 1e-7)))

(defn- circle->segments
  "Polygonise a sketch-circle into `n` straight segments. A circle is not
  exactly representable as a planar-faced prism, so this is an explicit
  approximation, not a silent one: callers asking for STEP get a prism whose
  facet count they chose."
  [{:keys [center radius]} n]
  (let [[cx cy] center]
    (mapv (fn [i]
            (let [t0 (* 2.0 Math/PI (/ (double i) n))
                  t1 (* 2.0 Math/PI (/ (double (inc i)) n))]
              {:kind :line
               :start [(+ cx (* radius (Math/cos t0))) (+ cy (* radius (Math/sin t0)))]
               :end   [(+ cx (* radius (Math/cos t1))) (+ cy (* radius (Math/sin t1)))]}))
          (range n))))

(defn sketch->ring
  "Chain a sketch's straight segments into one ordered closed ring of 2D
  points, or nil when the entities do not form exactly one closed loop.

  Returns nil (never a guess) when: there are no segments, an endpoint has
  no continuation, or the walk does not return to its start. `circle-facets`
  (default 32) controls circle polygonisation. Arcs and splines are not
  chained — they would need tessellation this function deliberately does not
  invent."
  ([sketch] (sketch->ring sketch 32))
  ([sketch circle-facets]
   (let [ents (:entities sketch)
         segs (vec (concat (filter #(= :line (:kind %)) ents)
                           (mapcat #(circle->segments % circle-facets)
                                   (filter #(= :circle (:kind %)) ents))))]
     (when (>= (count segs) 3)
       (loop [ring [(:start (first segs))]
              cur (:end (first segs))
              remaining (set (range 1 (count segs)))]
         (cond
           ;; closed the loop and used every segment
           (and (pt2= cur (first ring)) (empty? remaining)) ring
           ;; closed early, leaving segments over -> more than one loop
           (pt2= cur (first ring)) nil
           (empty? remaining) nil
           :else
           (if-let [i (first (filter (fn [i]
                                       (let [sg (nth segs i)]
                                         (or (pt2= cur (:start sg)) (pt2= cur (:end sg)))))
                                     remaining))]
             (let [sg (nth segs i)
                   nxt (if (pt2= cur (:start sg)) (:end sg) (:start sg))]
               (recur (conj ring cur) nxt (disj remaining i)))
             nil)))))))

(defn- to-3d
  "Lift a 2D sketch point onto the sketch plane."
  [plane [u v]]
  (case (:kind plane)
    :xy [(double u) (double v) 0.0]
    :xz [(double u) 0.0 (double v)]
    :yz [0.0 (double u) (double v)]
    :custom (let [n (k/v-normalize (:normal plane))
                  o (:origin plane)
                  ;; any vector not parallel to n, made orthonormal
                  a (if (< (Math/abs (nth n 0)) 0.9) [1.0 0.0 0.0] [0.0 1.0 0.0])
                  e1 (k/v-normalize (k/v-cross n a))
                  e2 (k/v-cross n e1)]
              (k/v+ o (k/v+ (k/v-scale e1 (double u)) (k/v-scale e2 (double v)))))))

(defn- newell-normal
  "Face normal from its own vertex loop (Newell's method) — robust to the
  ring's winding, unlike assuming CCW."
  [pts]
  (let [n (count pts)]
    (k/v-normalize
     (reduce (fn [acc i]
               (let [[x0 y0 z0] (nth pts i)
                     [x1 y1 z1] (nth pts (mod (inc i) n))]
                 (k/v+ acc [(* (- y0 y1) (+ z0 z1))
                            (* (- z0 z1) (+ x0 x1))
                            (* (- x0 x1) (+ y0 y1))])))
             [0.0 0.0 0.0] (range n)))))

(defn- centroid [pts]
  (k/v-scale (reduce k/v+ [0.0 0.0 0.0] pts) (/ 1.0 (count pts))))

(defn extrude-prism
  "Build a prism solid by extruding a closed 2D `ring` on `plane` along
  `dir` by `distance`. Returns `[solid edges vertices]` in `make-box`'s shape.

  All faces are `:plane` and all edges are `:line`, so the result is
  STEP-exportable via `brep.step/write-step`. Face normals are computed per
  face with Newell's method and then flipped to point away from the solid
  centroid, so the winding of the incoming ring does not matter."
  [id ring plane dir distance]
  (let [n (count ring)
        ext (k/v-scale (k/v-normalize dir) (double distance))
        base (mapv #(to-3d plane %) ring)
        top (mapv #(k/v+ % ext) base)
        all (into base top)
        verts (mapv (fn [i p] (k/brep-vertex (inc i) p)) (range) all)
        vid (fn [i] (:id (nth verts i)))
        pt (fn [i] (nth all i))
        mk-edge (fn [eid a b]
                  (k/brep-edge eid
                               (k/line-curve (pt a) (k/v-normalize (k/v- (pt b) (pt a))))
                               (vid a) (vid b)
                               [0.0 (k/v-distance (pt a) (pt b))]))
        base-edges (mapv (fn [i] (mk-edge (+ 100 i) i (mod (inc i) n))) (range n))
        top-edges (mapv (fn [i] (mk-edge (+ 200 i) (+ n i) (+ n (mod (inc i) n)))) (range n))
        vert-edges (mapv (fn [i] (mk-edge (+ 300 i) i (+ n i))) (range n))
        edges (into (into base-edges top-edges) vert-edges)
        solid-c (centroid all)
        mk-face (fn [fid eids face-pts]
                  (let [nrm (newell-normal face-pts)
                        c (centroid face-pts)
                        outward (if (neg? (k/v-dot nrm (k/v- c solid-c)))
                                  (k/v-scale nrm -1.0) nrm)]
                    (k/brep-face fid (k/plane-surface c outward)
                                 [(vec eids)] :forward)))
        bottom-face (mk-face 200 (map :id base-edges) base)
        top-face (mk-face 201 (map :id top-edges) top)
        side-faces (mapv (fn [i]
                           (let [j (mod (inc i) n)]
                             (mk-face (+ 210 i)
                                      [(:id (nth base-edges i)) (:id (nth vert-edges j))
                                       (:id (nth top-edges i)) (:id (nth vert-edges i))]
                                      [(pt i) (pt j) (pt (+ n j)) (pt (+ n i))])))
                         (range n))
        shell (k/brep-shell 300 (into [bottom-face top-face] side-faces) :forward)]
    [(k/brep-solid id [shell]) edges verts]))

(defn evaluate
  "Evaluate the feature tree, producing a BREP solid. Handles two base-
  feature cases:

  - `:extrude` with `:operation :new` generates a box-like prism from a
    unit-square cross-section along the extrusion direction (matches the
    original kami-cad Rust restoration — sketch entities are not yet
    consumed here, only :direction/:distance).
  - `:revolve` generates a real solid of revolution (brep.kernel/
    make-cylinder) ONLY for a full 2π turn of a single axis-parallel
    sketch-line profile (see revolve-profile-cylinder) — a cone/frustum
    profile (an angled line) or a partial-angle revolve returns
    `:error`, not a silently wrong shape; no cone tessellation exists
    yet to render either correctly.

  boolean/fillet/chamfer/sweep/loft/shell/pattern evaluation remain
  future work, tracked as TODOs in the source (a general boolean CSG
  kernel and arbitrary-profile revolve/sweep both need real geometry
  algorithms — polygon/polyhedron clipping, cone/freeform tessellation —
  that are correctness-critical enough not to rush).

  Returns `[:ok [solid edges verts]]` or `[:error msg]`."
  [tree]
  (let [sketches (sketch-by-id (:entries tree))]
    (loop [entries (:entries tree)
           result nil]
      (if (empty? entries)
        (if result [:ok result] [:error "feature tree produced no solid"])
        (let [{:keys [feature suppressed]} (first entries)]
          (if suppressed
            (recur (rest entries) result)
            (case (:kind feature)
              :extrude
              (case (:operation feature)
                :new
                (let [sk (get sketches (:sketch-ref feature))
                      ring (when sk (sketch->ring sk))]
                  (cond
                    ring
                    (recur (rest entries)
                           (extrude-prism 1 ring (:plane sk)
                                          (:direction feature) (:distance feature)))

                    ;; A sketch was named but is not one closed loop of
                    ;; segments -> refuse rather than fall back to a
                    ;; unit square that silently discards the profile.
                    sk
                    [:error (str "extrude: sketch " (:sketch-ref feature)
                                 " is not a single closed loop of straight"
                                 " segments (arcs/splines/multi-loop profiles"
                                 " are not implemented)")]

                    ;; No sketch reference at all: the documented legacy
                    ;; unit-square prism, kept for callers that relied on it.
                    :else
                    (let [half 0.5
                          ext (k/v-scale (k/v-normalize (:direction feature)) (:distance feature))
                          vmin [(- half) (- half) 0.0]
                          vmax (k/v+ [half half 0.0] ext)]
                      (recur (rest entries) (k/make-box 1 vmin vmax)))))

                (:add :cut :intersect)
                (if (nil? result)
                  [:error "no base solid to apply boolean to"]
                  ;; Previously this returned `result` unchanged — a boolean
                  ;; feature was a silent no-op that looked like success.
                  ;; Booleans need mesh CSG and produce a mesh, not a BREP,
                  ;; so they cannot be served by this BREP-returning fn.
                  [:error (str "evaluate: boolean operation " (:operation feature)
                               " on feature " (:id feature)
                               " cannot produce a BREP solid — use"
                               " brep.feature/evaluate-mesh, which applies it"
                               " via brep.mesh-csg and returns a triangle mesh")]))

              :revolve
              (let [full-turn? (>= (:angle feature) (- (* 2.0 Math/PI) 1e-6))
                    sketch (get sketches (:sketch-ref feature))
                    profile (when sketch (revolve-profile-cylinder sketch))]
                (cond
                  (not full-turn?)
                  [:error "revolve: only a full 2π turn is implemented (partial-angle pie-slice revolve is not yet supported)"]
                  (nil? sketch)
                  [:error (str "revolve: sketch-ref " (:sketch-ref feature) " not found")]
                  (nil? profile)
                  [:error "revolve: only a single axis-parallel sketch-line profile is implemented (an angled line -> cone/frustum, or any other profile shape, is not yet supported)"]
                  :else
                  (let [axis (k/v-normalize (:axis feature))
                        base (k/v+ (sketch-plane-origin (:plane sketch))
                                   (k/v-scale axis (:axial-min profile)))
                        height (- (:axial-max profile) (:axial-min profile))]
                    (recur (rest entries) (k/make-cylinder 1 base axis (:radius profile) height)))))

              :sketch (recur (rest entries) result)

              (recur (rest entries) result))))))))


;; ---------------------------------------------------------------------------
;; Boolean evaluation (2026-07-27).
;;
;; `evaluate` returns a BREP `[solid edges verts]`. A boolean of two BREPs is
;; not a BREP this kernel can build — there is no BREP-level clipper here —
;; but `brep.mesh-csg/mesh-boolean` does real CSG on triangle meshes. So
;; booleans get their own entry point with an honest return type: a mesh.
;;
;; Consequence to keep in view: a mesh cannot be written by
;; `brep.step/write-step` (which needs :plane faces and :line edges). A tree
;; containing cuts therefore yields a mesh for inspection/volume/rendering,
;; and the STEP path is only available for the pure-extrude subset. That is a
;; real limit of this kernel, stated rather than papered over.
;; ---------------------------------------------------------------------------

(defn- solid->mesh [[solid edges verts]]
  (let [[positions indices] (tess/tessellate-solid solid edges verts)]
    {:positions positions :indices indices}))

(defn evaluate-mesh
  "Evaluate a feature tree to a **triangle mesh**, applying boolean features
  via `brep.mesh-csg/mesh-boolean`.

  Base features are evaluated as in `evaluate` (extrude from a real sketch
  profile, revolve for the supported cylinder case) and then tessellated;
  `:add`/`:cut`/`:intersect` map to `:union`/`:difference`/`:intersection`.

  Each boolean's operand is the feature's own base geometry: an `:extrude`
  with `:operation :cut` extrudes its sketch, then subtracts that prism from
  the accumulated result. Returns `[:ok {:positions :indices}]` or
  `[:error msg]`.

  Not STEP-exportable — see the note above."
  [tree]
  (let [sketches (sketch-by-id (:entries tree))
        base-of (fn [feature]
                  (case (:kind feature)
                    :extrude
                    (let [sk (get sketches (:sketch-ref feature))
                          ring (when sk (sketch->ring sk))]
                      (if ring
                        [:ok (extrude-prism 1 ring (:plane sk)
                                            (:direction feature) (:distance feature))]
                        [:error (str "evaluate-mesh: extrude " (:id feature)
                                     " needs a sketch that is a single closed loop"
                                     " of straight segments")]))
                    [:error (str "evaluate-mesh: feature kind " (:kind feature)
                                 " is not supported as a boolean operand")]))]
    (loop [entries (:entries tree) acc nil]
      (if (empty? entries)
        (if acc [:ok acc] [:error "feature tree produced no geometry"])
        (let [{:keys [feature suppressed]} (first entries)]
          (cond
            suppressed (recur (rest entries) acc)
            (= :sketch (:kind feature)) (recur (rest entries) acc)
            :else
            (let [op (:operation feature)]
              (if (and (nil? acc) (contains? #{:add :cut :intersect} op))
                [:error "no base geometry to apply boolean to"]
                (let [[st base] (base-of feature)]
                  (if (= :error st)
                    [st base]
                    (let [m (solid->mesh base)]
                      (if (nil? acc)
                        (recur (rest entries) m)
                        (let [csg-op (case op :add :union :cut :difference
                                           :intersect :intersection :union)]
                          (recur (rest entries)
                                 (csg/mesh-boolean csg-op acc m)))))))))))))))
