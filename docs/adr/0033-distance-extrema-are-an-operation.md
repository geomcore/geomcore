# Distance extrema are an operation; projection belongs to the shape

The local minima of the distance between two curves lived inside the numeric intersection: the search in `intersect/numeric/curve_curve.rs`, its B-spline half in `curve_curve_hull.rs`, its results in `intersect/results/`. They now form their own operation, `operations/extrema/`, laid out like `discretization/` (ARCHITECTURE rule 3):

| File | Holds |
|---|---|
| `extrema/curve_curve.rs` | `curve_curve_extrema`, its seeds and refinement, `CurveCurveExtrema`, `CurveCurveExtremum` |
| `extrema/curve_curve_pieces.rs` | the minima search on the Bézier pieces of two B-splines |

**Point projection stays with each shape.** Projecting a point is asked for the point and its parameters; the distance comes with it. On every elementary shape it is the closed-form inversion and one evaluation: `parameters_of` then `eval_point` on a plane, cylinder, cone, sphere or torus, `parameter_of` then `eval_point` on a line or circle (the cone picks the nearer of its two generator lines in the meridian plane, the torus the nearer of its two tube circles, so the inversion of any point is its closest point). The inversion is core to the shape, used by intersections and UV curves, and `contains` is read from it. So `PointExtrema` (curves) and `project_point` (a required method of `ParametricSurface`) are implemented in the shape's file, as on surfaces; only their result types are shared, in `operations/projection.rs`. A first attempt moved the curve impls to `operations/`; it split each curve's inversion from its shape for no gain and left curves and surfaces organised differently, so it was undone.

**Extrema and intersection share the Bézier-piece machinery.** The numeric curve–curve intersection is read from the minima at zero distance (`search_curve_curve`, `hits_at_minima`), so it calls `extrema::curve_curve`. The other way round, the extrema search runs on the subdivision driver (ADR 0019), the B-spline pair layer (`pairs_with_bspline`) and the pair tests of `curve_curve_hull.rs`, which stay in `intersect/numeric/` because both searches use them, and its result carries the intersection's `Completeness` and `CurveCurveRegion`. The two modules depend on each other; moving that shared machinery to a neutral home is left until a third search needs it.

**Not done here.** Curve–surface and surface–surface extrema, and a local search refining one starting parameter next to the global one, have no caller yet. They go in `extrema/` when one appears (`curve_surface.rs`, `surface_surface.rs`).

The level-function extrema that seed marching are a different notion; their module is `marching/level_extrema.rs`, so that "extrema" in a path means distance extrema between shapes only.
