# Distance extrema are an operation; projection belongs to the shape

The local minima of the distance between two curves lived inside the numeric intersection: the search in `intersect/numeric/curve_curve.rs`, its B-spline half in `curve_curve_hull.rs`, its results in `intersect/results/`. They now form their own operation, `operations/extrema/`, laid out like `discretization/` (ARCHITECTURE rule 3):

| File | Holds |
|---|---|
| `extrema/curve_curve.rs` | `curve_curve_extrema`, its seeds and refinement, `CurveCurveExtrema`, `CurveCurveExtremum` |
| `extrema/curve_curve_pieces.rs` | the minima search on the Bézier pieces of two B-splines |

**Point projection stays with each shape.** Projecting a point is asked for the point and its parameters; the distance comes with it. On every elementary shape it is the closed-form inversion and one evaluation: `parameters_of` then `point_at` on a plane, cylinder, cone, sphere or torus, `parameter_of` then `point_at` on a line or circle (the cone picks the nearer of its two generator lines in the meridian plane, the torus the nearer of its two tube circles, so the inversion of any point is its closest point). The inversion is core to the shape, used by intersections and UV curves, and `contains` is read from it. So `PointExtrema` (curves) and `project_point` (a required method of `ParametricSurface`) are implemented in the shape's file, as on surfaces; only their result types are shared, in `operations/projection.rs`. A first attempt moved the curve impls to `operations/`; it split each curve's inversion from its shape for no gain and left curves and surfaces organised differently, so it was undone.

**Extrema calls the shared Bézier-piece machinery, never the reverse.** The numeric curve–curve intersection is read from the minima at zero distance (`search_curve_curve`, `hits_at_minima`), so it calls `extrema::curve_curve`, and its own grid fallback is `sampled_curve_curve_extrema`, a named entry point rather than a flag on the search. The extrema search runs on machinery the intersection searches share and that stays in `intersect/numeric/`: the subdivision driver (ADR 0019), the Bézier pieces, a curve as a B-spline over a range with the map back (`pairs_with_bspline`), and the curve's parameter domain (`parameter_domain.rs`; its seeding for the minima search stays in extrema). The seeding of the curve-pair minima (`curve_curve_seeds`) lives in extrema, so that machinery never calls extrema. Two result types still cross: the intersection's meetings are reported as `CurveCurveExtremum`s, and the extrema carry the intersection's `Completeness` and `CurveCurveRegion`.

That machinery serves four searches (curve–curve hits and minima, curve–surface hits, the patch pairs that seed marching), not only the intersection. Moving it to a home of its own is left until curve–surface extrema need it too. The release build uses the default codegen units, so module layout changes what is inlined, and a split once cost 2.3 % on a benchmark: such a move is benchmarked first.

**Not done here.** Curve–surface and surface–surface extrema, and a local search refining one starting parameter next to the global one, have no caller yet. They go in `extrema/` when one appears (`curve_surface.rs`, `surface_surface.rs`).

The level-function extrema that seed marching are a different notion; their module is `marching/level_extrema.rs`, so that "extrema" in a path means distance extrema between shapes only.

## Amendment: one Python method for both extrema (#201)

In Python, a 3D curve's `extrema(other)` takes a point or a curve. A point gives its projections, a list of `CurveProjection`. A curve gives the `CurveCurveExtrema` of the pair, which `geomcore.curve_curve_extrema(first, second)` returned before; that function now warns and goes in the next release. Users tried `curve.extrema(other_curve)` first and got a type error about `Point`. The Rust split stays: `PointExtrema::extrema` with each shape, and the free function `curve_curve_extrema` in `operations/extrema/`. A 2D curve's `extrema` still takes a point alone, since there is no 2D curve-pair extrema.

The point case keeps its plain list. A list has `len` and iteration, and a point's projections have no completeness to report, so a result class would only wrap them. Point extrema on the elementary and swept surfaces, which only B-spline and Bézier surfaces have today, are a separate issue.
