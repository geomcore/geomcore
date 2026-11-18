# Design decisions

Architecture decision records (ADRs) explain why `geomcore` is shaped the
way it is. They are historical: a later ADR may supersede an earlier one,
and the text says so where it happens. For the current API, use the
[reference](../geometry/elementary.md); for the vocabulary, use
[`CONTEXT.md`](https://github.com/geomcore/geomcore/blob/main/CONTEXT.md).

## Foundations

- [0002](0002-no-batch-transform-architecture.md) — No batch-transform architecture
- [0006](0006-angles-derive-from-the-confusion-tolerance.md) — Angles derive from the confusion tolerance and the geometry's extent
- [0009](0009-equality-is-exact-coincidence-takes-a-tolerance.md) — Equality is exact; coincidence takes a tolerance
- [0022](0022-each-public-item-has-one-path.md) — Each public item has one path, at the crate root
- [0027](0027-one-test-sample-list-per-shape-enum.md) — One test sample list per shape enum
- [0051](0051-a-result-does-not-depend-on-the-models-size-or-placement.md) — A result does not depend on the model's size or placement

## Types and the shape model

- [0005](0005-ranges-are-trimmed-curves.md) — Parameter ranges are given by trimmed curves
- [0008](0008-bspline-poles-stay-cartesian.md) — B-spline poles stay Cartesian, their weights apart
- [0010](0010-surface-normals-are-unit-point-variant-deferred.md) — Surface normals are unit, with second-derivative and step-off recovery
- [0013](0013-only-2d-types-carry-a-dimension-suffix.md) — Only 2D types carry a dimension suffix
- [0014](0014-the-2d-curve-tier-mirrors-the-3d-one.md) — The 2D curve tier mirrors the 3D one
- [0015](0015-shared-behaviour-lives-on-the-traits.md) — Shared behaviour lives on the traits
- [0016](0016-extrusion-and-revolution-surfaces-are-first-class.md) — Extrusion and revolution surfaces are first-class; offset is an operation
- [0017](0017-each-shape-enum-lists-its-variants-once.md) — Each shape enum lists its variants once
- [0024](0024-transform-and-reversal-are-shape-traits.md) — Transform and reversal are shape traits
- [0054](0054-a-unit-vector-is-a-direction-by-type.md) — A unit vector is a direction, by type

## Results and queries

- [0012](0012-results-describe-themselves.md) — Results describe themselves
- [0018](0018-exact-results-convert-straight-to-family-results.md) — Exact results convert straight to family results
- [0033](0033-distance-extrema-are-an-operation.md) — Distance extrema are an operation; projection belongs to the shape
- [0037](0037-evaluation-past-the-bounds-clamps-a-spline.md) — Evaluation past the bounds clamps a spline
- [0053](0053-a-result-pickles-with-its-operands-in-its-own-version.md) — A result pickles with its operands, in its own version

## Discretization, tessellation and fitting

- [0003](0003-discretization-deviations-are-proven.md) — Discretization deviations are proven bounds
- [0004](0004-discretization-fills-caller-owned-buffers.md) — Discretization fills caller-owned buffers
- [0056](0056-a-tessellation-is-a-parameter-grid-within-a-proven-deviation.md) — A tessellation is a parameter grid within a proven deviation
- [0036](0036-an-exact-quartic-is-fitted-within-confusion.md) — An exact quartic is fitted within confusion
- [0055](0055-a-fit-refuses-points-the-confusion-tolerance-cannot-tell-from-a-lower-shape.md) — A fit refuses points the confusion tolerance cannot tell from a lower shape

## Intersections

- [0001](0001-intersections-favour-speed-over-ergonomics.md) — Intersections favour speed over ergonomics (superseded in part by 0011 and 0018)
- [0007](0007-numeric-intersections-report-their-completeness.md) — Numeric intersections report their completeness
- [0011](0011-intersections-route-to-the-exact-solution-else-numeric.md) — Intersections route to the exact solution, else numeric, and say which ran
- [0019](0019-the-subdivision-searches-share-one-driver.md) — The subdivision searches share one driver
- [0023](0023-one-closed-form-pair-table-owns-routing.md) — One closed-form pair table owns routing
- [0025](0025-marching-is-a-surface-pair-a-tracer-and-a-seeding-policy.md) — Marching is a surface pair, a tracer and a seeding policy
- [0026](0026-any-surface-searches-as-patches-with-a-map-back.md) — Any surface searches as patches with a map back
- [0034](0034-a-tangent-contact-within-tolerance-is-complete.md) — A tangent contact within tolerance is one complete hit
- [0039](0039-a-trimmed-surface-clips-the-intersection-of-its-basis.md) — A trimmed surface clips the intersection of its basis
- [0040](0040-a-cone-and-a-cylinder-are-read-along-the-cylinder.md) — A cone and a cylinder are read along the cylinder
- [0041](0041-two-cones-are-read-along-the-first-cone.md) — Two cones are read along the first cone
- [0042](0042-an-unbounded-quartic-branch-is-an-exact-curve.md) — An unbounded quartic branch is an exact curve
- [0043](0043-touching-quadrics-cross-on-the-walk.md) — Touching quadrics cross on the walk
- [0044](0044-an-apex-on-the-other-surface-is-a-crossing-on-the-walk.md) — An apex on the other surface is a crossing on the walk
- [0045](0045-each-apex-on-the-other-cone-is-a-line-and-a-cubic.md) — Each apex on the other cone is a line and a cubic
- [0046](0046-a-march-seeds-around-a-singular-point-and-reports-a-seed-it-cannot-leave.md) — A march seeds around a singular point and reports a seed it cannot leave
- [0047](0047-circles-do-not-seed-a-march-against-a-torus.md) — Circles do not seed a march against a torus
- [0048](0048-tube-circles-share-the-quartic-walk.md) — Tube circles share the quartic walk
- [0049](0049-tube-circles-meet-a-cylinder-or-a-cone-in-up-to-four-points.md) — Tube circles meet a cylinder or a cone in up to four points
- [0050](0050-a-second-torus-is-read-along-the-same-tube-circles.md) — A second torus is read along the same tube circles
- [0052](0052-a-boundary-curve-lying-on-the-other-surface-is-one-whole-curve.md) — A boundary curve lying on the other surface is one whole curve

## Curves on surfaces

- [0030](0030-a-uv-curve-maps-to-an-exact-curve-or-is-refused.md) — A UV curve maps to an exact curve or is refused
- [0031](0031-the-surface-maps-a-3d-curve-to-its-uv-curve.md) — The surface maps a 3D curve to its UV curve, and curves do not depend on surfaces
- [0032](0032-a-uv-curve-on-a-spline-surface-is-composed-exactly.md) — A UV curve on a spline surface is composed exactly
- [0035](0035-a-uv-curve-maps-within-a-deviation-to-one-b-spline.md) — A UV curve maps within a deviation to one B-spline

## Python bindings

- [0021](0021-python-intersect-always-releases-the-interpreter-lock.md) — Python `intersect` always releases the interpreter lock
- [0029](0029-python-stubs-inherit-the-method-sets.md) — Python stubs inherit the method sets
- [0038](0038-the-python-package-is-flat-like-the-crate-root.md) — The Python package is flat, like the crate root
