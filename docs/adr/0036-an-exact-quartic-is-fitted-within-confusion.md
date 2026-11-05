# An exact quartic is fitted within confusion

Two cylinders whose axes are not parallel meet in a quartic, the bicylindric curve. It is not a conic, and no rational curve parametrizes it in general, so it has no `Curve`. It does have a closed form along the generators of the first cylinder (Dupont, Lazard, Lazard and Petitjean, *Near-optimal parameterization of the intersection of quadrics*, JSC 2008). The generator at angle `t` meets the second cylinder at heights `z = (±√Δ(t) − β(t)) / A`. Here `A = 1 − (a₁·a₂)²` depends only on the axes, and `Δ` factors as `A r₁² (upper − cos τ)(cos τ − lower)`, with `τ` the angle from the common normal of the axes. Where `Δ ≥ 0` is read off `upper` and `lower`, which are the distance between the axes and the radii, so that distance alone decides the outcome:

- farther than `r₁ + r₂`: empty;
- `r₁ + r₂` within confusion: the cylinders touch at one point;
- `|r₁ − r₂|` within confusion: the curve is nodal, and it goes to the numeric intersection;
- otherwise one loop, or two, each closed. A loop runs around the first cylinder when every generator meets the second one; otherwise it runs over an arc of generators, there and back on the two roots.

`intersect_exact` returns this closed form, a `BicylindricCurve`, with its exact points. `intersect` fits each loop with a cubic Hermite B-spline, through points of the closed form, together with both curves in parameters on the same knots. Each span is split until its midpoints and quarter points lie within a quarter of confusion of both cylinders, and until the curves in parameters map within that distance of the curve. The result is `Exact` and `Complete`.

## Why `Exact`

`IntersectionMethod` says how the outcome was found, and `Completeness` says whether every piece was found. Both answers come from the closed form here. The number of loops, the point of contact, and the arc each loop covers are decided by comparing one distance with the radii. Nothing is searched, so no loop can be missed and no tangent contact can be split. The fit only represents a curve whose every point is known exactly, as `circle_on_surface` represents a known circle. Its error, a quarter of confusion, is within the distance at which `Tolerance` already calls two points equal, so a consumer cannot tell it apart from an exact curve.

## Considered options

- **Report it `Numeric`, as the marched curve was**: rejected. It would tell the caller that a search ran and might be incomplete. It would also send the caller to the marched result's checks, such as `IncompleteReason`, which never apply here.
- **A `BicylindricCurve` variant in `Curve`**: rejected for now (YAGNI). Every consumer of `Curve` would need to handle a curve no file format stores, and a B-rep wants a B-spline edge anyway. `intersect_exact` gives the closed form to a caller that wants it.
- **A rational parametrization**: rejected. The quartic of two cylinders is rational only in the nodal case. In general its parametrization needs the square root of `Δ`, which no B-spline reproduces exactly.
- **Fit to the caller's deviation**: rejected. `intersect` takes a `Tolerance` and no deviation, as the marcher does. Confusion is the deviation every intersection curve already promises.
- **Check only the midpoint of each span**: rejected. An error that is antisymmetric about the midpoint vanishes there. The midpoint and both quarter points cost two more evaluations per span and catch it.
