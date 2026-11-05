# An exact quartic is fitted within confusion

Two cylinders whose axes are not parallel meet in a quartic, the bicylindric curve. So do a sphere and a cylinder whose axis misses its centre (Viviani's curve among them), and a sphere and a cone whose axis misses its centre. None of these is a conic, and in general no rational curve parametrizes them, so they have no `Curve`. Each does have a closed form along the generators of the ruled quadric, the cylinder or the cone (Dupont, Lazard, Lazard and Petitjean, *Near-optimal parameterization of the intersection of quadrics*, JSC 2008). The generator at angle `t`, `B(t) + s D(t)`, meets the other quadric where a quadratic in `s` vanishes. Its quarter discriminant `Δ` factors into terms `|cos τ − c|`, where `τ` is the angle from a reference generator:

- two cylinders: `Δ = A r₁² (upper − cos τ)(cos τ − lower)`, with `A = 1 − (a₁·a₂)²`, and `τ` measured from the common normal of the axes;
- a cylinder and a sphere: `Δ = 2 r e (cos τ − c)`, where `e` is the distance from the centre to the axis;
- a cone and a sphere: `Δ = (w·D)² − (|w|² − R²)`, where `w` runs from the apex to the centre and `w·D = p + q cos τ`. With the apex outside the sphere this factors into two terms.

The generators where `Δ ≥ 0`, and with them the outcome, follow from comparing one or two distances with the radii.

For two cylinders, the distance `d` between the axes decides:

- `d > r₁ + r₂`: empty;
- `d = r₁ + r₂` within confusion: the cylinders touch at one point;
- `d = |r₁ − r₂|` within confusion: the curve is nodal, and it goes to the numeric intersection;
- otherwise one loop or two, each closed. A loop runs around the first cylinder when every generator meets the second one; otherwise it runs over an arc of generators, there and back on the two roots.

For a sphere of radius `R` and a cylinder of radius `r` whose axis is `e` from its centre, `R` decides:

- `R > r + e`: two loops around the cylinder;
- `R = r + e`: Viviani's curve, a loop that crosses itself where the sphere touches the cylinder from inside;
- `|r − e| < R < r + e`: one loop over an arc of generators;
- `R = |r − e|`: a tangent point;
- `R < |r − e|`: empty.

For a sphere and a cone, the apex decides first. Inside the sphere, there is a loop on each nappe. On the sphere, there is one loop through the apex, the other root of each generator. Outside the sphere, each nappe that faces the centre is read alone, from the distances of the centre to its nearest and farthest generators in the plane of the axis and the centre: two loops, a loop that crosses itself, one loop over an arc, a tangent point, or nothing.

`intersect_exact` returns this closed form, a `QuarticCurve`, with its exact points. `intersect` fits each loop with a cubic Hermite B-spline through points of the closed form, together with the curves in both surfaces' parameters on the same knots. Each span is split until its midpoints and quarter points lie within a quarter of confusion of both surfaces, and until the curves in parameters map within that distance of the curve. The result is `Exact` and `Complete`. The point where the surfaces touch, where a loop crosses itself, or the apex on the sphere is in `contact_points`.

## Why `Exact`

`IntersectionMethod` says how the outcome was found, and `Completeness` says whether every piece was found. Both answers come from the closed form here. The number of loops, the point of contact, and the arc each loop covers are decided by comparing distances with the radii. Nothing is searched, so no loop can be missed and no tangent contact can be split. The fit only represents a curve whose every point is known exactly, as `circle_on_surface` represents a known circle. Its error, a quarter of confusion, is within the distance at which `Tolerance` already calls two points equal, so a consumer cannot tell it apart from an exact curve.

## Touching within confusion is touching

When the sphere's radius is within confusion of the radius at which it touches the ruled quadric, the sphere is read with that radius: `r + e` for a cylinder, or the distance from its centre to the generator it touches, or to the apex, for a cone. The curve then crosses itself exactly where the surfaces touch. Its points lie on the sphere of that radius, which is within confusion of the given one, and the fit is checked against it. The other choice, two loops a hair apart or one loop pinched almost shut, would be fitted through a square root that vanishes nowhere but changes by orders of magnitude within confusion. The tangent point of `Tolerance` is a point, so the crossing is one too.

## No curve in the sphere's parameters through a pole

A loop through a pole of the sphere has no continuous longitude: the angle jumps by `π` there. Viviani's curve in the sphere's own frame passes through both poles. There the curve in the sphere's parameters is `None`, as it is for a curve whose parameters cannot be fitted, and the ruled quadric's is kept. A loop that only comes near a pole is fitted, since its longitude turns fast but continuously.

## Considered options

- **Report it `Numeric`, as the marched curve was**: rejected. It would tell the caller that a search ran and might be incomplete. It would also send the caller to the marched result's checks, such as `IncompleteReason`, which never apply here.
- **A `QuarticCurve` variant in `Curve`**: rejected for now (YAGNI). Every consumer of `Curve` would need to handle a curve no file format stores, and a B-rep wants a B-spline edge anyway. `intersect_exact` gives the closed form to a caller that wants it.
- **One closed-form type per pair** (`BicylindricCurve`, a Viviani type, a cone type): rejected. The three share the walk along generators, the components and the fit; only the generator, the other quadric and the factors of `Δ` differ. One type keeps one fit and one set of tests for it.
- **A rational parametrization**: rejected. The quartic is rational only where it crosses itself. In general its parametrization needs the square root of `Δ`, which no B-spline reproduces exactly.
- **Fit to the caller's deviation**: rejected. `intersect` takes a `Tolerance` and no deviation, as the marcher does. Confusion is the deviation every intersection curve already promises.
- **Check only the midpoint of each span**: rejected. An error that is antisymmetric about the midpoint vanishes there. The midpoint and both quarter points cost two more evaluations per span and catch it.
