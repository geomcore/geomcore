# Intersection

One entry point, two methods. `a.intersect(&b, tol)` works for any two
curves or surfaces, B-splines, Bézier types and trimmed curves
included, in either order. It solves the pair exactly whenever a closed form exists, and
runs a numeric search (curve–curve, curve–surface, marching
surface–surface) otherwise. The result says which method ran and
whether every hit was found, so a switch from nanoseconds to
milliseconds is never hidden. Adapted from
[Intersection](https://avenirgeometry.com/articles/intersection/).

## What an intersection returns

One result per family: `CurveCurveIntersection` (and
`CurveCurveIntersection2D`, the plane's instance of the same
`GenericCurveCurveIntersection`), `CurveSurfaceIntersection` and
`SurfaceSurfaceIntersection`. Each one carries:

- the **hits**: the point, the multiplicity (1 for a crossing, 2+ for a
  tangency) and the parameters on both operands, always filled: where a
  closed form solves for one operand's parameter only, the other's is
  recovered in closed form from the point. Hits are sorted by the first
  operand's parameter;
- the **overlaps**, where a curve lies on the other shape: a parameter
  range on each curve plus `same_direction`, or the curve's ranges on a
  surface;
- for two surfaces, the **curves** (each with a `tangent` flag and its
  curve in the parameters of each surface, see below), the **contact
  points** and whether the surfaces are **coincident**;
- the **method** (`Exact` or `Numeric`) and the **completeness**: an
  exact result is always complete; a numeric one is complete, incomplete
  (with the regions where hits may be missing) or sampled. The method
  says how the result was found, not what type its curves are: an exact
  `Circle` can come out of a numeric search, as on the swept surfaces
  below.

Hits on circles and ellipses keep their parameters inside `bounds()`,
and a tangency is decided by distance, not by root separation:

<figure class="geom-fig">
<img src="../../assets/figures/intersection/01-a-tangency-is-a-distance.svg" alt="Tangency decided by distance">
<figcaption>A tangency is a distance: a closest approach within confusion is one tangent hit, not two hits or none.</figcaption>
</figure>

## Exact: one call, either order

=== "Python"

    ```python
    result = plane.intersect(sphere)    # SurfaceSurfaceIntersection
    circle = result.curves[0].curve     # Circle
    hits = line.intersect(sphere).hits  # [CurveSurfaceHit], in either order
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{IntersectionMethod, Tolerance};

    let tol = Tolerance::DEFAULT;
    let result = plane.intersect(&sphere, tol);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert_eq!(sphere.intersect(&plane, tol), result);
    ```

A pair of concrete types calls its exact solver directly, with no
run-time dispatch; only `Curve`, `Surface`, `ElementaryCurve` or
`ElementarySurface` operands pay for a `match`. For hot loops over a
pair of concrete types, `a.intersect_exact(&b, tol)` (the `IntersectExact`
trait) returns the pair's own enum, such as `PlaneSphereIntersection`,
on the stack, recovers no parameter it does not need, and reports
`NotAnalytic` where no closed form exists instead of searching.

Covered pairs include plane/quadric/torus surfaces against each other
in their symmetric configurations — among them the tee of two equal
cylinders whose axes meet, two ellipses that cross at the two points
where the cylinders touch, both listed in `contact_points` — and any two
cylinders whose axes cross, a sphere and a cylinder or a cone whose
axis misses its centre, a cone and a cylinder off its axis, or two
cones in general position. Those meet in a quartic, in one loop or two:
the bicylindric curve of two cylinders, or, for a sphere, Viviani's
curve among others. Where a generator of one cone runs along one of the
other's, the curve runs off to infinity there: each piece between two
such generators comes back whole as an exact `QuarticBranch`, unbounded
both ways, the loops that stay bounded fitted as usual. A cylinder whose
axis runs along a generator of a cone meets it once on each of its own
generators, on one loop that is bounded, or runs off to infinity in
branches like those. A cone and a cylinder, or two cones, that touch
at a point meet in a curve crossing itself there, one loop through the
point twice, or at that point alone; so does a cone whose apex lies on
the other surface, at its apex. Where each apex of two cones lies on
the other, the line through both lies on both, and the rest is a cubic
through both apexes: `intersect` returns the line and the cubic's
pieces between its points at infinity as `QuarticBranch` curves, with
both apexes as contact points. So does a cylinder along a generator of
a cone through its apex: that generator and a cubic through the apex,
in one branch, the apex the contact point. Surfaces that touch twice
and meet nowhere else, such as a pipe beside a cone's apex touching
both nappes, give those two points as `contact_points`.
`intersect_exact` returns its closed form
(`QuarticCurve`, exact points on both surfaces), and `intersect` fits
each loop within confusion and reports the result exact and complete,
with the point where the surfaces touch, if any, in `contact_points`:
two cylinders with axes `|r1 − r2|` apart meet in one loop that crosses
itself there. Two cones with one apex and axes apart split into lines
instead: up to four through the apex, the generators they share, with the apex in
`contact_points`. A cone and a cylinder, or two cones, whose axes meet
at a point as far from both, so that both touch one sphere centred
there, split into two plane conics (the ellipses of a reducer), with the
points where they touch in `contact_points`. Two cones with parallel
axes and equal half-angles, each the other moved sideways, meet in one
plane conic: the hyperbola of two neighbouring chamfers.
A torus meets a plane or a sphere in any position in closed form too,
read along its tube circles instead of a quadric's generators. A plane
cuts a spiric section: two loops around the axis, one loop or two
(Cassini's ovals, parallel to the axis at the minor radius from it), or
one crossing itself where the plane touches the inside of the tube, that
point in `contact_points`. A sphere off the axis meets it in a quartic of the
same shapes. A sphere touching it twice meets it in two circles through
both points, and one centred in the torus's plane can hold two tube
circles. A plane or a sphere that touches the torus and meets it nowhere
else gives that point. A horn or spindle torus, whose tube circles cross
on its axis, runs the numeric search away from the special positions.
A cylinder or a cone off the torus's axis meets it in a curve of degree
eight, an `OcticCurve`, read along the same tube circles: each holds up
to four of its points, the roots of a quartic. A pipe through the tube
cuts two loops, one along a diameter four, and one touching the torus
crosses itself there, that point in `contact_points`, as does a cone
whose apex lies on the torus. One touching it and meeting it nowhere
else gives that point. A horn or spindle torus, or a cylinder or a cone
holding a tube circle, runs the numeric search.
Two tori off one axis meet in an `OcticCurve` too, read along the first
one's tube circles: two loops for rings linked like a chain, loops
crossing where tori side by side touch, a point where they touch and
part. Pipe bends joined end to end share a tube circle and touch all
along it: it comes back as a tangent curve (`OcticCurve::shared_circle`),
with the loops of the rest of the curve. A circle of one torus lying on
the other otherwise, a latitude or Villarceau circle, runs the numeric
search.
Lines/circles/conics against each
other and against the elementary surfaces, except an ellipse, a
parabola or a hyperbola against a torus: substituted into the torus they
give an equation of degree 8, and that pair runs the numeric search. A
circle against a torus gives one of degree 4, the torus being
bicircular: up to four hits, or the whole circle where it lies on the
torus (a tube, latitude or Villarceau circle). Any two conics —
coplanar (up to four hits) or in crossing planes (met on the shared
line). The 2D pairs are exact too: any two of `Line2D`, `Circle2D`,
`Ellipse2D`, `Parabola2D` and `Hyperbola2D`, whole or trimmed, solved
by the coplanar solvers with the same tangency multiplicities and
overlaps; only a pair holding a `BSplineCurve2D` or a `BezierCurve2D`
runs the numeric search. A Bézier curve or surface is searched as the
B-spline it is.

A plane meets an extrusion or a revolution in closed form in four
positions: generator lines where it holds the extrusion's direction, the
profile carried onto it otherwise, circles where it is perpendicular to
the revolution's axis, and the profile turned into it where it holds the
axis. The curves are exact types, a `Line`, a `Circle` or a curve of the
profile's kind. They sit where the profile crosses the plane, though,
and that crossing is exact only for an elementary profile. On a
B-spline or Bézier profile it comes from the numeric curve–plane search,
so a horizontal cut of a spline revolution is a list of `Circle`s with
`method` `Numeric` and completeness `Complete`. A sphere cut by a coaxial cylinder, computed live below —
two exact rings in ~0.4 µs on a Linux x86-64 desktop, single thread:

<iframe src="https://geomcore.github.io/geomcore-viewer/?scene=https://geomcore.github.io/geomcore/assets/scenes/intersect-exact.json" width="100%" height="480" style="border: 1px solid #333; border-radius: 8px;" allow="fullscreen" title="Sphere and cylinder with their two exact intersection rings"></iframe>

Two more landmarks, as figures:

<figure class="geom-fig">
<img src="../../assets/figures/intersection/15-plane-and-plane.svg" alt="Plane and plane intersection">
<figcaption>Plane against plane: a line, or Parallel / Coincident classification by tolerance.</figcaption>
</figure>

<figure class="geom-fig">
<img src="../../assets/figures/intersection/16-plane-and-sphere.svg" alt="Plane and sphere intersection">
<figcaption>Plane against sphere: a circle, a tangent point, or empty — classified, not sampled.</figcaption>
</figure>

<iframe src="https://geomcore.github.io/geomcore-viewer/?scene=https://geomcore.github.io/geomcore/assets/scenes/plane-sphere.json" width="100%" height="480" style="border: 1px solid #333; border-radius: 8px;" allow="fullscreen" title="A plane against a sphere at three distances: a circle, a tangent point, nothing"></iframe>

## Numeric: curves

The routed `a.intersect(&b, tol)` takes the numeric search for pairs
without a closed form (a B-spline, or an ellipse against a torus), and
trimmed pieces are searched over their own range. A curve or surface
type you define yourself, by implementing `ParametricCurve` or
`ParametricSurface`, goes through `intersect_curve_curve(a, b, tol)` and
`intersect_curve_surface(curve, surface, tol)`: they always run the
numeric search and return the same result types as the routed call.
On B-splines there is no sampling grid: Bézier pieces subdivide where their control boxes overlap until each
pair holds at most one hit, then Newton solves it — crossings a grid
would miss between samples are found, and tangencies/overlaps fall
back to the sampled search instead of burning the budget.

Where a curve lies on the other operand (a planar B-spline in a plane, a
straight B-spline on a line, two B-splines along one line, a circle and
its own `to_bspline()`), the pieces never separate. The hits and regions
left along the curve are then joined into an **overlap** wherever the
curve stays on the other operand between them (samples within
`tol.confusion`, tangents along it), extended by halving to where the
curve leaves it, within `tol.confusion` in space of the true end. The
hits inside drop out, and the result stays `Incomplete`, with an
`Unseparated` region over the overlap: samples found it, not a proof. A
stretch where the operands only touch (they part as the square of the
distance from the contact) is one hit at the closest approach, with
multiplicity 2, and the result stays `Complete`; an overlap between curves whose
range on the second curve would cross that curve's seam is not joined,
and keeps its hits and region.

The search reports whether it finished: the result's `method()` is
`Numeric`, and its `completeness()` is `Complete` when every piece pair
was proven apart or isolated, `Incomplete` with the parameter regions
where hits may be missing otherwise, and `Sampled` when the hits come
from samples (a pair without a B-spline, such as an ellipse against a
torus). Between two curves the sampled search halves its grid's spans
wherever their chords come within the curves' offset from them, so
crossings closer together than the grid spacing are each found.

=== "Python"

    ```python
    result = spline_curve.intersect(spline_surface)  # CurveSurfaceIntersection
    for hit in result.hits:
        hit.curve_parameter, hit.u, hit.v, hit.point
    complete = result.is_complete()
    intersection = curve_a.intersect(curve_b)        # CurveCurveIntersection
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;

    let result = spline.intersect(&surface, tol);
    for h in result.hits() {
        let _ = (h.curve_parameter, h.surface_parameters, h.point);
    }
    let complete = result.is_complete();
    let intersection = a.intersect(&b, tol);
    ```

## Numeric: marching surfaces

For a pair without a closed form, the routed `intersect` of two
surfaces seeds from the points where iso-curves of an elementary side
(rulings and circles of a cylinder or cone, circles of a sphere or
torus, lines of a plane) meet the other surface exactly, so an
unbounded surface is searched whole. A B-spline side adds its own
seeds, and two bounded surfaces add a parameter grid. It then traces
with curvature-following steps corrected by Newton on both surfaces,
until the trace closes or leaves its box: the smaller bounded side's
box, or, when neither side is bounded, the box of the seeds grown by
the larger surface size. A branch that leaves that seeded box is cut
there and reported as unbounded. The trace fits each
loop or arc as a `BSplineCurve`: where the surfaces cross, cubic
Hermite segments through the traced points and their exact tangents,
each checked at its midpoint, so the whole curve stays within
confusion of both surfaces. Two B-spline surfaces additionally
seed from overlapping patch pairs, so isolated loops smaller than the
seed spacing are found there. A marched result is `Sampled`: a loop
smaller than the iso-curve spacing can be missed. Where a branch of two
unbounded surfaces was cut at its box, the result is `Incomplete` with
`Unbounded` as its reason, and its regions mark where the cuts are.

The confusion also sets the size of a fitted curve: each segment must
stay within it of both surfaces, so a looser confusion takes fewer,
longer segments. Two wavy B-spline patches cut at the default `1e-7`
give curves of 296 and 194 poles; at `1e-3` they give one curve of 42
poles. The same value decides what counts as a hit or a tangency, so
loosen it for a lighter fit only where a point within it of both
surfaces is on both. In Python, `tol=1e-3` is short for
`Tolerance(confusion=1e-3)`.

<figure class="geom-fig">
<img src="../../assets/figures/intersection/31-marching-round-a-cylinder-and-a-sphere.svg" alt="Marching round a cylinder and a sphere">
<figcaption>Marching: seed, trace along the tangent, correct back onto both surfaces, fit the trail.</figcaption>
</figure>

<figure class="geom-fig">
<img src="../../assets/figures/intersection/33-the-loop-no-boundary-sees.svg" alt="An isolated loop">
<figcaption>The hard case: an isolated loop no seed grid line crosses. Patch-pair seeding finds these on B-spline pairs.</figcaption>
</figure>

=== "Python"

    ```python
    # Exact curves (Line, Circle, ...) where a closed form exists,
    # BSplineCurve where marched.
    result = surface_a.intersect(surface_b)
    curves = [piece.curve for piece in result.curves]
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;

    let result = a.intersect(&b, tol);
    let curves = result.curves();   // BSpline curves where marched
    ```

### The curves in each surface's parameters

Each curve also comes as `first_uv_curve` and `second_uv_curve`: the
curve in the `(u, v)` parameters of each surface, in the same parameter
as `curve`, so the surface's point at the `(u, v)` of `t` is the curve's
point at `t`. A marched curve always has both: the march walked the
parameters of every point it traced, and they are fitted on the curve's
own knots, continuous across the seams (a loop around a cylinder ends a
period away in `u`). An exact curve has the exact one where a closed
form exists (a circle in a plane's parameters, a ruling or a parallel
circle as a line in a cylinder's) and `None` otherwise, as an ellipse
on a cylinder: [`curve_3d_to_uv_within`](parametrization.md) maps those
within a deviation.

=== "Python"

    ```python
    from geomcore import Cylinder, Plane, Point, Vector

    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 1.0)
    (piece,) = Plane(Point(0.0, 0.0, 2.0), Vector.unit_z()).intersect(cylinder).curves
    on_plane = piece.first_uv_curve       # Circle2D
    on_cylinder = piece.second_uv_curve   # Line2D: the iso-line v = 2
    uv = on_cylinder.point_at(0.7)
    assert cylinder.point_at(uv.x, uv.y).distance(piece.curve.point_at(0.7)) < 1e-12
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{Cylinder, Plane, Point, Tolerance, Vector};

    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let plane = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
    let result = plane.intersect(&cylinder, Tolerance::DEFAULT);
    let piece = &result.curves()[0];
    let on_cylinder = piece.second_uv_curve.as_ref().unwrap(); // the line v = 2
    ```

Before seeding, a surface with bounded parameters is checked for lying
on the other: a 7 × 7 grid of its points, and the centres of the
grid's cells, must lie within confusion of the other surface, with the
normals parallel within the tilt that distance allows between samples
(about twice the confusion over their spacing), so a fitted or slightly
moved copy still counts. Surfaces lying on each
other, such as a sphere and its own `to_bspline()`, come back coincident
(`Sampled`, since samples found it). A patch lying on a larger or
unbounded surface, such as a cylinder and its `to_bspline_patch`, comes
back `Incomplete` with one `Unseparated` region covering the patch, and
no curves. Where the surfaces touch along a curve (a plane tangent to a
cylinder patch along a ruling, a sphere inside a cylinder of its radius)
the curve is traced along the direction in which they keep touching and
reported with `tangent` set; a touch at a point is a contact point,
at a sphere's pole too.

Marching is approximate by design. Its limits:
- the curves are fitted, to about 1e-4 at the default tolerance;
- at a point where a surface has no tangent plane (a cone's apex, a
  spindle torus's points on its axis), no trace can step: the branches
  leaving it are seeded a first step away, where the other surface
  crosses the ring of points around it. A branch along which the
  surfaces only touch is kept where the other surface comes within
  confusion of that ring. A seed that no trace leaves and no curve
  passes through makes the result `Incomplete` with an `Unseparated`
  region at it, rather than dropping it;
- the coincidence check sees only a whole surface lying on the other:
  two patches that overlap in part are marched; the seeds where they lie
  on each other make the result `Incomplete` with an `Unseparated`
  region around the overlap, but spurious curves may still cross it;
- where crossing curves touch each other, traces stop rather than turn
  onto the other branch;
- a loop that no seeding iso-curve crosses and that leaves no distance
  minimum on either grid is missed: a plane cutting a small circle
  (radius 0.014) off a B-spline sphere of radius 1 near its equator or
  pole comes back `Sampled` with no curves.

A wavy B-spline patch through a sphere, marched live below — one fitted
loop in ~8 ms on the same machine (the loop arrives in the viewer as
a sampled polyline until native B-spline types land,
[viewer#5](https://github.com/geomcore/geomcore-viewer/issues/5);
tunable sampling density is
[viewer#12](https://github.com/geomcore/geomcore-viewer/issues/12)):

<iframe src="https://geomcore.github.io/geomcore-viewer/?scene=https://geomcore.github.io/geomcore/assets/scenes/intersect-marching.json" width="100%" height="480" style="border: 1px solid #333; border-radius: 8px;" allow="fullscreen" title="B-spline patch and sphere with their marching intersection loop"></iframe>

## Trimmed surfaces

A face lies on a patch of its surface, not the whole of it.
`surface.trimmed(u_range, v_range)` restricts a surface to a rectangle
of its parameters, a `TrimmedSurface`; a periodic direction's range may
cross the seam, as `(1.5π, 2.5π)` on a cylinder. `intersect` takes it
with any surface or another trimmed surface, in either order: the bases
are intersected, exactly where a closed form exists, and each curve is
cut to the pieces inside the ranges of the trimmed sides. A piece keeps
the curve of the bases, and its `parameter_range` says which part lies
on both patches; an untrimmed pair's is the curve's own bounds. A
curve cut on both sides of its seam keeps the piece across it whole,
its range running past the period.

=== "Python"

    ```python
    import math
    from geomcore import Cylinder, Plane, Point, Vector

    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 1.0)
    half = cylinder.trimmed((0.0, math.pi), (-5.0, 5.0))
    plane = Plane(Point.origin(), Vector(0.0, 1.0, 1.0))
    (piece,) = plane.intersect(half).curves   # half of the ellipse
    start, end = piece.parameter_range
    assert math.isclose(end - start, math.pi)
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{Cylinder, Plane, Point, Tolerance, Vector};
    use std::f64::consts::PI;

    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let half = cylinder.trimmed((0.0, PI), (-5.0, 5.0)).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0)).unwrap();
    let result = plane.intersect(&half, Tolerance::DEFAULT);
    let (start, end) = result.curves()[0].parameter_range;
    assert!((end - start - PI).abs() < 1e-9);
    ```

A march between trimmed surfaces is held to the box of the patches, so
two patches of unbounded surfaces whose intersection runs to infinity
and has no closed form come back `Sampled` with curves ending on the
ranges instead of cut and `Incomplete`. An exact curve that runs to
infinity, a line or a `QuarticBranch`, is cut where it leaves the
ranges, and the result stays exact. Contact points outside the ranges are dropped, and two
patches of one surface are coincident only where they overlap. A curve
and a trimmed surface do not intersect yet: intersect the curve with
the `basis()`.

All the scenes on this page were computed by the kernel itself —
see `docs/assets/scenes/gen_scenes.py` — not drawn by hand.
