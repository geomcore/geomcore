# Intersection

Call `a.intersect(b)` in Python or `a.intersect(&b, tol)` in Rust.
Either operand order works for the supported pairs. Read `method` to
see how the kernel found the result and `completeness` to see whether the
search may have missed intersections.

| Operands | Support |
|---|---|
| Two curves | Use two 2D curves or two 3D curves. Trimmed curves are supported. |
| Curve and surface | Use a 3D curve and a surface. A trimmed surface clips the result to its parameter ranges. |
| Two surfaces | Both operands can be trimmed. |

Adapted from [Intersection](https://avenirgeometry.com/articles/intersection/).

## Exact: one call, either order

The Python examples on this page share these shapes:

```python
from geomcore import BSplineCurve, BSplineSurface, Line, Plane, Point, Sphere, Vector

plane = Plane(Point(0.0, 0.0, 0.0), Vector.unit_z())
sphere = Sphere(Point(0.0, 0.0, 0.5), 1.0)
line = Line(Point(-2.0, 0.0, 0.5), Vector.unit_x())

# A dome over the unit square, a curve diving through it, and a second
# curve through the first one's middle.
dome = [
    [Point(x, y, 1.0 if x == y == 0.5 else 0.0) for y in (0.0, 0.5, 1.0)]
    for x in (0.0, 0.5, 1.0)
]
spline_surface = BSplineSurface(2, 2, dome, [0.0, 1.0], [3, 3], [0.0, 1.0], [3, 3], False, False)
spline_curve = BSplineCurve.interpolate(
    [Point(0.1, 0.2, 1.0), Point(0.4, 0.5, -0.5), Point(0.7, 0.4, 1.0), Point(0.9, 0.8, -0.5)], 3
)
curve_a = spline_curve
curve_b = BSplineCurve.interpolate(
    [Point(0.0, 0.5, 0.0), curve_a.point_at(0.5), Point(1.0, 0.0, 0.5)], 2
)
surface_a = spline_surface
surface_b = Plane(Point(0.0, 0.0, 0.25), Vector(0.0, 0.2, 1.0))
```

=== "Python"

    ```python
    result = plane.intersect(sphere)    # SurfaceSurfaceIntersection
    assert result.is_complete()
    circle = result.curves[0].curve     # Circle
    hits = line.intersect(sphere).hits  # [CurveSurfaceHit], in either order
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{IntersectionMethod, Plane, Point, Sphere, Tolerance, Vector};

    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let sphere = Sphere::new(Point::new(0.0, 0.0, 0.5), 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    let result = plane.intersect(&sphere, tol);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    assert_eq!(result.curves().len(), 1);
    ```

## What an intersection returns

The result type depends on the operands. Curve results give isolated
hits with coordinates, parameters and multiplicity. A simple crossing
has multiplicity 1. A tangency has multiplicity 2 or higher. Overlaps
carry parameter ranges. Surface–surface results give curves, contact
points and coincidence information. Each intersection curve also has a
tangency flag and its UV curve on each surface where available.

Read these fields before using the result:

| Field | Meaning |
|---|---|
| `method` | The route used to solve the pair, exact or numeric. An exact route can fit a closed-form loop as a B-spline. A numeric route can return a circle. |
| `completeness` | Complete, incomplete with unresolved regions, or sampled. Use `is_complete()` before relying on the search to have found every hit. |
| `deviation` | Per surface–surface intersection curve. Zero for a closed-form curve. For a fitted curve, the reported distance from the surfaces. |

Hits on circles and ellipses keep their parameters inside `bounds()`.
The kernel uses distance to classify tangencies within tolerance.

<figure class="geom-fig">
<img src="../../assets/figures/intersection/01-a-tangency-is-a-distance.svg" alt="Tangency decided by distance">
<figcaption>A tangency is a distance: a closest approach within confusion is one tangent hit, not two hits or none.</figcaption>
</figure>

## Pair-specific exact solvers

In Rust, a pair of concrete types calls its exact solver directly, with
no run-time dispatch; only `Curve`, `Surface`, `ElementaryCurve` or
`ElementarySurface` operands pay for a `match`. For hot loops over a
pair of concrete types, `a.intersect_exact(&b, tol)` (the `IntersectExact`
trait, Rust only) returns the pair's own enum, such as `PlaneSphereIntersection`,
on the stack, recovers no parameter it does not need, and reports
`NotAnalytic` where no closed form exists instead of searching.

See the [closed-form pair reference](../geometry/intersection-pairs.md)
for quadric and torus sections, conics, and plane cuts of swept surfaces.

A coaxial cylinder cuts a sphere in two rings.

```viewer
scene: assets/scenes/intersect-exact.json
title: Sphere and cylinder with their two exact intersection rings
```

Two more landmarks, as figures:

<figure class="geom-fig">
<img src="../../assets/figures/intersection/15-plane-and-plane.svg" alt="Plane and plane intersection">
<figcaption>Plane against plane: a line, or Parallel / Coincident classification by tolerance.</figcaption>
</figure>

<figure class="geom-fig">
<img src="../../assets/figures/intersection/16-plane-and-sphere.svg" alt="Plane and sphere intersection">
<figcaption>Plane against sphere: a circle, a tangent point, or empty — classified, not sampled.</figcaption>
</figure>

```viewer
scene: assets/scenes/plane-sphere.json
title: A plane against a sphere at three distances: a circle, a tangent point, nothing
```

## Numeric: curves

The routed `intersect` takes the numeric search for pairs
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
`Unseparated` region over the overlap: samples found it, not a proof.
The reason stays `Unseparated` where the search ran out of halvings
along a curved overlap, since no budget would separate its pieces; it is
`BudgetExhausted` only where pieces beside the overlap were left
unsettled. A
stretch where the operands only touch (they part as the square of the
distance from the contact) is one hit at the closest approach, with
multiplicity 2, and the result stays `Complete`; an overlap between curves whose
range on the second curve would cross that curve's seam is not joined,
and keeps its hits and region.

The search reports whether it finished: the result's method is
`Numeric`, and its completeness is `Complete` when every piece pair
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
smaller than the iso-curve spacing can be missed. Two bounded surfaces
whose boxes are apart are proved not to meet: the result is empty and
`Complete`. The boxes of the control points are checked before the
march, and the tight boxes (see `tight_bounding_box`) after a march that
found nothing. Where a branch of two
unbounded surfaces was cut at its box, the result is `Incomplete` with
`Unbounded` as its reason, and its regions mark where the cuts are.
Where a curve's `deviation` exceeds the confusion, the result is
`Incomplete` with `LooseFit` (`LOOSE_FIT` in Python), and a region holds
each such curve: it is returned, but lies only within its deviation of
the surfaces. That happens where a surface's parameters degenerate on
the curve (a sphere's pole), where their rounding is coarser than the
confusion, and along a touch.

The confusion also sets the size of a fitted curve: each segment must
stay within it of both surfaces, so a looser confusion takes fewer,
longer segments. The steps never grow past a thousandth of the smaller
surface's size, so a looser confusion finds the same curves. Two wavy
B-spline patches cut at the default `1e-7` give curves of 298 and 198
poles; at `1e-3` the same two curves take 52 and 34 poles. The same
value decides what counts as a hit or a tangency, so loosen it for a
lighter fit only where a point within it of both surfaces is on both. In Python, `tol=1e-3` is short for
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

### The parameter of a fitted curve

A fitted curve, marched or fitted from a closed form (a loop where two
cylinders cross), is a cubic `BSplineCurve` whose parameter is its arc
length, to within the fit: its knots are the lengths of its segments,
each measured as a circular arc between its end tangents, so its speed
stays close to 1 (within 1e-4 on the loops of two crossing cylinders, a
torus and a sphere, and a plane and a B-spline sphere) and a piece's
length is about the length of its range. A fitted closed loop is
periodic, its period its length: a piece of it can run across its start,
as an arc of a circle runs across its seam, and a `TrimmedCurve` takes
the piece's range. Where the surfaces touch along a marched curve, it is
an interpolant whose parameter runs from 0 to 1, spaced by the square
roots of the distances between the traced points, and a loop is not
periodic.

### The curves in each surface's parameters

Each curve also comes as `first_uv_curve` and `second_uv_curve`: the
curve in the `(u, v)` parameters of each surface, in the same parameter
as `curve`, so the surface's point at the `(u, v)` of `t` is the curve's
point at `t`. A marched curve always has both: the march walked the
parameters of every point it traced, and they are fitted on the curve's
own knots, continuous across the seams (a loop around a cylinder ends a
period away in `u`). An exact curve has the exact one where a closed
form exists (a circle in a plane's parameters, a ruling or a parallel
circle as a line in a cylinder's) and `None` otherwise: an ellipse on a
cylinder, a circle of two spheres that is no parallel of either, a
hyperbola of two cones, the ellipses where two equal cylinders cross, and
any curve running through a sphere's pole or a cone's apex, where its
`u` jumps. A meridian clipped to a trimmed sphere between its poles has
its line again. [`curve_3d_to_uv_within`](parametrization.md) maps the
others within a deviation.

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
loop in ~8 ms on the same machine. The loop is an exact B-spline in the
scene, so the viewer refines it as you zoom:

```viewer
scene: assets/scenes/intersect-marching.json
title: B-spline patch and sphere with their marching intersection loop
```

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
patches of one surface are coincident only where they overlap with area.
Two that only touch report where: a shared edge is a `tangent` curve
with its curves in both patches' parameters, and a shared corner a
contact point. Two faces of one cylinder stacked at height 1 meet along
the circle there:

=== "Python"

    ```python
    import math
    from geomcore import Cylinder, Point, Vector

    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 1.0)
    low = cylinder.trimmed((0.0, 2 * math.pi), (0.0, 1.0))
    high = cylinder.trimmed((0.0, 2 * math.pi), (1.0, 2.0))
    result = low.intersect(high)
    assert not result.is_coincident()
    (edge,) = result.curves               # the circle at height 1
    assert edge.tangent
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{Cylinder, Point, Tolerance, Vector};
    use std::f64::consts::TAU;

    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let low = cylinder.trimmed((0.0, TAU), (0.0, 1.0)).unwrap();
    let high = cylinder.trimmed((0.0, TAU), (1.0, 2.0)).unwrap();
    let result = low.intersect(&high, Tolerance::DEFAULT);
    assert!(!result.is_coincident());
    assert_eq!(result.curves().len(), 1); // the circle at height 1
    ```

This needs the two bases to parametrize the surface alike, up to
shifted, reversed or swapped parameters, as two frames of one cylinder
on its axis do; patches of bases parametrized otherwise are still
coincident when they only touch (ADR 0039).

A curve meets a trimmed surface, in either order, as it meets the
basis, clipped to the ranges: a hit is kept when its surface parameters
lie in them, and reports them there, so a hit at a cylinder's angle 0
on a face from 3π/2 to 5π/2 reports `u = 2π`. An overlap, a stretch of
the curve lying on the basis, is cut to its pieces inside the ranges.
Filtering the basis's hits on `u_range()` by hand would miss that hit,
one period below the range.

=== "Python"

    ```python
    import math
    from geomcore import Cylinder, Line, Point, Vector

    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 1.0)
    face = cylinder.trimmed((1.5 * math.pi, 2.5 * math.pi), (0.0, 5.0))
    line = Line(Point(-5.0, 0.0, 1.0), Vector.unit_x())
    (hit,) = line.intersect(face).hits     # only the side x >= 0
    assert math.isclose(hit.u, 2 * math.pi)
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{Cylinder, Line, Point, Tolerance, Vector};
    use std::f64::consts::{PI, TAU};

    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let face = cylinder.trimmed((1.5 * PI, 2.5 * PI), (0.0, 5.0)).unwrap();
    let line = Line::new(Point::new(-5.0, 0.0, 1.0), Vector::X).unwrap();
    let result = line.intersect(&face, Tolerance::DEFAULT);
    assert_eq!(result.hits().len(), 1); // only the side x >= 0
    assert!((result.hits()[0].surface_parameters.0 - TAU).abs() < 1e-12);
    ```

A face moves with `transformed`: the basis maps through the transform
and each range with the basis's parametrization, scaled where the
parameter is a length (a cylinder's height) and kept where it is an
angle, so the image holds the image of every point of the face. Moving
the basis and trimming it to the same ranges is wrong under a scaling.
`face.trimmed(u_range, v_range)` trims a face again, over the part of
each range the new one covers; a periodic direction's range is read
modulo the period.

All the scenes on this page were computed by the kernel itself —
see `docs/assets/scenes/gen_scenes.py` — not drawn by hand.
