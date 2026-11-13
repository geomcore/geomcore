# Projection & containment

Three questions, one distance: *is the point on the shape*
(`contains`), *where is the nearest point* (`project_point`), and
*where are all the local minima* (`extrema`). Every query takes a
`Tolerance` and decides by distance against `tol.confusion`. Adapted
from [Inversion and projection](https://avenirgeometry.com/articles/projection/).

In Rust the curve queries are methods of the `PointExtrema` and
`PointExtrema2D` traits, built on `extrema`, and the surface queries
methods of `ParametricSurface`; `use geomcore::prelude::*` brings them
all into scope. `parameter_of(point, tol)` (`parameters_of` on a
surface) is the inverse of `point_at` for a point on the shape: the
projection's parameter by default, a closed form on the elementary
geometries.

## The foot-point condition

A foot point is where the curve meets its own tangent at a right angle
— the distance is stationary exactly where $(C(u) - P) \cdot C'(u) =
0$. On the ten [elementary](../geometry/elementary.md) geometries this
stationarity equation has a closed form (quartic/cubic on the conics,
an `atan2` on the rest), so projection is exact. On B-splines there is
no closed form: the kernel subdivides Bézier pieces with a bounding-box
hierarchy and solves each convex piece with bounded Newton — every
minimum found, however close two basins sit.

<figure class="geom-fig">
<img src="../../assets/figures/projection/01-four-feet-on-an-ellipse.svg" alt="Four foot points on an ellipse">
<figcaption>Four feet on an ellipse: one query point, four stationary points, one nearest.</figcaption>
</figure>

An ellipse already shows why `extrema` (all minima) and
`project_point` (the nearest) are different calls: a single point can
have four feet, and the nearest is only one of them. Computed live
below — `extrema` on a 3×2 ellipse from a near-centre query, four feet:

<iframe src="https://geomcore.github.io/geomcore-viewer/?scene=https://geomcore.github.io/geomcore/assets/scenes/projection-feet.json" width="100%" height="480" style="border: 1px solid #333; border-radius: 8px;" allow="fullscreen" title="Ellipse with the four projection feet of one query point"></iframe>

<figure class="geom-fig">
<img src="../../assets/figures/projection/02-three-candidates-on-a-b-spline.svg" alt="Candidate feet on a B-spline">
<figcaption>On a B-spline the candidates come from every span — the search must cover all of them to never miss the nearest.</figcaption>
</figure>

## contains

A point *lies on* a curve or surface when it is within `tol.confusion`
of it, tested by re-evaluating the inverse parameter — never by
pointer or topology.

The Python examples on this page share these shapes:

```python
import numpy
from geomcore import Circle, Ellipse, Line, Point, Sphere, Vector

circle = Circle(Point(0, 0, 0), Vector.unit_z(), 2.0)
sphere = Sphere(Point(0, 0, 0), 2.0)
ellipse = Ellipse(Point(0, 0, 0), Vector.unit_z(), Vector.unit_x(), 3.0, 2.0)
line = Line(Point(0, 0, 3), Vector.unit_x())
us = numpy.linspace(0.0, 6.0, 7)
vs = numpy.linspace(-1.0, 1.0, 7)
```

=== "Python"

    ```python
    circle.contains(Point(2, 0, 0))          # True, on the radius-2 circle
    sphere.contains(Point(0, 0, 2.01))       # False, off the surface
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::Tolerance;

    let tol = Tolerance::DEFAULT;
    assert!(circle.contains(Point::new(2.0, 0.0, 0.0), tol));
    assert!(!sphere.contains(Point::new(0.0, 0.0, 2.01), tol));
    ```

## project_point

Returns the closest point: its parameter(s), the point itself and the
distance. A distance within `tol.confusion` snaps to `0.0`, so projection
agrees with `contains`: with the default tolerance (`confusion` `1e-7`) a
point `1e-8` off a cylinder projects at distance `0.0`, and one `1e-7` off
at its measured distance. Pass a smaller `Tolerance` to read smaller
distances. In Python the result is a `CurveProjection`,
`CurveProjection2D` or `SurfaceProjection` that also holds the `target`
it was computed on and the `query_point`. Batch variants
(`project_points`) share the setup and return a list.

=== "Python"

    ```python
    foot = circle.project_point(Point(3, 0, 0))
    foot.parameter, foot.point, foot.distance, foot.target is circle
    sfoot = sphere.project_point(Point(0, 0, 5))
    sfoot.u, sfoot.v, sfoot.point, sfoot.distance
    ```

=== "Rust"

    ```rust
    let foot = circle.project_point(Point::new(3.0, 0.0, 0.0), tol);
    assert_eq!((foot.parameter, foot.distance), (0.0, 1.0));
    assert_eq!(foot.point, Point::new(2.0, 0.0, 0.0));

    let sfoot = sphere.project_point(Point::new(0.0, 0.0, 5.0), tol);
    ```

## extrema

All local minima of the distance, nearest first, and nothing else: a
farthest point is never listed, and an end of a bounded curve only
where the distance grows from it into the curve. What `project_point`
returns first, `extrema` returns whole.

=== "Python"

    ```python
    minima = ellipse.extrema(Point(4, 1, 0))   # [CurveProjection, ...] nearest first
    ```

=== "Rust"

    ```rust
    let minima = ellipse.extrema(Point::new(4.0, 1.0, 0.0), tol);
    ```

On seams the parameters are normalised into `bounds()`; a hit within
confusion of the seam reports 0 whichever side it was found from, so
callers can test results against `bounds()` and trimmed ranges without
losing hits.

## Many points at once

In Python, `project_points_array` takes an `(N, 3)` array-like (`(N, 2)`
for a planar curve) and returns numpy arrays instead of one projection
object per point; `points_at_array` does the same for evaluation. The
computation runs with the interpreter lock released, and numpy is
imported only when these are called (`pip install geomcore[numpy]`).

```python
import numpy

queries = numpy.random.default_rng(0).normal(size=(100_000, 3))
parameters, points, distances = circle.project_points_array(queries)
# parameters (N,), points (N, 3), distances (N,); a surface gives (N, 2) parameters
grid = sphere.points_at_array(numpy.column_stack([us, vs]))   # (N, 3)
```

## Curves against curves

`curve.extrema(other)` with another curve reports the local minima of
the distance between the two curves over their domains: the primitive
behind closest-approach queries and trimmed-arc handling. With a point,
the same method gives that point's projections, as above. In Rust the
curve pair is the free function `curve_curve_extrema` (ADR 0033).

=== "Python"

    ```python
    minima = line.extrema(circle).minima
    distances = [m.distance for m in minima]
    ```

=== "Rust"

    ```rust
    use geomcore::curve_curve_extrema;

    let minima = curve_curve_extrema(&line, &circle, tol);
    let closest = minima.minima().iter().map(|m| m.distance);
    ```
