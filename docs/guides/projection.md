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

At an interior stationary point, the vector from the query point to
the curve is perpendicular to the tangent. Its parameter satisfies
$(C(u) - P) \cdot C'(u) = 0$. On the ten [elementary](../geometry/elementary.md) geometries this
stationarity equation has a closed form (quartic/cubic on the conics,
an `atan2` on the rest), so projection is exact. On B-splines there is
no closed form: the kernel subdivides Bézier pieces with a bounding-box
hierarchy and solves each convex piece with bounded Newton — every
minimum found, however close two basins sit.

<figure class="geom-fig">
<img src="../../assets/figures/projection/01-four-feet-on-an-ellipse.svg" alt="Four foot points on an ellipse">
<figcaption>Four feet on an ellipse: one query point, four stationary points, one nearest.</figcaption>
</figure>

The diagram shows four stationary points. For the query below, the
3 × 2 ellipse has two distance minima. `extrema` returns both,
nearest first, and `project_point` returns the nearest.

```viewer
scene: assets/scenes/projection-feet.json
title: Two distance minima on an ellipse for one query point
```

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
from geomcore import Circle, Ellipse, Line, Point, Sphere, Vector

circle = Circle(Point(0, 0, 0), Vector.unit_z(), 2.0)
sphere = Sphere(Point(0, 0, 0), 2.0)
ellipse = Ellipse(Point(0, 0, 0), Vector.unit_z(), Vector.unit_x(), 3.0, 2.0)
line = Line(Point(0, 0, 3), Vector.unit_x())
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
    query = Point(0.5, 0.2, 0.0)
    minima = ellipse.extrema(query)
    assert len(minima) == 2
    assert minima[0].point.distance(ellipse.project_point(query).point) < 1e-7
    ```

=== "Rust"

    ```rust
    let minima = ellipse.extrema(Point::new(0.5, 0.2, 0.0), tol);
    assert_eq!(minima.len(), 2);
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

us = numpy.linspace(0.0, 6.0, 7)
vs = numpy.linspace(-1.0, 1.0, 7)
queries = numpy.random.default_rng(0).normal(size=(100_000, 3))
parameters, points, distances = circle.project_points_array(queries)
# parameters (N,), points (N, 3), distances (N,); a surface gives (N, 2) parameters
grid = sphere.points_at_array(numpy.column_stack([us, vs]))   # (N, 3)
```

## Closest approach between shapes

Use shape-to-shape extrema to measure clearance. The result lists
local distance minima in increasing order, with the witness points on
both operands and the search's `completeness`. Use `minimum_distance`
in Python when you only need the closest pair found. Inspect
`extrema(other).completeness` when you need to know whether the search
proved the global minimum.

Python supports two curves in the same dimension, a 3D curve and a
surface, or two surfaces. In Rust, use `curve_curve_extrema`,
`curve_surface_extrema` or `surface_surface_extrema`.

=== "Python"

    ```python
    closest = line.minimum_distance(circle)
    assert closest is not None and abs(closest.distance - 3.0) < 1e-7

    result = sphere.extrema(line)
    closest = result.closest()
    assert closest is not None and abs(closest.distance - 1.0) < 1e-7
    ```

=== "Rust"

    ```rust
    use geomcore::curve_curve_extrema;

    let minima = curve_curve_extrema(&line, &circle, tol);
    let closest = minima.minima().iter().map(|m| m.distance);
    ```
