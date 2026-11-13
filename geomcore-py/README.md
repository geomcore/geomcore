# geomcore (Python bindings)

Python bindings for [geomcore](https://github.com/geomcore/geomcore), a
geometric kernel for CAD-grade curves and surfaces, written in Rust:
intersection, projection and extrema, B-spline fitting, evaluation,
transformations and curve-on-surface parametrization. Every result says
how it was found and whether it is complete.

```python
from geomcore import BSplineCurve, Circle, Line, Plane, Point, Vector

circle = Circle(Point(0, 0, 0), Vector(0, 0, 1), 2.0)
plane = Plane(Point(1, 0, 0), Vector(1, 0, 0))

hits = circle.intersect(plane)  # closed form: EXACT, COMPLETE
print([hit.point for hit in hits])  # (1, ±1.732, 0)

projection = circle.project_point(Point(3, 1, 0))
print(projection.point, projection.distance)

line = Line(Point(0, 0, 3), Vector(1, 1, 0))
print(circle.extrema(line).closest().distance)  # 3.0

curve = BSplineCurve.interpolate(
    [Point(0, 0, 0), Point(1, 2, 0), Point(3, 2, 1), Point(4, 0, 0)], 3
)
print(curve.point_at(0.5))
```

Documentation: <https://geomcore.github.io/geomcore/>, and the
[Python guide](https://geomcore.github.io/geomcore/python/).

## Install

```sh
pip install geomcore          # or geomcore[numpy] for the array calls
```

Pre-built wheels are published for Linux (x86_64, aarch64), macOS, and
Windows, for Python 3.11+ (abi3: one wheel covers every version).

## What's inside

- **Shapes.** Lines, circles, ellipses, parabolas, hyperbolas, B-spline
  curves (rational and periodic), Bézier curves and their 2D twins;
  planes, cylinders, cones, spheres, tori, B-spline and Bézier surfaces,
  and surfaces extruded or revolved from a curve. `TrimmedCurve` and
  `TrimmedSurface` restrict any of them to a range; `QuarticBranch` is a
  piece of the closed-form quartic two quadrics meet in.
- **Intersections** with `a.intersect(b)` between any two curves or
  surfaces, in either order. The result holds the hits with their
  parameters, point and multiplicity, the overlaps, the surface–surface
  curves and contact points, the `method` that ran (`EXACT` closed form
  or `NUMERIC` search) and its `completeness`. A result is a sequence of
  its primary items: `len`, iteration and indexing run over them.
- **Projection and extrema.** `shape.project_point(p)` gives the closest
  point, its parameters and the distance (`signed_distance` on a
  surface); `a.extrema(b).closest()` the closest approach of two curves.
- **Evaluation and local properties.** `point_at`, `points_at` and
  `derivative`; curve `tangent`, `curvature` and `curvature_center`;
  surface `normal` and `curvatures`; `length`, `discretize`, bounding
  boxes, closedness and continuity; `span_index` and `basis_functions`
  on B-splines.
- **Fitting.** `BSplineCurve.interpolate`, least-squares `approximate`,
  and `approximate_within` a deviation, curves and surfaces.
- **Curve-on-surface parametrization.** `surface.curve_3d_to_uv(curve)`
  and `surface.uv_curve_to_3d(uv_curve)` in closed form, and
  `surface.curve_3d_to_uv_within(curve, max_deviation)` and
  `surface.uv_curve_to_3d_within(uv_curve, max_deviation)` as one
  B-spline within a deviation, such as a helix on a cylinder.
- **Transformations** with one `Transform` type, applied with
  `transformed`.
- **numpy** (the `[numpy]` extra): `points_at_array` and the other array
  calls take and return arrays, and `numpy.array(discretization)` is the
  `(N, 3)` array of its points.
- **Pickling and copying** of every value, curve and surface.
- **Type stubs** (PEP 561) for editors and type checkers.

Evaluation, transformations and curve-on-surface parametrization are
validated against golden fixtures at CAD-kernel-grade tolerances (1e-7).
See the [repository](https://github.com/geomcore/geomcore) for the Rust
crate, the [changelog](https://github.com/geomcore/geomcore/blob/main/CHANGELOG.md)
and the [roadmap](https://github.com/geomcore/geomcore/blob/main/ROADMAP.md).

## License

Dual-licensed under MIT or Apache-2.0, at your option.
