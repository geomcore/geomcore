# geomcore (Python bindings)

Python bindings for [geomcore](https://github.com/geomcore/geomcore), a
pure, standalone geometric kernel for CAD-grade curves and surfaces —
parametric evaluation, rigid transformations, and analytic curve-on-surface
parametrization, implemented in Rust.

```python
from geomcore import Point, Vector
from geomcore.curves import Circle
import math

circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
point = circle.eval_point(math.pi / 4)
```

## Install

```sh
pip install geomcore
```

Pre-built wheels are published for Linux (x86_64, aarch64), macOS, and
Windows, for Python 3.9+ (abi3 — one wheel covers all versions).

## What's inside

- Evaluation of lines, circles, ellipses, parabolas, hyperbolas and
  B-spline curves (rational and periodic), plus planes, cylinders, cones,
  spheres, tori and B-spline surfaces, with bulk `eval_points` in a single
  native call per batch.
- Rigid transformations via a single `Transform` type, applied with
  `transformed`.
- Containment on every curve and surface; `project_point` and `extrema`
  return projection objects with the parameters, the closest point and
  the distance.
- Intersections with `a.intersect(b)` between any two curves or surfaces
  (B-splines and trimmed curves included), in either order. The result is
  an object: hits with their parameters, point and multiplicity, overlaps,
  curves and contact points, the `method` that ran (`EXACT` closed form or
  `NUMERIC` search) and its `completeness`. Numeric `curve_curve_extrema`
  returns its distance minima and completeness.
- Local properties: curve `tangent`, `curvature` and
  `curvature_center` (signed in 2D); surface `normal` and `curvatures`
  (principal, mean and Gaussian, with principal directions), each `None`
  where it degenerates. Closedness and continuity queries, 2D bounding
  boxes (`BoundingBox2D`), surface reversal, and one-sided B-spline
  evaluation at knots and seams (`Side`).
- Curve-on-surface parametrization: `parametrize_on` (analytic) and
  `parametrize_numeric`.
- B-spline `interpolate` and least-squares `approximate`.
- Type stubs (PEP 561) for IDEs and type checkers.

Evaluation, transformations and curve-on-surface parametrization are
validated against golden fixtures at CAD-kernel-grade tolerances (1e-7).
See the [repository](https://github.com/geomcore/geomcore) for full
documentation, the [roadmap](https://github.com/geomcore/geomcore/blob/main/ROADMAP.md),
and the Rust crate.

## License

Dual-licensed under MIT or Apache-2.0, at your option.
