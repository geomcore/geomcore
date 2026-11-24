# Build and place geometry

Use a trimmed curve for an arc or segment, a Bézier curve for a profile
with control points, and a swept surface when the same profile runs
along a direction or around an axis.

## Build a bounded profile

`TrimmedCurve` keeps the basis curve and its parameter range. On a circle,
the bounds are angles in radians. Use `TrimmedCurve.segment` for a
straight edge between two points and `arc_through` when you know three
points on the arc.

```python
import math

from geomcore import (
    Axis, BezierCurve, Circle, ExtrusionSurface, Point, RevolutionSurface,
    Tolerance, Transform, TrimmedCurve, Vector,
)

circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
arc = TrimmedCurve(circle, 0.0, math.pi / 2)
assert arc.end_point().distance(Point(0.0, 2.0, 0.0)) < 1e-12

profile = BezierCurve(2, [Point(2, 0, 0), Point(3, 0, 1), Point(2, 0, 2)])
assert profile.bounds() == (0.0, 1.0)
```

A quadratic Bézier uses three control points. Its first and last points
lie on the curve. The middle point shapes the bend. Use
[interpolation](../nurbs/editing-fitting.md#fitting) when the curve must
pass through the supplied points.

## Sweep the profile

On an extrusion, `u` follows the profile and `v` measures distance
along the sweep. On a revolution, `u` is an angle in radians and `v`
follows the profile. Call `trimmed` to keep a rectangle in either
surface's parameter space.

```python
extrusion = ExtrusionSurface(profile, Vector.unit_y())
strip = extrusion.trimmed(profile.bounds(), (0.0, 4.0))
assert strip.point_at(0.0, 4.0) == Point(2, 4, 0)

axis = Axis(Point.origin(), Vector.unit_z())
revolution = RevolutionSurface(profile, axis)
quarter = revolution.trimmed((0.0, math.pi / 2), profile.bounds())
assert quarter.point_at(math.pi / 2, 0.0).distance(Point(0, 2, 0)) < 1e-12
```

In Rust, use `ExtrusionSurface::new(profile.into(), direction)` and
`RevolutionSurface::new(profile.into(), axis)` with
`geomcore::prelude::*` for the shared surface methods.

## Place the profile and query it

A [frame](../geometry/elementary.md#foundations) gives a shape a local
origin and axes. A transform moves existing geometry. In both Python
and Rust, `rotate.then(shift)` applies the rotation first. Use
`transformed` on the whole trimmed shape to carry its bounds through
the transform.

```python
axis = Axis(Point.origin(), Vector.unit_z())
rotate = Transform.rotation(axis, math.pi / 2)
shift = Transform.translation(Vector(10, 0, 0))
placed = arc.transformed(rotate.then(shift))

probe = Point(10, 2, 0.0001)
assert placed.contains(probe, tol=Tolerance(confusion=0.001))
assert not placed.contains(probe)
```

| Tolerance | Default | What it controls |
|---|---:|---|
| `confusion` | `1e-7` | Distance in model units for treating locations as coincident |
| `angular` | `1e-12` | Angular comparisons such as parallelism |
| `parametric` | `1e-9` | Comparisons in curve or surface parameter space |

Python accepts a number for `tol` as shorthand for setting
`confusion`. In Rust, start with `Tolerance::DEFAULT` or pass a
`Tolerance::new(confusion, angular, parametric)` validated at construction.
