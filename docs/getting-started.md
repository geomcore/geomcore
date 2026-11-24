# Getting started

Evaluate a point on a circle, move the circle, and find where it
crosses an ellipse. Run the Python steps in order. The last step shows
how to inspect the result in a browser from Python or Rust.

## 1. Install

=== "Python"

    ```sh
    pip install geomcore
    ```

    Use Python 3.11 or later.

=== "Rust"

    ```sh
    cargo add geomcore
    ```

    Add the dependency to a Cargo project. See the
    [Rust API reference](https://docs.rs/geomcore).

## 2. Evaluate a point

A circle is stored as a centre, a normal and a radius. Its parameter is
the angle from the circle's own x direction, so `π/4` is the point at
45°:

```python
import math

from geomcore import Circle, Point, Vector

circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
p = circle.point_at(math.pi / 4)
print(p)  # Point(1.4142135623730951, 1.414213562373095, 0.0)
```

The inverse map goes back from a point to a parameter:

```python
print(circle.parameter_of(p))  # 0.7853981633974483
assert circle.contains(p)
```

This forward/inverse pair is the shape of nearly every query in the
kernel. The [elementary geometry reference](geometry/elementary.md)
lists the domain and parameter of each of the ten closed forms.

## 3. Move it

Rigid transformations apply to any curve or surface and keep the
parameter meaningful. `transformed` returns a new shape; nothing is
mutated:

```python
from geomcore import Transform, Vector

moved = circle.transformed(Transform.translation(Vector(0.5, 0.0, 0.0)))
assert moved.contains(moved.point_at(0.0))
```

The moved circle still has radius 2; only its centre changed.

## 4. Ask where two curves meet

Build a 3 × 2 ellipse in the world frame and intersect it with the moved
circle. `intersect` works in either order and decides the method for
you. This pair has a closed form, so the result is exact and proved
complete:

```python
from geomcore import Ellipse, Frame, Transform, Vector

moved = circle.transformed(Transform.translation(Vector(0.5, 0.0, 0.0)))
ellipse = Ellipse.from_frame(Frame.world(), 3.0, 2.0)
result = moved.intersect(ellipse)

print(result.method)        # IntersectionMethod.EXACT
print(result.is_complete()) # True: every hit was found
for hit in result.hits:
    print(hit.point, hit.first_parameter, hit.second_parameter, hit.multiplicity)
```

The circle crosses the ellipse four times. Each hit carries the point
and the parameter on both curves; `multiplicity` is 1 for a crossing
and 2 or more for a tangency, so a tangent contact is not mistaken for
two hits. The [intersection guide](guides/intersection.md) covers
the numeric fallback and how to read `method` and `completeness`.

## 5. See it

Here is the same circle and ellipse, computed by the kernel and drawn by
the hosted viewer. Drag to pan, scroll to zoom, hover an object for its
parameters, and press `2` or `3` to switch between the plane and 3D
views. Click a hit to copy its coordinates.

```viewer
scene: assets/scenes/getting-started.json
title: A circle and an ellipse with their four exact crossings
```

To draw from your own script, follow the
[viewer installation and compatibility instructions](viewer.md#install).
Use a separate environment if the Python viewer requires an older kernel
than your project. For Rust, start the viewer server before running the
example.

=== "Python"

    ```python
    import geomcore_viewer as gv
    from geomcore import Ellipse, Frame, Transform, Vector

    moved = circle.transformed(Transform.translation(Vector(0.5, 0.0, 0.0)))
    ellipse = Ellipse.from_frame(Frame.world(), 3.0, 2.0)
    result = moved.intersect(ellipse)

    view = gv.Viewer("A circle and an ellipse")
    view.add(ellipse, label="E")
    view.add(moved, label="C")
    view.add(result, label="C ∩ E")   # hits, with the operands faded
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{Circle, Ellipse, Frame, Point, Tolerance, Transform, Vector};
    use geomcore_viewer::{Geom, Opts, Viewer};

    let mut view = Viewer::new("A circle and an ellipse")?;
    let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 2.0)?;
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0)?
        .transformed(Transform::translation(Vector::new(0.5, 0.0, 0.0)))?;
    view.add(&ellipse, Opts::new().label("E"))?;
    view.add(&circle, Opts::new().label("C"))?;
    let result = circle.intersect(&ellipse, Tolerance::DEFAULT);
    view.add_all(
        result.hits().iter().map(|hit| Geom::from(hit.point)),
        Opts::new().label("C ∩ E"),
    )?;
    ```

If you would rather not run anything, `view.save("scene.json")` writes
the scene and the hosted viewer opens it from a URL or by drag and drop.
See [Visualize with the viewer](viewer.md) for the embed pattern used
on this page and the full list of controls.

## 6. Where next

You have installed the kernel, evaluated and moved a curve, intersected
two curves exactly, and drawn the result. From here:

- **Do a task** — [Projection & containment](guides/projection.md),
  [Intersection](guides/intersection.md),
  [Discretize a curve](guides/discretization.md), or
  [Map a curve onto a surface](guides/parametrization.md).
- **Look something up** — the
  [elementary geometry reference](geometry/elementary.md), the
  [NURBS reference](nurbs/curves.md), or the
  [Python API notes](python.md).
- **Understand a design** — [NURBS: from Bézier to
  B-spline](nurbs/index.md), [Why ten special
  cases](explanation/why-ten-special-cases.md), or the
  [design decisions](adr/index.md).
