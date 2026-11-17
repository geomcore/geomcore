# Getting started

This tutorial takes you from an empty environment to a curve you can
inspect in the browser. You will install `geomcore`, evaluate a point on
a circle, move the circle, intersect it with an ellipse, and draw the
result with `geomcore-viewer`. It takes about ten minutes.

`geomcore` is a geometry kernel for CAD-grade curves and surfaces —
parametric evaluation, rigid transformations, projection, intersection,
fitting, and curves on surfaces — in Rust with Python bindings. It is
not a CAD system: there is no topology, B-rep or boolean operation. It
is the geometry layer such a system would stand on, usable on its own.

`geomcore-viewer` is the companion viewer. It draws the exact objects
the kernel computes — not hand-sampled drawings — either live from a
script or notebook, or as a static scene embedded in a page like this
one. You do not need it for the kernel itself, but it is the fastest way
to see whether the geometry is what you meant.

The code below is Python. The Rust API uses the same names; the install
and viewer steps show both, and every reference page has Rust tabs.

You need Python 3.11+ or a Rust toolchain, `pip` or `cargo`, and a
browser for the viewer step.

## 1. Install

=== "Python"

    ```sh
    pip install geomcore
    ```

    Wheels cover Linux (x86_64, aarch64), macOS and Windows on Python
    3.11+. The package ships type stubs, so an editor shows the full API
    and every docstring.

=== "Rust"

    ```sh
    cargo add geomcore
    ```

    Add `geomcore` to a binary or library crate as usual. The API is at
    [docs.rs](https://docs.rs/geomcore).

The viewer is a separate install. The Python package currently draws
geomcore 0.44 through 0.47, so install the pair explicitly and let pip
resolve a version both accept:

=== "Python"

    ```sh
    pip install "geomcore==0.47.*" geomcore-viewer
    ```

=== "Rust"

    ```sh
    cargo add geomcore-viewer --features geomcore
    ```

    The Rust client 0.5 builds against geomcore 0.48. It does not start
    the viewer server itself: run `geomcore-viewer serve` once from the
    Python package. See [Visualize with the viewer](viewer.md) for the
    compatibility notes and the server setup.

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

To get the same picture from your own script, install the viewer as in
step 1 and send the objects to a local tab. The first `add` starts a
background server and opens the browser; re-running the script updates
the scene in place.

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
