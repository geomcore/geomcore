# Curve on surface

A curve that lies on a surface has an exact 2D image in the surface's
parameter space — the *UV curve*, also called a pcurve.
`surface.curve_3d_to_uv` computes it in closed form where one exists
(lines and circles on the five elementary surfaces, any curve on a
plane); `curve_3d_to_uv_within` gives one B-spline within a deviation
anywhere else. A surface intersection hands over the UV curves of its
curves on both surfaces already (see
[Intersection](intersection.md#the-curves-in-each-surfaces-parameters)).
Adapted from
[Parametric domain](https://avenirgeometry.com/articles/parametric-domain/).

## Why the rectangle matters

Kernels work in $(u, v)$ because the rectangle is flat: offsetting,
trimming and marching are all easier there than in space. But a uniform
grid in the rectangle is not uniform on the surface — the first
fundamental form stretches it — so every downstream algorithm must know
what the parameters measure (see [Elementary](../geometry/elementary.md#what-the-parameters-measure)).

<figure class="geom-fig">
<img src="../../assets/figures/parametric-domain/01-a-uniform-grid-drawn-on-the.svg" alt="A uniform grid drawn on the surface">
<figcaption>A uniform grid in (u, v) is not uniform on the surface: the metric stretches it.</figcaption>
</figure>

<iframe src="https://geomcore.github.io/geomcore-viewer/?scene=https://geomcore.github.io/geomcore/assets/scenes/param-grid.json" width="100%" height="480" style="border: 1px solid #333; border-radius: 8px;" allow="fullscreen" title="A uniform (u, v) grid on a sphere, one cell at the equator and one against the pole"></iframe>

## Analytic UV curves

A circle coaxial with a cylinder is a horizontal line in $(u, v)$; a
line in a plane is a line. The result is a `Line2D`/`Circle2D` (or a
`Curve2D` run-time value) that evaluates through the surface back onto
the 3D curve within `tol.confusion` — that check is what "lies on"
means, and anything farther off is `CurveNotOnSurface`.

=== "Python"

    ```python
    uv_curve = cylinder.curve_3d_to_uv(circle)   # Line2D in (u, v)
    ```

=== "Rust"

    ```rust
    let uv_curve = cylinder.curve_3d_to_uv(&circle.into(), tol).unwrap();
    ```

A few images across the surfaces:

<figure class="geom-fig">
<img src="../../assets/figures/parametric-domain/02-three-curves-on-a-cylinder.svg" alt="Curves on a cylinder">
<figcaption>On a cylinder: coaxial circles are horizontal lines, generators vertical ones.</figcaption>
</figure>

<iframe src="https://geomcore.github.io/geomcore-viewer/?scene=https://geomcore.github.io/geomcore/assets/scenes/param-cylinder.json" width="100%" height="480" style="border: 1px solid #333; border-radius: 8px;" allow="fullscreen" title="Three curves on a cylinder and their images in (u, v)"></iframe>

<figure class="geom-fig">
<img src="../../assets/figures/parametric-domain/04-three-curves-on-a-cone.svg" alt="Curves on a cone">
<figcaption>On a cone the same construction winds with the growing radius.</figcaption>
</figure>

<iframe src="https://geomcore.github.io/geomcore-viewer/?scene=https://geomcore.github.io/geomcore/assets/scenes/param-cone.json" width="100%" height="480" style="border: 1px solid #333; border-radius: 8px;" allow="fullscreen" title="Three curves on a cone and their images in (u, v)"></iframe>

<figure class="geom-fig">
<img src="../../assets/figures/parametric-domain/05-three-curves-on-a-sphere.svg" alt="Curves on a sphere">
<figcaption>On a sphere: parallels are horizontal, meridians vertical — and both pinch at the poles.</figcaption>
</figure>

<iframe src="https://geomcore.github.io/geomcore-viewer/?scene=https://geomcore.github.io/geomcore/assets/scenes/param-sphere.json" width="100%" height="480" style="border: 1px solid #333; border-radius: 8px;" allow="fullscreen" title="Three curves on a sphere and their images in (u, v)"></iframe>

<figure class="geom-fig">
<img src="../../assets/figures/parametric-domain/06-three-curves-on-a-torus.svg" alt="Curves on a torus">
<figcaption>On a torus both families close: the rectangle glues into the surface itself.</figcaption>
</figure>

<iframe src="https://geomcore.github.io/geomcore-viewer/?scene=https://geomcore.github.io/geomcore/assets/scenes/param-torus.json" width="100%" height="480" style="border: 1px solid #333; border-radius: 8px;" allow="fullscreen" title="Three curves on a torus and their images in (u, v)"></iframe>

## Within a deviation

Where no closed form exists, `curve_3d_to_uv_within` gives one UV curve:
a degree-5 B-spline in the 3D curve's own parameter, whose image on the
surface lies within the asked deviation of the 3D curve at the same
parameter. It does not jump at a seam: a curve that goes once around a
cylinder ends one period of `u` from where it starts. It mirrors
`uv_curve_to_3d_within`. A curve that starts or ends at a sphere's pole
or a cone's apex, where every `u` is the same point, takes the limit of
`u` from inside the curve: a meridian keeps its own `u` up to the pole.

In Python both are surface methods, beside the exact `curve_3d_to_uv`
and `uv_curve_to_3d`: `surface.curve_3d_to_uv_within(curve, d)` and
`surface.uv_curve_to_3d_within(uv_curve, d)`. In Rust they are the
generic free functions `curve_3d_to_uv_within(&surface, &curve, d, tol)`
and `uv_curve_to_3d_within(&surface, &uv_curve, d, tol)`.

=== "Python"

    ```python
    approximation = surface.curve_3d_to_uv_within(curve, 1e-6)
    uv_curve = approximation.curve        # BSplineCurve2D
    approximation.deviation               # <= 1e-6
    approximation.is_proven()             # False: sampled
    ```

=== "Rust"

    ```rust
    use geomcore::curve_3d_to_uv_within;

    let approximation = curve_3d_to_uv_within(&surface, &curve, 1e-6, tol)?;
    let uv_curve = approximation.curve();
    ```

### Recipe: a helix

A helix is a straight line in the UV space of a cylinder: one turn is
2π in `u` and one pitch in `v`. `uv_curve_to_3d_within` turns that line
into one B-spline within the deviation, so a helix needs no type of its
own:

```python
import math

from geomcore import Cylinder, Line2D, Point, Point2D, TrimmedCurve2D, Vector, Vector2D

radius, pitch, turns = 5.0, 2.0, 3.0
cylinder = Cylinder(Point.origin(), Vector.unit_z(), radius)
# The line's parameter is its length in the UV plane.
uv_line = Line2D(Point2D(0.0, 0.0), Vector2D(2 * math.pi, pitch))
uv_piece = TrimmedCurve2D(uv_line, 0.0, turns * math.hypot(2 * math.pi, pitch))
approximation = cylinder.uv_curve_to_3d_within(uv_piece, 1e-6)
helix = approximation.curve  # degree 9, 60 poles
approximation.is_proven()  # True: within 4.2e-7, proven
helix.end_point()  # Point(5.0, ~0, 6.0): three turns up
```

A negative `u` direction turns the other way, a left-handed helix. The
B-spline is an approximation; a caller that needs the exact curve, for
its exact length or projection, would need a `Helix` type of its own.

## Viewing UV curves

The viewer draws a UV curve twice: on the surface in 3D and in its own
$(u, v)$ pane. Pass the 2D curve with the already-added surface's
handle.

=== "Python"

    ```python
    s = view.add(cylinder, label="S")
    view.add(circle, label="C")
    view.pcurve(cylinder.curve_3d_to_uv(circle), s, label="C in (u, v)")
    ```

=== "Rust"

    ```rust
    use geomcore_viewer::{Geom, Opts};

    let s = view.add(Geom::from(&cylinder), Opts::new().label("S"))?;
    view.add(&circle, Opts::new().label("C"))?;
    view.pcurve(&cylinder.curve_3d_to_uv(&circle.into(), tol)?, &s, Opts::new().label("C in (u, v)"))?;
    ```
