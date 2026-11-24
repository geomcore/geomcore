# Visualize with the viewer

[`geomcore-viewer`](https://github.com/geomcore/geomcore-viewer) draws
the curves and surfaces `geomcore` computes. Each scene opens in 2D
or 3D. When you draw a curve in a surface's `(u, v)` parameters, its
parametric plot sits beside the 3D view or below it in narrow embeds.
The viewer evaluates analytic geometry directly, so curves and
surfaces refine as you zoom.

Drag to orbit the demo scene, right-drag to pan, and scroll to zoom.
Hover an object to inspect it.

```viewer
title: The viewer's demo scene
height: 640
```

## Install

=== "Python"

    ```sh
    pip install geomcore-viewer
    ```

    The Python package draws objects from geomcore 0.44 through 0.47 and
    declares that range, so pip resolves a compatible geomcore with it.
    For the Python viewer examples on this site, use a separate
    environment with the supported kernel version. The other examples
    in the documentation use geomcore 0.50.

    ```sh
    pip install "geomcore==0.47.*" geomcore-viewer
    ```

    The viewer has no runtime dependency of its own beyond geomcore and
    ships the built frontend in the wheel.

=== "Rust"

    ```toml
    [dependencies]
    geomcore = "0.48"
    geomcore-viewer = { version = "0.6", features = ["geomcore"] }
    ```

    The Rust client 0.6 builds its `geomcore` feature against geomcore
    0.48. It never starts the server itself: run `geomcore-viewer serve`
    from the Python package once, or start any Python script that uses
    the viewer. If no server answers, the client warns once and keeps the
    objects locally, so `view.save("scene.json")` still works.

## Draw from your own code

`view.add` accepts curves, surfaces and points. In Python it also takes
the result objects of projections, intersections and extrema. The first
`add` starts a background server and opens a browser tab; later runs
reuse it, and re-running a script swaps the scene in place. From a
debugger or a notebook every `gv.add(x)` lands in the same tab.

=== "Python"

    ```python
    from geomcore import Circle, Ellipse, Frame, Point, Vector
    import geomcore_viewer as gv

    view = gv.Viewer("Ellipse ∩ circle")
    E = Ellipse.from_frame(Frame.world(), 3.0, 2.0)
    C = Circle(Point(2, 1, 0), Vector.unit_z(), 1.5)
    view.add(E, label="E")
    view.add(C, label="C")
    view.add(C.intersect(E), label="P")   # the hits, operands faded
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{Circle, Cylinder, Point, Tolerance, Vector};
    use geomcore_viewer::{Geom, Opts, Viewer};

    let mut view = Viewer::new("A circle on a cylinder")?;
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0)?;
    let ring = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 2.0)?;
    let s = view.add(
        Geom::from(&cylinder).window(None, Some([-2.0, 2.0])),
        Opts::new().label("S"),
    )?;
    view.add(&ring, Opts::new().label("C"))?;
    view.pcurve(
        &cylinder.curve_3d_to_uv(&ring.into(), Tolerance::DEFAULT)?,
        &s,
        Opts::new().label("C in (u, v)"),
    )?;
    ```

The viewer also draws a UV curve twice: on the surface in 3D and in its
own `(u, v)` pane. Pass the 2D curve with the already-added surface's
handle (`view.pcurve(...)` in both clients). B-spline surfaces, trimmed
surfaces and trimmed curves are sent exactly; `wire`, `hidden` and
`silhouette` give textbook line art. The Rust client does not map a
whole result object to a scene: add the points and curves it holds, as
`view.add_all(result.hits().iter().map(|hit| Geom::from(hit.point)), ...)`.

## Embed a scene

`view.save("scene.json")` writes a scene without starting a server. The
hosted viewer opens one from `?scene=<url>` or by drag and drop:

```text
https://geomcore.github.io/geomcore-viewer/?scene=https://example.org/scene.json
```

The scene URL must be reachable by the browser, and the server it comes
from must allow cross-origin reads (GitHub Pages does).

On this site, use the `viewer` block instead of hand-writing an iframe.
It keeps the scene path, accessible title, caption and full-screen link
in one place:

````markdown
```viewer
scene: assets/scenes/projection-feet.json
title: Two distance minima on an ellipse for one query point
height: 480
```
````

`scene` is relative to the published docs root (or an absolute URL),
`title` is both the iframe title and the visible caption, and `height`
is optional. The hook in `docs_hooks/hooks.py` expands the block into a
lazy-loaded `<figure class="gc-viewer">` pointing at the viewer URL
configured under `extra.viewer_url` in `mkdocs.yml`. `extra.scene_base_url`
is the published root relative scenes resolve against; keeping it
separate from `site_url` lets `mkdocs serve` preview the deployed scenes.
Point `extra.viewer_url` at your own deployment to use a self-hosted
build. The full-screen link opens the same scene in a new tab, which is
also the fallback on a phone or when the embed is blocked.

## What it draws

- The elementary curves and surfaces, 2D conics, and points.
- B-spline and Bézier curves and surfaces, exactly: degree, poles,
  weights, knots, multiplicities and periodicity. Trimmed B-spline
  pieces keep their range.
- Trimmed curves and surfaces, extrusions and revolutions, and bounding
  boxes.
- In Python, the result objects of `project_point`, `extrema`, `intersect`
  and `curve_curve_extrema`, with their operands faded and a note saying
  which method ran.
- Curves in a surface's parameters, via `pcurve`.
- Point clouds, meshes and polylines for anything outside the kernel's
  types.

Infinite surfaces (planes, cylinders, cones) take a window; pass `u=`
and `v=` to choose it, or add the finite objects first and let the
viewer fit the default. A B-spline surface draws its own knot domain
unless `u=`/`v=` narrows it. The full option list is in the
[viewer README](https://github.com/geomcore/geomcore-viewer#api).

## Controls

| Key | Does |
|---|---|
| `F` | fit the scene, and keep following it as it changes |
| `P` | presentation mode: header and panel hidden, lines and labels larger |
| `T` | light / dark theme, remembered by the browser |
| `U` | show or hide the `(u, v)` pane |
| `←` / `→` | previous / next lesson step |
| `Esc` | close the object list in narrow embeds, leave presentation mode, or clear the selection |

In a narrow embed, select **Objects** to open the object list.

In the plane view, hover an object for its type, label, parameter `t`
and `(x, y)`; click a point to copy its coordinates. In 3D, drag to
orbit, right-drag to pan and scroll to zoom.

Scene format and protocol details:
[protocol/SCENE.md](https://github.com/geomcore/geomcore-viewer/blob/main/protocol/SCENE.md).
