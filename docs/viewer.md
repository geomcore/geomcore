# Viewer

Live [geomcore-viewer](https://github.com/geomcore/geomcore-viewer) instance
below — Babylon.js 3D pane plus a PixiJS 2D pane for parameter-space curves.
The demo scene shows an ellipse, a circle, a hyperbola and a line, with the points where they meet.

<iframe src="https://geomcore.github.io/geomcore-viewer/" width="100%" height="640" style="border: 1px solid #333; border-radius: 8px;" allow="fullscreen" title="geomcore-viewer live demo"></iframe>

## Install

=== "Python"

    **geomcore-viewer 0.4 draws geomcore 0.43 at the latest.** It reads
    shapes through the names geomcore had up to 0.43 (`eval_point`,
    `semi_angle`, `first`…). With geomcore 0.44 or later, `view.add`
    raises `AttributeError` for every shape but a point, and refuses a
    `TrimmedSurface`.
    [geomcore-viewer#38](https://github.com/geomcore/geomcore-viewer/issues/38)
    tracks the viewer release that follows the renames. Until it is out,
    draw from an environment pinned to geomcore 0.43:

    ```sh
    pip install "geomcore==0.43.*" geomcore-viewer
    ```

    geomcore 0.43 imports the shapes from `geomcore.curves` and
    `geomcore.surfaces` and spells several methods differently, so the
    viewer examples on this site, written for the current release, do not
    run there as they are.

    geomcore-viewer 0.4 has no dependencies and ships the built viewer
    in the wheel. It reads geomcore objects through their public
    getters rather than declaring a geomcore version, so pip does not
    warn about a geomcore it cannot draw. It works with geomcore 0.33
    through 0.43, under both the class names from 0.39 on (`Point`,
    `Circle`…) and the older ones (`Point3D`, `Circle3D`…);
    geomcore-viewer 0.3 cannot draw objects from geomcore 0.39 or later.
    Drawing result objects directly (projections, intersections,
    extrema) needs geomcore 0.36.

=== "Rust"

    ```toml
    [dependencies]
    geomcore = "0.36"
    geomcore-viewer = { version = "0.4", features = ["geomcore"] }
    ```

    The `geomcore` feature of geomcore-viewer 0.4 is built against
    geomcore 0.36. Paired with another geomcore minor version, Cargo
    pulls in both versions and the viewer does not accept the other
    version's types, so keep geomcore at 0.36 in code that uses the
    viewer until a viewer release supports a later one. The Rust client never starts the server:
    run `geomcore-viewer serve` from the Python package once.

## From your own code

The embedded page above is static. From a script, a debugger or a notebook, `view.add(...)` sends geomcore objects to a local viewer tab that updates live:

=== "Python"

    ```python
    from geomcore import Circle, Ellipse, Frame, Point, Vector
    import geomcore_viewer as gv

    view = gv.Viewer("Ellipse ∩ circle")          # starts the local viewer and opens a tab
    E = Ellipse.from_frame(Frame.world(), 3.0, 2.0)
    C = Circle(Point(2, 1, 0), Vector.unit_z(), 1.5)
    view.add(E, label="E")
    view.add(C, label="C")
    view.add(C.intersect(E), label="P")           # the result draws directly: hits P, operands faded
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{Circle, Cylinder, Point, Tolerance, Vector};
    use geomcore_viewer::{Geom, Opts, Viewer};

    let mut view = Viewer::new("A circle on a cylinder")?;
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0)?;
    let ring = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 2.0)?;
    let s = view.add(Geom::from(&cylinder).window(None, Some([-2.0, 2.0])), Opts::new().label("S"))?;
    view.add(&ring, Opts::new().label("C"))?;
    view.pcurve(&cylinder.curve_3d_to_uv(&ring.into(), Tolerance::DEFAULT)?, &s, Opts::new().label("C in (u, v)"))?;
    ```

Protocol details: [protocol/SCENE.md](https://github.com/geomcore/geomcore-viewer/blob/main/protocol/SCENE.md).
