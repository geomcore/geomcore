# Python

```sh
pip install geomcore
```

Pre-built wheels are published for Linux (x86_64, aarch64), macOS and
Windows, for Python 3.9+ (abi3: one wheel covers every version). The
package ships type stubs (PEP 561), so IDEs and type checkers see the full
API.

=== "Python"

    ```python
    from geomcore import Point3D, Vector3D
    from geomcore.curves import Circle3D
    import math

    circle = Circle3D(Point3D.origin(), Vector3D.unit_z(), 2.0)
    point = circle.eval_point(math.pi / 4)
    ```

=== "Rust"

    ```rust
    use geomcore::{Point3D, Vector3D};
    use geomcore::curves::Circle3D;

    let circle = Circle3D::new(Point3D::ORIGIN, Vector3D::Z, 2.0).unwrap();
    let point = circle.eval_point(std::f64::consts::FRAC_PI_4);
    ```

## Local properties and one-sided evaluation

Curves and surfaces answer the same differential queries as in Rust,
with the tolerance as the keyword `tol=None` (`Tolerance.default()`):

- 3D curves: `tangent(u)`, `curvature(u)` and `curvature_center(u)`,
  `is_closed()`, `continuity(u)` and `min_continuity()`.
- Surfaces: `normal(u, v)` (`Plane.normal()` takes no parameters),
  `curvatures(u, v)` returning a `SurfaceCurvatures`, `is_u_closed()` /
  `is_v_closed()`, `continuity_u(u)` / `continuity_v(v)` and their
  `min_continuity_u()` / `min_continuity_v()`, and `reversed_u()` /
  `reversed_v()` with `u_reversed_parameter(u)` /
  `v_reversed_parameter(v)` mapping a parameter onto the reversed
  surface.
- B-splines: `eval_d1_sided`, `eval_d2_sided` and
  `eval_derivative_sided`, taking a `Side` (`Side.BELOW` or `Side.ABOVE`)
  per direction, to read each side of a knot or periodic seam.

Where Rust returns `None` (no tangent at a cusp, no normal or curvatures
at a cone's apex), Python returns `None`. A continuity query returns
`SMOOTH_CONTINUITY` (`2**32 - 1`, Rust's `u32::MAX`) where the shape is
smooth.

=== "Python"

    ```python
    from geomcore import SMOOTH_CONTINUITY, Point3D, Side
    from geomcore.curves import BSplineCurve3D

    poles = [Point3D(0.0, 0.0, 0.0), Point3D(1.0, 0.0, 0.0), Point3D(1.0, 1.0, 0.0)]
    kinked = BSplineCurve3D(1, poles, [0.0, 1.0, 2.0], [2, 1, 2], False)
    assert kinked.continuity(1.0) == 0
    assert kinked.continuity(0.5) == SMOOTH_CONTINUITY
    _, incoming = kinked.eval_d1_sided(1.0, Side.BELOW)  # (1, 0, 0)
    _, outgoing = kinked.eval_d1_sided(1.0, Side.ABOVE)  # (0, 1, 0)
    ```

=== "Rust"

    ```rust
    use geomcore::{BSplineCurve3D, Point3D, Side, Tolerance};

    let poles = vec![
        Point3D::new(0.0, 0.0, 0.0),
        Point3D::new(1.0, 0.0, 0.0),
        Point3D::new(1.0, 1.0, 0.0),
    ];
    let kinked = BSplineCurve3D::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
    assert_eq!(kinked.continuity(1.0, Tolerance::DEFAULT), 0);
    assert_eq!(kinked.continuity(0.5, Tolerance::DEFAULT), u32::MAX);
    let (_, incoming) = kinked.eval_d1_sided(1.0, Side::Below); // (1, 0, 0)
    let (_, outgoing) = kinked.eval_d1_sided(1.0, Side::Above); // (0, 1, 0)
    ```

## Reprs

Reprs print floats as Python does and show coordinates below `1e-15` as
`0.0`, so rounding noise such as `cos(pi / 2)` does not clutter them. A repr is
therefore for reading, not an exact round trip of such tiny coordinates.

A Python API reference generated from docstrings, and tutorials, are on
the roadmap. Until then, the Rust items on
[docs.rs](https://docs.rs/geomcore) mirror the Python surface one-to-one.
