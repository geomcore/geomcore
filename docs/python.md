# Python

```sh
pip install geomcore
```

Pre-built wheels are published for Linux (x86_64, aarch64), macOS and
Windows, for Python 3.9+ (abi3: one wheel covers every version). The
package ships type stubs (PEP 561), so IDEs and type checkers see the full
API, with the documentation of the methods every curve or surface shares.

=== "Python"

    ```python
    from geomcore import Point, Vector
    from geomcore.curves import Circle
    import math

    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    point = circle.point_at(math.pi / 4)
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{Circle, Point, Vector};

    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let point = circle.point_at(std::f64::consts::FRAC_PI_4);
    ```

## Local properties and one-sided evaluation

Curves and surfaces answer the same differential queries as in Rust,
with the tolerance as the keyword `tol=None` (`Tolerance.default()`).
Every `tol` also takes a number, read as `Tolerance(confusion=number)`,
so `p.coincides(q, 1e-3)` loosens the distance tolerance alone:

- 3D curves: `tangent(u)`, `curvature(u)` and `curvature_center(u)`,
  `is_closed()`, `continuity(u)` and `min_continuity()`, `reversed()`
  with `reversed_parameter(u)` mapping a parameter onto the reversed
  curve, and `transformed(transform)`.
- 2D curves: the same queries but `transformed`; `curvature(u)` is
  signed (positive where the tangent turns left). `bounding_box()` and
  `tight_bounding_box()` return a `BoundingBox2D`, or `None` when the
  curve is unbounded.
- Surfaces: `normal(u, v)` (`Plane.normal()` takes no parameters),
  `curvatures(u, v)` returning a `SurfaceCurvatures`, `is_u_closed()` /
  `is_v_closed()`, `continuity_u(u)` / `continuity_v(v)` and their
  `min_continuity_u()` / `min_continuity_v()`, `reversed_u()` /
  `reversed_v()` with `u_reversed_parameter(u)` /
  `v_reversed_parameter(v)` mapping a parameter onto the reversed
  surface, and `transformed(transform)`.
- B-splines: `first_derivatives_sided`, `second_derivatives_sided` and
  `derivative_sided`, taking a `Side` (`Side.BELOW` or `Side.ABOVE`)
  per direction, to read each side of a knot or periodic seam.

Where Rust returns `None` (no tangent at a cusp, no normal or curvatures
at a cone's apex), Python returns `None`. A continuity query returns the
order as an `int`, or `math.inf` where the shape is smooth (Rust's
`Continuity::Order(k)` and `Continuity::Smooth`), so `c >= 2` and
`min(...)` need no special case.

=== "Python"

    ```python
    import math

    from geomcore import Point, Side
    from geomcore.curves import BSplineCurve

    poles = [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(1.0, 1.0, 0.0)]
    kinked = BSplineCurve(1, poles, [0.0, 1.0, 2.0], [2, 1, 2], False)
    assert kinked.continuity(1.0) == 0
    assert kinked.continuity(0.5) == math.inf
    _, incoming = kinked.first_derivatives_sided(1.0, Side.BELOW)  # (1, 0, 0)
    _, outgoing = kinked.first_derivatives_sided(1.0, Side.ABOVE)  # (0, 1, 0)
    ```

=== "Rust"

    ```rust
    use geomcore::{BSplineCurve, Continuity, Point, Side, Tolerance};

    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
    ];
    let kinked = BSplineCurve::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
    assert_eq!(kinked.continuity(1.0, Tolerance::DEFAULT), Continuity::Order(0));
    assert_eq!(kinked.continuity(0.5, Tolerance::DEFAULT), Continuity::Smooth);
    let (_, incoming) = kinked.first_derivatives_sided(1.0, Side::Below); // (1, 0, 0)
    let (_, outgoing) = kinked.first_derivatives_sided(1.0, Side::Above); // (0, 1, 0)
    ```

## Values, numpy and pickling

Points and vectors are sequences of their coordinates. They unpack
(`x, y, z = p`), index (`p[0]`, `p[-1]`) and convert to numpy, so
`numpy.array(p)` is a float array of shape `(3,)` and
`numpy.array(points)` one of shape `(N, 3)`. Vectors also divide by a
number (`v / 2`) and give their length to `abs(v)`.

A `Transform` reads and writes 4×4 homogeneous matrices, so a robot pose
or a placement from another tool comes in and goes out directly.
`a @ b` composes like the matrix product, `b` first, so `a @ b` is
`b.then(a)`:

```python
from geomcore import Transform
import numpy

pose = Transform.from_matrix(matrix)   # four rows of four numbers
back = pose.inverse()                  # pose @ back is the identity, within rounding
numpy.array((pose @ back).to_matrix())
```

`from_matrix` raises `ConstructionError` for a matrix that shears or
scales the directions unequally, since a curve or surface keeps its shape
only under a uniform scale of a rotation or a mirror.

Every value and every curve and surface pickles, and copies with `copy` and
`copy.deepcopy`, so shapes cross to `multiprocessing` workers and into
caches. The copy is equal to the original, bit for bit. The intersection,
projection and approximation results do not pickle, because they keep a
reference to the shapes they were computed on. Values stay unhashable
(ADR 0009): `==` compares coordinates exactly.

## Reprs

Reprs print floats as Python does and show coordinates below `1e-15` as
`0.0`, so rounding noise such as `cos(pi / 2)` does not clutter them. A repr is
therefore for reading, not an exact round trip of such tiny coordinates.

A Python API reference generated from docstrings, and tutorials, are on
the roadmap. Until then, the Rust items on
[docs.rs](https://docs.rs/geomcore) mirror the Python surface one-to-one.
