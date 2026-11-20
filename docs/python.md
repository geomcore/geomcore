# Python

```sh
pip install geomcore
```

Pre-built wheels are published for Linux (x86_64, aarch64), macOS and
Windows, for Python 3.11+ (abi3: one wheel covers every version). The
package ships type stubs (PEP 561), so IDEs and type checkers see the full
API, with the docstring of every class, method and property.

Every name is at the package root, as in the Rust crate: the shapes
(`Circle`, `Plane`, `BSplineSurface`), the values (`Point`, `Frame`,
`Tolerance`), the results (`CurveCurveIntersection`), the enums and the
errors. Three more names stand for a whole kind of shape: `Curve`,
`Curve2D` and `Surface` are the unions of the classes of the Rust enums
of the same names, both for annotations and for `isinstance`.
`typing.get_args` gives the tuple of classes. As in Rust, a `TrimmedCurve` is a piece of a curve rather than
one of them, so code that takes either says so:

```python
from geomcore import Circle, Curve, Point, Surface, TrimmedCurve, Vector

def clearance(curve: Curve | TrimmedCurve, surface: Surface) -> float: ...

shape = Circle(Point(0, 0, 0), Vector(0, 0, 1), 2.0)
assert isinstance(shape, Curve)  # False for a TrimmedCurve
```

=== "Python"

    ```python
    from geomcore import Circle, Point, Vector
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

## Extrema and minimum distance

`extrema(other)` gives the local minima of the distance between two
shapes, closest first, with both parameters and both points of each:
a `CurveCurveExtrema`, `CurveSurfaceExtrema` or `SurfaceSurfaceExtrema`
by pair family, whichever side the call is made on. Its `completeness`
says whether the search proved there is no other minimum: a pair with a
closed form (a line or a circle against a plane, a line against a sphere
or a cylinder, two planes, two spheres, two parallel cylinders) is
complete; any other pair is searched numerically and reports
`SAMPLED`. `minimum_distance(other)` is the closest of them, or `None`.

=== "Python"

    ```python
    from geomcore import Line, Plane, Point, Sphere, Vector

    sphere = Sphere(Point(0, 0, 0), 1.0)
    floor = Plane(Point(0, 0, -3), Vector(0, 0, 1))
    clearance = sphere.minimum_distance(floor)
    assert abs(clearance.distance - 2.0) < 1e-12
    assert clearance.first_point.distance(Point(0, 0, -1)) < 1e-12   # on the sphere
    assert clearance.second_point.distance(Point(0, 0, -3)) < 1e-12  # on the floor

    rail = Line(Point(0, 0, 3), Vector(1, 0, 0))
    extrema = rail.extrema(sphere)                     # a CurveSurfaceExtrema
    assert extrema.is_complete()
    (closest,) = extrema
    assert closest.curve_point.distance(Point(0, 0, 3)) < 1e-12
    ```

=== "Rust"

    ```rust
    use geomcore::{Line, Plane, Point, Sphere, Tolerance, Vector};
    use geomcore::{curve_surface_extrema, surface_surface_extrema};

    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let floor = Plane::new(Point::new(0.0, 0.0, -3.0), Vector::Z).unwrap();
    let extrema = surface_surface_extrema(&sphere, &floor, Tolerance::DEFAULT);
    let clearance = extrema.closest().unwrap();
    assert!((clearance.distance - 2.0).abs() < 1e-12);

    let rail = Line::new(Point::new(0.0, 0.0, 3.0), Vector::X).unwrap();
    assert!(curve_surface_extrema(&rail, &sphere, Tolerance::DEFAULT).is_complete());
    ```

Two 2D curves have their extrema in the plane: a `CurveCurveExtrema2D`
of `CurveCurveExtremum2D` minima, with the same search and completeness
as two curves in space. A 2D curve against a 3D curve or a surface raises
`TypeError`.

=== "Python"

    ```python
    from geomcore import Circle2D, Line2D, Point2D, Vector2D

    wheel = Circle2D(Point2D(1, 3), 1.0)
    ground = Line2D(Point2D(0, 0), Vector2D(1, 0))
    gap = wheel.minimum_distance(ground)
    assert abs(gap.distance - 2.0) < 1e-12
    assert gap.second_point.distance(Point2D(1, 0)) < 1e-9
    ```

=== "Rust"

    ```rust
    use geomcore::{Circle2D, Line2D, Point2D, Tolerance, Vector2D, curve_curve_extrema};

    let wheel = Circle2D::new(Point2D::new(1.0, 3.0), 1.0).unwrap();
    let ground = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let gap = curve_curve_extrema(&wheel, &ground, Tolerance::DEFAULT);
    assert!((gap.closest().unwrap().distance - 2.0).abs() < 1e-12);
    ```

A point is an operand of a curve's `extrema`, in either tier, and of a
B-spline or Bézier surface's; other surfaces raise `TypeError` for a
point, as `project_point` gives its nearest point.

## Local properties and one-sided evaluation

Curves and surfaces answer the same differential queries as in Rust,
with the tolerance as the keyword `tol=None` (`Tolerance.default()`).
Every `tol` also takes a number, read as `Tolerance(confusion=number)`,
so `p.coincides(q, 1e-3)` loosens the distance tolerance alone:

- 3D curves: `tangent(u)`, `curvature(u)` and `curvature_center(u)`,
  `is_closed()`, `continuity(u)` and `min_continuity()`, `reversed()`
  with `reversed_parameter(u)` mapping a parameter onto the reversed
  curve, and `transformed(transform)`.
- 2D curves: the same queries, `transformed` taking a `Transform2D`;
  `curvature(u)` is signed (positive where the tangent turns left). `bounding_box()` and
  `tight_bounding_box()` return a `BoundingBox2D`, or `None` when the
  curve is unbounded. `offset(distance, max_deviation, tol=None)` moves
  the curve along its left normal and returns a `CurveOffset2D`: a line or
  a circle (or a piece of one) exactly, any other curve as a
  `BSplineCurve2D` within `max_deviation`, with the parameters of its
  cusps in `cusp_parameters`; it raises `CurveOffset2DError`.
- Surfaces: `normal(u, v)` (on a `Plane` the same at every point,
  with the normal axis also as `Plane.axis()`),
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

    from geomcore import BSplineCurve, Point, Side

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

For many points, keep them in an array and call the array forms:
`points_at_array` evaluates an array of parameters (of `(u, v)` rows on a
surface), `normals_array` gives a surface's normals at `(u, v)` rows,
`project_points_array` projects an `(N, 3)` array and
`Transform.apply_points` moves one, each returning numpy arrays without
building a Python object per point. They
and the list forms (`points_at`, `project_points`) let other Python
threads run during the work, so a thread pool over chunks of points
scales. A list of `Point` objects you already hold converts with
`numpy.array(points)`, or two to three times faster with
`numpy.array([p.components() for p in points])`.

A point argument takes a `Point` or a sequence of its coordinates: a
tuple, a list or a row of an array, so `circle.project_point((3, 0, 0))`
and `Line(array[0], direction)` work. A list of points takes a list of
those or an array of shape `(N, 3)` (`(N, 2)` in 2D), empty included, and
a grid of poles a list of rows or an array of shape `(rows, columns, 3)`.
A vector, a string or a sequence of the wrong length raises `TypeError`.
`Point` itself takes three numbers: write `Point(*row)` for a row.
Operators are stricter: `p - (1, 2, 3)` raises `TypeError`, since a tuple
does not say whether it is a point or a vector.

numpy's operators step aside for points and vectors, so a numpy number
scales a vector from either side (`numpy.float64(2) * v` is a `Vector`).
An array on either side of `+`, `-` or `*` raises `TypeError`, and `==`
against an array is `False`. Convert with `numpy.array(p)` to compute on
the coordinates as an array.

A `Discretization` converts to its points: `numpy.array(curve.discretize(0.01))`
is a float array of shape `(N, 3)` (`(N, 2)` for a `Discretization2D`).
Intersection and extrema results iterate, index and reverse over what
they hold, but do not convert: they hold hits, curves and minima, not
points. Take `numpy.array([hit.point for hit in result])` for the points
of the hits.

A `Transform` reads and writes 4×4 homogeneous matrices, so a robot pose
or a placement from another tool comes in and goes out directly.
`a @ b` composes like the matrix product, `b` first, so `a @ b` is
`b.then(a)`:

```python
from geomcore import Transform
import numpy

matrix = [                             # a quarter turn about z, then a move
    [0.0, -1.0, 0.0, 1.0],
    [1.0, 0.0, 0.0, 2.0],
    [0.0, 0.0, 1.0, 3.0],
    [0.0, 0.0, 0.0, 1.0],
]
pose = Transform.from_matrix(matrix)   # four rows of four numbers
back = pose.inverse()                  # pose @ back is the identity, within rounding
numpy.array((pose @ back).to_matrix())
```

`from_matrix` raises `ConstructionError` for a matrix that shears or
scales the directions unequally, since a curve or surface keeps its shape
only under a uniform scale of a rotation or a mirror.

`Transform2D`, its plane twin, moves the 2D curves, points, vectors, axes
and frames, and reads and writes 3×3 homogeneous matrices the same way.
Its `rotation(center, angle)` turns counterclockwise about a point, and
`mirror_axis(axis)` reflects across a line, which turns a counterclockwise
`Circle2D` clockwise (`is_counterclockwise()` is then false).

Each 2D class answers what its 3D twin does, under the same names, except
where the plane has nothing to say. A 2D shape has no z coordinate or
direction, and no normal: a 2D circle or arc picks its turning direction
with `is_counterclockwise()` and the `counterclockwise=` argument of
`TrimmedCurve2D.arc`, where `TrimmedCurve.arc` takes a normal. A 2D
frame is built from its x direction (`Frame2D.from_x`), a 2D angle is
signed without a reference (`Vector2D.signed_angle`), and `Transform2D`
mirrors across a line (`mirror_axis`) rather than a plane. Surfaces and
every result that involves one live in space only. `test_parity.py` keeps
this list row by row.

Every value and every curve and surface pickles, and copies with `copy` and
`copy.deepcopy`, so shapes cross to `multiprocessing` workers and into
caches. The copy is equal to the original, bit for bit. Every result
pickles and copies too: the intersections, projections, extrema,
curvatures, approximations and discretizations, and the hits, regions and
curves inside them. A result that holds its operands (`first`,
`second`, `target`) pickles them with it, so a worker can return an
intersection. A shape's pickle calls its public constructor, so it loads
in a later geomcore as long as that constructor takes the same
arguments. A result's pickle loads only in the geomcore version that
wrote it, since its parts change between releases; another version
refuses it with a `ValueError` naming both (ADR 0053). Values stay unhashable
(ADR 0009): `==` compares coordinates exactly.

## Reprs

A repr is Python that rebuilds the object: `eval(repr(x))` (with the
`geomcore` names in scope) gives an object equal to `x`. A shape prints
as its constructor called with keywords, `Circle(center=Point(0.0, 0.0,
0.0), normal=Vector(0.0, 0.0, 1.0), radius=2.0)`, when those keywords
rebuild it. When they cannot say where its parametrization starts, as
for a tilted or mirrored circle, it prints as the static constructor
pickling uses instead: `Circle.from_frame(Frame.from_axes(...), 2.0)`. A
B-spline prints its poles, knots and multiplicities. A list longer than
ten items prints its first eight, `...` and its last, so a long pole list
stays readable, and its repr fails to paste rather than paste as another
shape.

Reprs print floats as Python does and show coordinates below `1e-15` as
`0.0`, so rounding noise such as `cos(pi / 2)` does not clutter them. Such
a coordinate is the one place a repr pastes back to a slightly different
object. A tolerance keeps its tiny values.

A Python API reference generated from docstrings is on the roadmap.
Until then, the Rust items on
[docs.rs](https://docs.rs/geomcore) mirror the Python surface one-to-one,
and the [getting-started tutorial](getting-started.md) plus the
[how-to guides](viewer.md) show the calls in context.
