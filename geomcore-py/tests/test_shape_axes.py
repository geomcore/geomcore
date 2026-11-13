"""Every shape with an axis spells it the same way (#196): the
constructors take `axis_direction` (the circle keeps its `normal`), and
`axis()` returns an `Axis`."""

from geomcore import Axis, Circle, Cone, Cylinder, Point, Torus, Vector


def test_every_rotational_shape_has_an_axis():
    center, up = Point(1.0, 2.0, 0.0), Vector(0.0, 0.0, 2.0)
    shapes = [
        Circle(center, up, 1.0),
        Cylinder(center, axis_direction=up, radius=1.0),
        Cone(center, axis_direction=up, half_angle=0.3, reference_radius=1.0),
        Torus(center, axis_direction=up, major_radius=3.0, minor_radius=1.0),
    ]
    for shape in shapes:
        assert shape.axis() == Axis(center, up), type(shape).__name__
    assert Circle(center, -up, 1.0).axis().direction() == Vector(0.0, 0.0, -1.0)
