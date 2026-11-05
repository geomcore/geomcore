"""The predicates and measures of the value classes (#198): angles,
parallel and perpendicular tests, components, the axis queries and the
plane's signed distance. The tolerance comes from the caller, as an
angular bound or a `Tolerance`, with the default one when omitted."""

import math

from geomcore import Axis, Axis2D, Point, Point2D, Tolerance, Vector, Vector2D
from geomcore.surfaces import Plane


def test_vectors_answer_parallel_and_perpendicular_for_any_length():
    z = Vector(0.0, 0.0, 3.0)
    assert z.is_parallel(Vector(0.0, 0.0, -0.5))
    assert not z.is_parallel(Vector(1e-6, 0.0, 1.0))
    assert z.is_parallel(Vector(1e-6, 0.0, 1.0), angular=1e-5)
    assert z.is_perpendicular(Vector(2.0, 5.0, 0.0))
    assert not Vector.zero().is_parallel(z) and not Vector.zero().is_perpendicular(z)
    assert Vector2D(2.0, 0.0).is_parallel(Vector2D(-1.0, 0.0))
    assert Vector2D(0.0, 2.0).is_perpendicular(Vector2D(1.0, 0.0), Tolerance.default().angular)


def test_vector_angles_and_components():
    assert math.isclose(Vector.unit_x().angle(Vector(-1.0, 1.0, 0.0)), 3 * math.pi / 4)
    v = Vector(1.0, 2.0, 3.0)
    assert v.component_along(Vector(0.0, 0.0, 2.0)) == Vector(0.0, 0.0, 3.0)
    assert v.component_perpendicular_to(Vector.unit_z()) == Vector(1.0, 2.0, 0.0)
    a, b = Vector2D(1.0, 0.0), Vector2D(0.0, -1.0)
    assert math.isclose(a.signed_angle(b), -math.pi / 2)
    assert math.isclose(a.angle(b), math.pi / 2)
    assert (a.rotated(-math.pi / 2) - b).magnitude() < 1e-15


def test_axis_queries():
    axis = Axis(Point(1.0, 0.0, 0.0), Vector(0.0, 0.0, 2.0))
    p = Point(4.0, 4.0, -2.0)
    assert axis.point_at(-2.0) == axis.closest_point(p) == Point(1.0, 0.0, -2.0)
    assert axis.parameter_of(p) == -2.0
    assert axis.distance(p) == 5.0
    assert axis.contains(Point(1.0, 0.0, 9.0)) and not axis.contains(p)
    assert axis.is_coaxial(Axis(Point(1.0, 0.0, 5.0), Vector(0.0, 0.0, -1.0)))
    assert not axis.is_coaxial(Axis(Point(2.0, 0.0, 0.0), Vector.unit_z()))
    assert axis.is_parallel(Axis(Point(2.0, 0.0, 0.0), Vector.unit_z()))
    assert axis.reversed().direction() == -Vector.unit_z()
    flat = Axis2D(Point2D(0.0, 1.0), Vector2D(1.0, 0.0))
    assert flat.distance(Point2D(5.0, -2.0)) == 3.0
    assert flat.parameter_of(Point2D(5.0, -2.0)) == 5.0
    assert flat.contains(Point2D(-4.0, 1.0), Tolerance.default())
    assert flat.is_parallel(Axis2D(Point2D(0.0, 0.0), Vector2D(-1.0, 0.0)))
    assert flat.reversed().point_at(1.0) == Point2D(-1.0, 1.0)


def test_plane_signed_distance_follows_the_orientation():
    # The metrology task (M-05): how far each probe point is above the datum.
    datum = Plane(Point(0.0, 0.0, 1.0), Vector.unit_z())
    probes = [Point(0.0, 0.0, 1.5), Point(3.0, -2.0, 0.25)]
    assert [datum.signed_distance(p) for p in probes] == [0.5, -0.75]
    assert datum.reversed_u().signed_distance(probes[0]) == -0.5
