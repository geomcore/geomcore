"""The predicates and measures of the value classes (#198): angles,
parallel and perpendicular tests, components, the axis queries, the
plane's signed distance, frame placement of directions, the transform
queries and the box measures. The tolerance comes from the caller, as an
angular bound or a `Tolerance`, with the default one when omitted."""

import math

from geomcore import (
    Axis,
    Axis2D,
    BoundingBox,
    BoundingBox2D,
    ConstructionError,
    Frame,
    Frame2D,
    Point,
    Point2D,
    Tolerance,
    Transform,
    Vector,
    Vector2D,
)
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


def test_points_give_midpoints_components_and_images():
    a, b = Point(1.0, -2.0, 3.0), Point(3.0, 4.0, -1.0)
    assert a.midpoint(b) == Point(2.0, 1.0, 1.0)
    assert a.components() == (1.0, -2.0, 3.0)
    shift = Transform.translation(Vector(0.0, 0.0, 1.0))
    assert a.transformed(shift) == shift.apply_point(a)
    assert Point2D(1.0, 2.0).midpoint(Point2D(3.0, 0.0)) == Point2D(2.0, 1.0)
    assert Point2D(1.0, 2.0).components() == (1.0, 2.0)


def test_frames_place_directions_and_move():
    frame = Frame.from_z(Point(5.0, 5.0, 5.0), Vector(0.0, 0.0, 2.0))
    assert frame.vector_at(0.0, 0.0, 2.0) == Vector(0.0, 0.0, 2.0)
    assert frame.local_components(Vector(0.0, 0.0, 2.0))[2] == 2.0
    moved = frame.transformed(Transform.translation(Vector(1.0, 0.0, 0.0)))
    assert moved.origin() == Point(6.0, 5.0, 5.0)
    try:
        frame.transformed(Transform.scaling(Point.origin(), 0.0))
    except ConstructionError:
        pass
    else:
        raise AssertionError("a collapsing transform has no image frame")

    plane = Frame2D.from_x(Point2D(1.0, 1.0), Vector2D(0.0, 1.0))
    assert plane.point_at(2.0, 0.0) == Point2D(1.0, 3.0)
    assert plane.local_coordinates(Point2D(1.0, 3.0)) == (2.0, 0.0)
    assert plane.vector_at(1.0, 0.0) == Vector2D(0.0, 1.0)
    assert plane.local_components(Vector2D(0.0, 1.0)) == (1.0, 0.0)
    assert plane.x_axis().direction() == plane.x_direction()
    assert plane.y_axis().origin() == plane.origin()


def test_transform_between_frames_and_its_queries():
    source = Frame.from_z(Point(1.0, 0.0, 0.0), Vector.unit_z())
    target = Frame(Point(0.0, 2.0, 0.0), Vector.unit_x(), Vector.unit_y())
    t = Transform.between_frames(source, target)
    image = t.apply_point(source.point_at(1.0, 2.0, 3.0))
    assert image.distance(target.point_at(1.0, 2.0, 3.0)) < 1e-12
    assert t.preserves_orientation()
    assert math.isclose(t.scale_factor(), 1.0)
    assert not Transform.mirror_plane(Frame.world()).preserves_orientation()
    assert Transform.scaling(Point.origin(), -2.0).scale_factor() == 2.0
    back = Transform.between_frames(target, source)
    assert t.then(back).is_identity()
    assert not t.is_identity(Tolerance(confusion=1.0))


def test_boxes_grow_intersect_and_measure():
    box = BoundingBox(Point.origin(), Point(2.0, 2.0, 2.0))
    assert Point(1.0, 1.0, 1.0) in box
    assert Point(3.0, 1.0, 1.0) not in box
    assert box.including(Point(5.0, 1.0, 1.0)).max() == Point(5.0, 2.0, 2.0)
    other = BoundingBox(Point(1.0, 1.0, 1.0), Point(3.0, 3.0, 3.0))
    assert box.intersection(other) == BoundingBox(Point(1.0, 1.0, 1.0), Point(2.0, 2.0, 2.0))
    assert box.intersection(BoundingBox(Point(4.0, 0.0, 0.0), Point(5.0, 1.0, 1.0))) is None
    assert box.distance(Point(5.0, 1.0, 1.0)) == 3.0
    assert box.distance(Point(1.0, 1.0, 1.0)) == 0.0

    square = BoundingBox2D(Point2D.origin(), Point2D(2.0, 2.0))
    assert Point2D(1.0, 1.0) in square
    assert square.distance(Point2D(5.0, 6.0)) == 5.0
    assert square.intersection(BoundingBox2D(Point2D(2.0, 0.0), Point2D(3.0, 1.0))) == BoundingBox2D(
        Point2D(2.0, 0.0), Point2D(2.0, 1.0)
    )
    assert square.including(Point2D(-1.0, 0.5)).min() == Point2D(-1.0, 0.0)


def test_bounded_curves_have_ends():
    from geomcore.curves import BSplineCurve, Circle, Line, TrimmedCurve

    poles = [Point(0.0, 0.0, 0.0), Point(1.0, 2.0, 0.0), Point(3.0, 0.0, 1.0)]
    spline = BSplineCurve(2, poles, [0.0, 2.0], [3, 3], False)
    assert spline.start_point() == Point(0.0, 0.0, 0.0)
    assert spline.end_point() == Point(3.0, 0.0, 1.0)
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    assert circle.start_point().coincides(circle.end_point())
    segment = TrimmedCurve(Line(Point.origin(), Vector.unit_x()), 1.0, 3.0)
    assert segment.end_point() == Point(3.0, 0.0, 0.0)
    assert not hasattr(Line(Point.origin(), Vector.unit_x()), "start_point")
