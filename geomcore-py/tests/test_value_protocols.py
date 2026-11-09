"""The Python protocols of the value classes (#182): iteration, numpy,
division, `Transform` matrices and inverse, and pickling of every value
and shape class."""

import copy
import math
import pickle

import pytest

from geomcore import (
    Axis,
    Axis2D,
    BSplineCurve,
    BSplineCurve2D,
    BSplineSurface,
    BezierCurve,
    BezierCurve2D,
    BezierSurface,
    BoundingBox,
    BoundingBox2D,
    Circle,
    Circle2D,
    Cone,
    ConstructionError,
    Cylinder,
    Ellipse,
    Ellipse2D,
    ExtrusionSurface,
    Frame,
    Frame2D,
    Hyperbola,
    Hyperbola2D,
    Line,
    Line2D,
    Parabola,
    Parabola2D,
    Plane,
    Point,
    Point2D,
    QuarticBranch,
    RevolutionSurface,
    Sphere,
    Tolerance,
    Torus,
    Transform,
    TrimmedCurve,
    TrimmedCurve2D,
    Vector,
    Vector2D,
)

# Directions of no particular length, so every frame is rotated.
_TILTED = Vector(0.3, -1.7, 2.2)
_HINT = Vector(1.1, 0.4, -0.9)
_CENTER = Point(1.5, -2.0, 0.25)
_FRAME = Frame(_CENTER, _TILTED, _HINT)
_MIRROR = Transform.mirror_plane(Frame(Point(0.1, 0.2, 0.3), Vector(1, 2, 3), _HINT))
_FRAME_2D = Frame2D.from_x(Point2D(0.5, -1.0), Vector2D(0.8, 1.3))


def _spline_curve() -> BSplineCurve:
    poles = [Point(0, 0, 0), Point(1, 2, 0.5), Point(3, 2, 1), Point(4, 0, 0.3), Point(5, 1, 2)]
    return BSplineCurve(3, poles, [0.0, 0.4, 1.0], [4, 1, 4], False)


def _rational_periodic_curve() -> BSplineCurve:
    poles = [Point(1, 0, 0), Point(0, 1, 0.2), Point(-1, 0, 0), Point(0, -1, -0.2)]
    return BSplineCurve.new_rational(
        2, poles, [1.0, 0.7, 1.3, 0.9], [0.0, 1.0, 2.0, 3.0, 4.0], [1, 1, 1, 1, 1], True
    )


def _spline_surface() -> BSplineSurface:
    poles = [[Point(i, j, (i * j) % 3) for j in range(3)] for i in range(3)]
    return BSplineSurface(2, 2, poles, [0.0, 1.0], [3, 3], [0.0, 1.0], [3, 3], False, False)


def _quartic_branch():
    """A branch of the quartic where two wide cones meet, their axes
    crossing at a right angle a little apart."""
    return QuarticBranch(
        Cone(Point(0, 0, 0), Vector(0, 0, 1), 1.0, 0.5),
        Cone(Point(0, 0.3, 0), Vector(1, 0, 0), 1.0, 0.5),
        1,
    )


def _cone_cylinder_branch():
    """A branch of the curve where a cylinder whose axis runs along a
    generator of a cone meets it."""
    cone = Cone(Point(0, 0, 1), Vector(0, 0, 1), math.pi / 4, 1.0)
    along = Vector(1, 0, 1).normalized()
    return QuarticBranch(cone, Cylinder(Point(0, 0.5, 0), along, 1.0), 0)


def _shapes():
    profile = _spline_curve()
    planar = BSplineCurve2D(2, [Point2D(0, 0), Point2D(1, 2), Point2D(3, 1)], [0.0, 1.0], [3, 3], False)
    return [
        Line(_CENTER, _TILTED),
        Line2D(Point2D(1, 2), Vector2D(0.3, -0.7)),
        Circle.from_frame(_FRAME, 2.5),
        Circle.from_frame(_FRAME, 2.5).transformed(_MIRROR),
        Circle2D.from_frame(_FRAME_2D, 1.5),
        Circle2D.from_frame(_FRAME_2D, 1.5).reversed(),
        Ellipse.from_frame(_FRAME, 3.0, 1.25),
        Ellipse2D.from_frame(_FRAME_2D, 3.0, 1.25),
        Hyperbola.from_frame(_FRAME, 2.0, 0.5),
        Hyperbola2D.from_frame(_FRAME_2D, 2.0, 0.5),
        Parabola.from_frame(_FRAME, 0.75),
        Parabola2D.from_frame(_FRAME_2D, 0.75),
        profile,
        _rational_periodic_curve(),
        planar,
        BSplineCurve2D.new_rational(
            2, [Point2D(0, 0), Point2D(1, 2), Point2D(3, 1)], [1.0, 0.5, 1.0], [0.0, 1.0], [3, 3], False
        ),
        BezierCurve(3, [Point(0, 0, 0), Point(1, 2, 0), Point(3, 2, 1), Point(4, 0, 1)]),
        BezierCurve.new_rational(2, [Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)], [1.0, 0.7, 1.0]),
        BezierCurve2D(2, [Point2D(0, 0), Point2D(1, 2), Point2D(3, 1)]),
        _quartic_branch(),
        _quartic_branch().reversed().transformed(_MIRROR),
        TrimmedCurve(_quartic_branch(), -1.0, 2.0),
        _cone_cylinder_branch(),
        TrimmedCurve(Circle.from_frame(_FRAME, 2.5), 0.25, 2.0),
        TrimmedCurve(profile, 0.1, 0.8),
        TrimmedCurve2D(planar, 0.2, 0.9),
        Plane.from_frame(_FRAME),
        Plane.from_frame(_FRAME).transformed(_MIRROR),
        Cylinder.from_frame(_FRAME, 2.0),
        Cylinder.from_frame(_FRAME, 2.0).reversed_u(),
        Cone.from_frame(_FRAME, 0.4, 1.5),
        Sphere.from_frame(_FRAME, 3.0),
        Torus.from_frame(_FRAME, 4.0, 1.0),
        ExtrusionSurface(profile, _TILTED),
        RevolutionSurface(profile, Axis(Point(-3, 0, 0), _TILTED)),
        _spline_surface(),
        BezierSurface(1, 1, [[Point(0, 0, 0), Point(0, 1, 0)], [Point(1, 0, 1), Point(1, 1, 0)]]),
    ]


def _values():
    return [
        Point(1.5, -2.0, 1e-300),
        Point2D(-0.1, 7.0),
        Vector(0.1, 0.2, 0.3),
        Vector2D(3.0, -4.0),
        Axis(_CENTER, _TILTED),
        Axis2D(Point2D(1, 2), Vector2D(0.3, -0.7)),
        _FRAME,
        _FRAME.flipped_y(),
        _FRAME_2D,
        Frame2D(Point2D(0, 0), Vector2D(0.6, 0.8), Vector2D(0.8, -0.6)),
        Tolerance(confusion=1e-5, angular=1e-9, parametric=1e-7),
        Transform.rotation(Axis(_CENTER, _TILTED), 0.9).then(_MIRROR),
        Transform.scaling(Point(0, 0, 0), 0.0),
        BoundingBox(Point(0, 0, 0), Point(1, 2, 3)),
        BoundingBox2D(Point2D(0, 0), Point2D(1, 2)),
    ]


_EVERYTHING = _values() + _shapes()


@pytest.mark.parametrize("original", _EVERYTHING, ids=lambda o: type(o).__name__)
def test_every_value_and_shape_pickles_to_an_equal_object(original):
    for protocol in range(2, pickle.HIGHEST_PROTOCOL + 1):
        restored = pickle.loads(pickle.dumps(original, protocol=protocol))
        assert type(restored) is type(original)
        assert restored == original
    assert copy.deepcopy(original) == original
    assert copy.copy(original) == original


def test_points_and_vectors_unpack_index_and_convert_to_numpy():
    numpy = pytest.importorskip("numpy")
    point = Point(1.0, 2.0, 3.0)
    x, y, z = point
    assert (x, y, z) == (1.0, 2.0, 3.0)
    assert (len(point), point[0], point[-1]) == (3, 1.0, 3.0)
    with pytest.raises(IndexError):
        point[3]
    assert list(Vector2D(4.0, 5.0)) == [4.0, 5.0]
    assert numpy.array(point).tolist() == [1.0, 2.0, 3.0]
    assert numpy.array(point).dtype == numpy.float64
    assert numpy.array([point, Point(4, 5, 6)]).shape == (2, 3)
    assert numpy.asarray(Point2D(1, 2), dtype=numpy.float32).dtype == numpy.float32


def test_vectors_divide_and_have_an_absolute_value():
    assert Vector(3.0, 6.0, 9.0) / 3.0 == Vector(1.0, 2.0, 3.0)
    assert Vector2D(1.0, 1.0) / 4.0 == Vector2D(0.25, 0.25)
    assert abs(Vector(3.0, 4.0, 0.0)) == 5.0
    assert abs(Vector2D(-3.0, 4.0)) == 5.0
    with pytest.raises(ZeroDivisionError):
        Vector(1.0, 0.0, 0.0) / 0.0


def test_a_transform_round_trips_its_matrix_and_composes_with_matmul():
    turn = Transform.rotation(Axis(_CENTER, _TILTED), 0.9)
    shift = Transform.translation(Vector(1, 2, 3))
    assert Transform.from_matrix(turn.to_matrix()) == turn
    assert (turn @ shift) == shift.then(turn)
    p = Point(0.5, -1.0, 2.0)
    composed = (turn @ shift).apply_point(p)
    assert composed.distance(turn.apply_point(shift.apply_point(p))) < 1e-14
    numpy = pytest.importorskip("numpy")
    product = numpy.array(turn.to_matrix()) @ numpy.array(shift.to_matrix())
    assert numpy.allclose(product, numpy.array((turn @ shift).to_matrix()), rtol=0, atol=1e-15)
    assert Transform.from_matrix(numpy.array(shift.to_matrix())) == shift
    assert repr(shift) == (
        "Transform.from_matrix([[1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 2.0], "
        "[0.0, 0.0, 1.0, 3.0], [0.0, 0.0, 0.0, 1.0]])"
    )


def test_a_transform_inverse_undoes_it():
    pose = Transform.rotation(Axis(_CENTER, _TILTED), 0.9).then(_MIRROR)
    back = pose.inverse()
    for p in [Point(0, 0, 0), Point(3, -1, 2), Point(-7, 0.5, 9)]:
        assert back.apply_point(pose.apply_point(p)).distance(p) < 1e-13
    with pytest.raises(ConstructionError):
        Transform.scaling(Point(0, 0, 0), 0.0).inverse()


def test_a_matrix_that_is_not_a_similarity_is_refused():
    sheared = [[1.0, 0.5, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    with pytest.raises(ConstructionError):
        Transform.from_matrix(sheared)
    projective = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.1, 1.0]]
    with pytest.raises(ConstructionError):
        Transform.from_matrix(projective)
    with pytest.raises(ValueError):
        Transform.from_matrix([[1.0, 0.0], [0.0, 1.0]])


def test_a_robot_pose_matrix_moves_a_shape():
    angle = math.pi / 6
    pose = Transform.from_matrix(
        [
            [math.cos(angle), -math.sin(angle), 0.0, 10.0],
            [math.sin(angle), math.cos(angle), 0.0, 0.0],
            [0.0, 0.0, 1.0, 5.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    moved = Circle(Point(0, 0, 0), Vector(0, 0, 1), 1.0).transformed(pose)
    assert moved.point_at(0.0).distance(Point(10 + math.cos(angle), math.sin(angle), 5)) < 1e-14
