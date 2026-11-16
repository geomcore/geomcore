"""One rule for point arguments (#322): a point argument takes a `Point`
or a sequence of its coordinates (a tuple, a list, a row of an array),
and a list of points takes a list of those or an `(N, 3)` array, empty
included."""

import numpy
import numpy.typing
import pytest

from geomcore import (
    BezierCurve,
    BezierSurface,
    BoundingBox,
    BSplineCurve,
    BSplineCurve2D,
    BSplineSurface,
    Circle,
    Line,
    Plane,
    Point,
    Point2D,
    Sphere,
    Transform,
    Vector,
)

_POLES = numpy.array([[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [3.0, 2.0, 0.0], [4.0, 0.0, 0.0]])
_LIST = [Point(*row) for row in _POLES]


def test_a_point_argument_takes_a_point_a_tuple_a_list_or_an_array_row() -> None:
    circle = Circle(Point(0, 0, 0), Vector(0, 0, 1), 1.0)
    expected = circle.project_point(Point(3, 0, 0))
    for point in ((3, 0, 0), [3.0, 0.0, 0.0], numpy.array([3.0, 0.0, 0.0]), numpy.arange(3, 6)[[0, 1, 1]] * [1, 0, 0]):
        assert circle.project_point(point).point == expected.point
    assert Line((0, 0, 0), Vector(1, 0, 0)) == Line(Point(0, 0, 0), Vector(1, 0, 0))
    plane = Plane((0, 0, 0), Vector(0, 0, 1))
    assert plane.signed_distance(numpy.array([0.0, 0.0, 2.0])) == 2.0
    assert Sphere((0, 0, 0), 1.0).contains((1, 0, 0))
    move = Transform.translation(Vector(1, 0, 0))
    assert move.apply_point((0, 0, 0)) == Point(1, 0, 0)
    assert (0.5, 0.5, 0.5) in BoundingBox((0, 0, 0), (1, 1, 1))


def test_a_point_argument_refuses_a_vector_and_a_wrong_length() -> None:
    circle = Circle(Point(0, 0, 0), Vector(0, 0, 1), 1.0)
    for refused in (Vector(3, 0, 0), (3, 0), "abc", Point2D(1, 2)):
        with pytest.raises(TypeError, match="expected a Point or a sequence of 3 numbers"):
            circle.project_point(refused)  # type: ignore[arg-type]


def test_a_list_of_points_takes_points_tuples_rows_or_an_array() -> None:
    expected = BSplineCurve.clamped_uniform(3, _LIST)
    for poles in (_POLES, _POLES.tolist(), [tuple(row) for row in _POLES], list(_POLES)):
        assert BSplineCurve(3, poles, [0.0, 1.0], [4, 4], False) == expected
        assert BSplineCurve.interpolate(poles, 3) == BSplineCurve.interpolate(_LIST, 3)
    assert BSplineCurve.new_rational(3, _POLES, [1.0] * 4, [0.0, 1.0], [4, 4], False).degree() == 3
    assert BezierCurve(3, _POLES) == BezierCurve(3, _LIST)
    assert BoundingBox.from_points(_POLES) == BoundingBox.from_points(_LIST)
    assert BSplineCurve2D.clamped_uniform(3, _POLES[:, :2]).degree() == 3
    projections = Circle(Point(0, 0, 0), Vector(0, 0, 1), 1.0).project_points(_POLES)
    assert len(projections) == 4


def test_an_empty_list_of_points_is_no_points() -> None:
    circle = Circle(Point(0, 0, 0), Vector(0, 0, 1), 1.0)
    no_points: list[Point] = []
    empties: list["list[Point] | numpy.typing.NDArray[numpy.float64]"] = [no_points, numpy.empty((0, 3))]
    for empty in empties:
        parameters, points, distances = circle.project_points_array(empty)
        assert len(parameters) == len(points) == len(distances) == 0
        assert circle.project_points(empty) == []


def test_a_grid_of_points_takes_rows_or_an_array() -> None:
    grid = numpy.array([[[i, j, i * j] for j in range(2)] for i in range(2)], dtype=float)
    expected = BSplineSurface(1, 1, grid, [0.0, 1.0], [2, 2], [0.0, 1.0], [2, 2], False, False)
    rows = [[Point(*point) for point in row] for row in grid]
    assert BSplineSurface(1, 1, rows, [0.0, 1.0], [2, 2], [0.0, 1.0], [2, 2], False, False) == expected
    assert BSplineSurface(1, 1, list(grid), [0.0, 1.0], [2, 2], [0.0, 1.0], [2, 2], False, False) == expected
    assert BezierSurface(1, 1, grid) == BezierSurface(1, 1, rows)
    uvs = numpy.array([[0.25, 0.5], [1.0, 1.0]])
    assert expected.points_at(uvs) == expected.points_at([(0.25, 0.5), (1.0, 1.0)])


def test_a_surface_rebuilds_from_its_flat_knots() -> None:
    grid = numpy.array([[[i, j, (i - 1) * (j - 1)] for j in range(3)] for i in range(3)], dtype=float)
    surface = BSplineSurface(2, 2, grid, [0.0, 1.0], [3, 3], [0.0, 1.0], [3, 3], False, False)
    assert surface.u_flat_knots() == [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    rebuilt = BSplineSurface.from_flat_knots(2, 2, grid, surface.u_flat_knots(), surface.v_flat_knots())
    assert rebuilt == surface
    weights = [[1.0, 2.0, 1.0]] * 3
    rational = BSplineSurface.from_flat_knots(2, 2, grid, surface.u_flat_knots(), surface.v_flat_knots(), weights)
    assert rational.is_rational()


def test_point_of_an_array_says_to_unpack_it() -> None:
    with pytest.raises(TypeError, match=r"Point\(\*row\)"):
        Point(numpy.array([1.0, 2.0, 3.0]))  # type: ignore[arg-type, call-arg]
    with pytest.raises(TypeError, match=r"Point2D\(\*row\)"):
        Point2D((1.0, 2.0))  # type: ignore[arg-type, call-arg]
    assert Point(*numpy.array([1.0, 2.0, 3.0])) == Point(1, 2, 3)


def test_operators_take_points_not_sequences() -> None:
    # #325: an operator reads its operand as a point or a vector only.
    with pytest.raises(TypeError):
        Point(1, 2, 3) - (1, 2, 3)  # type: ignore[operator]
