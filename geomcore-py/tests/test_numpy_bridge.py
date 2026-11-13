"""The numpy bridge: bulk queries that read and return float arrays."""

import math
from typing import List, Sequence, Union

import numpy
import pytest

from geomcore import (
    BSplineCurve,
    BSplineCurve2D,
    BSplineSurface,
    Circle,
    Circle2D,
    Point,
    Point2D,
    Sphere,
    Vector,
    Vector2D,
)


def _rows(points: Sequence[Union[Point, Point2D]]) -> List[List[float]]:
    """The coordinates of each point, as rows."""
    return [
        [point.x, point.y, point.z] if isinstance(point, Point) else [point.x, point.y]
        for point in points
    ]


def test_curve_points_at_array_matches_points_at() -> None:
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    us = numpy.linspace(0.0, 2.0 * math.pi, 7)
    array = circle.points_at_array(us)
    assert array.shape == (7, 3)
    assert array.dtype == numpy.float64
    numpy.testing.assert_allclose(array, _rows(circle.points_at(list(us))))


def test_curve_2d_points_at_array_matches_points_at() -> None:
    circle = Circle2D(Point2D.origin(), 1.5)
    array = circle.points_at_array([0.0, 1.0, 2.0])
    assert array.shape == (3, 2)
    numpy.testing.assert_allclose(array, _rows(circle.points_at([0.0, 1.0, 2.0])))


def test_surface_points_at_array_matches_points_at() -> None:
    sphere = Sphere(Point.origin(), 2.0)
    uvs = numpy.array([[0.0, 0.0], [1.0, 0.5], [2.0, -0.5]])
    array = sphere.points_at_array(uvs)
    assert array.shape == (3, 3)
    pairs = [(float(u), float(v)) for u, v in uvs]
    numpy.testing.assert_allclose(array, _rows(sphere.points_at(pairs)))


def test_curve_project_points_array_matches_project_points() -> None:
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    queries = numpy.array([[3.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 1.0, 2.0]])
    parameters, points, distances = circle.project_points_array(queries)
    projections = circle.project_points([Point(*row) for row in queries])
    assert parameters.shape == (3,)
    assert points.shape == (3, 3)
    assert distances.shape == (3,)
    numpy.testing.assert_allclose(parameters, [p.parameter for p in projections])
    numpy.testing.assert_allclose(points, _rows([p.point for p in projections]))
    numpy.testing.assert_allclose(distances, [p.distance for p in projections])


def test_curve_2d_project_points_array_has_planar_rows() -> None:
    circle = Circle2D(Point2D.origin(), 1.0)
    parameters, points, distances = circle.project_points_array([[2.0, 0.0], [0.0, 3.0]])
    assert points.shape == (2, 2)
    numpy.testing.assert_allclose(parameters, [0.0, math.pi / 2.0], atol=1e-12)
    numpy.testing.assert_allclose(points, [[1.0, 0.0], [0.0, 1.0]], atol=1e-12)
    numpy.testing.assert_allclose(distances, [1.0, 2.0])


def test_surface_project_points_array_matches_project_points() -> None:
    sphere = Sphere(Point.origin(), 2.0)
    queries = numpy.array([[3.0, 0.0, 0.0], [0.0, 1.0, 1.0]])
    parameters, points, distances = sphere.project_points_array(queries)
    projections = sphere.project_points([Point(*row) for row in queries])
    assert parameters.shape == (2, 2)
    numpy.testing.assert_allclose(parameters, [(p.u, p.v) for p in projections])
    numpy.testing.assert_allclose(points, _rows([p.point for p in projections]))
    numpy.testing.assert_allclose(distances, [p.distance for p in projections])


def _sample_points() -> numpy.ndarray:
    u = numpy.linspace(0.0, 1.0, 12)
    return numpy.column_stack([u, numpy.sin(3.0 * u), u * u])


def test_fitting_a_curve_from_an_array_equals_fitting_from_points() -> None:
    samples = _sample_points()
    from_array = BSplineCurve.interpolate(samples, 3)
    from_points = BSplineCurve.interpolate([Point(*row) for row in samples], 3)
    assert from_array.poles() == from_points.poles()
    assert from_array.knots() == from_points.knots()
    approximated = BSplineCurve.approximate(samples, 3, 6)
    assert approximated.poles() == BSplineCurve.approximate(
        [Point(*row) for row in samples], 3, 6
    ).poles()


def test_fitting_a_2d_curve_from_an_array_equals_fitting_from_points() -> None:
    samples = _sample_points()[:, :2]
    from_array = BSplineCurve2D.approximate_within(samples, 3, 1e-3)
    from_points = BSplineCurve2D.approximate_within(
        [Point2D(*row) for row in samples], 3, 1e-3
    )
    assert from_array.curve.poles() == from_points.curve.poles()


def test_fitting_a_surface_from_a_grid_equals_fitting_from_rows() -> None:
    u, v = numpy.meshgrid(numpy.linspace(0.0, 1.0, 6), numpy.linspace(0.0, 1.0, 5), indexing="ij")
    grid = numpy.stack([u, v, numpy.sin(u + v)], axis=-1)
    assert grid.shape == (6, 5, 3)
    rows = [[Point(*point) for point in row] for row in grid]
    from_array = BSplineSurface.approximate(grid, 3, 3, 5, 4)
    from_rows = BSplineSurface.approximate(rows, 3, 3, 5, 4)
    assert from_array.poles() == from_rows.poles()


def test_inputs_of_other_dtypes_and_layouts_are_converted() -> None:
    circle = Circle(Point.origin(), Vector.unit_z(), 1.0)
    as_integers = circle.project_points_array(numpy.array([[2, 0, 0], [0, 3, 0]]))
    transposed = numpy.asfortranarray(numpy.array([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0]]))
    for left, right in zip(as_integers, circle.project_points_array(transposed)):
        numpy.testing.assert_allclose(left, right)


def test_a_wrong_shape_raises_value_error() -> None:
    circle = Circle(Point.origin(), Vector.unit_z(), 1.0)
    with pytest.raises(ValueError, match=r"points must have shape \(N, 3\)"):
        circle.project_points_array(numpy.zeros((4, 2)))
    with pytest.raises(ValueError, match=r"us must have shape \(N\)"):
        circle.points_at_array(numpy.zeros((2, 2)))
    with pytest.raises(ValueError, match=r"points must have shape \(N, M, 3\)"):
        BSplineSurface.approximate(numpy.zeros((4, 4)), 3, 3, 4, 4)


def test_empty_inputs_give_empty_arrays() -> None:
    circle = Circle2D(Point2D.origin(), 1.0)
    assert circle.points_at_array([]).shape == (0, 2)
    parameters, points, distances = circle.project_points_array(numpy.zeros((0, 2)))
    assert (parameters.shape, points.shape, distances.shape) == ((0,), (0, 2), (0,))


def test_returned_arrays_are_fresh_and_writable() -> None:
    circle = Circle2D(Point2D.origin(), 1.0)
    first = circle.points_at_array([0.0])
    second = circle.points_at_array([0.0])
    first[0, 0] = 42.0
    assert second[0, 0] == pytest.approx(1.0)
    assert first.flags.writeable


@pytest.mark.parametrize(
    "point, vector",
    [(Point(1.0, 2.0, 3.0), Vector(1.0, 2.0, 3.0)), (Point2D(1.0, 2.0), Vector2D(1.0, 2.0))],
)
def test_numpy_on_the_left_defers_to_the_value_operators(point, vector) -> None:
    # #325: numpy read the values as sequences and its operators ran first.
    two = numpy.float64(2.0)
    doubled = type(vector)(*(2.0 * c for c in vector))
    for scaled in (two * vector, vector * two, vector / numpy.float64(0.5)):
        assert type(scaled) is type(vector) and scaled == doubled
    array = numpy.ones(len(point))
    for value in (point, vector):
        for operation in (
            lambda: value + array,
            lambda: array + value,
            lambda: value - array,
            lambda: array - value,
            lambda: value * array,
            lambda: array * value,
        ):
            with pytest.raises(TypeError):
                operation()
        same = numpy.array(value)
        assert (value == same) is False and (same == value) is False
        numpy.testing.assert_array_equal(numpy.asarray(value), list(value))
    assert type(point + vector) is type(point) and type(point - point) is type(vector)


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_a_non_finite_point_fails_every_fit(bad: float) -> None:
    # #307: one bad scan row built a NaN curve, a fit reported proven.
    from geomcore import ConstructionError

    row = numpy.array([[0.0, 0, 0], [1, 1, 0], [2, bad, 0], [3, 1, 0], [4, 0, 0], [5, 1, 0]])
    fits = [
        lambda: BSplineCurve.interpolate(row, 3),
        lambda: BSplineCurve.approximate(row, 3, 4),
        lambda: BSplineCurve.approximate_within(row, 3, 0.01),
        lambda: BSplineCurve2D.approximate_within(row[:, :2], 3, 0.01),
        lambda: BSplineSurface.approximate_within(numpy.stack([row, row + [0, 0, 1]]), 1, 3, 0.01),
        lambda: BSplineCurve.clamped_uniform(3, [Point(*p) for p in row]),
    ]
    for fit in fits:
        with pytest.raises(ConstructionError, match="NaN or infinite"):
            fit()


def test_a_discretization_converts_to_its_points_and_results_reverse() -> None:
    # #319: numpy saw a discretization as one object, and `reversed` raised.
    from geomcore import Line, Line2D, Plane

    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    planar = Circle2D(Point2D.origin(), 2.0)
    for discretization in (circle.discretize(0.05), planar.discretize(0.05)):
        array = numpy.array(discretization)
        assert array.shape == (len(discretization), len(discretization[0]))
        assert array.dtype == numpy.float64
        numpy.testing.assert_array_equal(array, _rows(discretization.points))
        assert numpy.asarray(discretization, dtype=numpy.float32).dtype == numpy.float32
        assert list(reversed(discretization)) == discretization.points[::-1]
    sphere = Sphere(Point.origin(), 2.0)
    line = Line(Point(-3.0, 0.5, 0.0), Vector.unit_x())
    results = [
        circle.intersect(line),
        planar.intersect(Line2D(Point2D(-3.0, 0.5), Vector2D(1.0, 0.0))),
        line.intersect(sphere),
        Plane(Point(0.0, 0.0, 0.5), Vector.unit_z()).intersect(sphere),
        circle.extrema(Circle(Point(5.0, 0.0, 0.0), Vector.unit_z(), 1.0)),
    ]
    for result in results:
        assert len(result) >= 1
        assert list(reversed(result)) == list(result)[::-1]
