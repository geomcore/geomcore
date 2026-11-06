"""Building a B-spline the textbook way and reading its shape (#184): flat
knot vectors, clamped uniform knots, single-pole and single-weight edits,
geometric continuity, and a repr that shows the knots of a small curve."""

import math

import pytest

from geomcore import BSplineCurve, BSplineCurve2D, ConstructionError, Point, Point2D, Tolerance

_POLES = [Point(0, 0, 0), Point(1, 2, 0.5), Point(3, 2, 1), Point(4, 0, 0.3), Point(5, 1, 2)]


def test_the_nurbs_guide_example_runs():
    poles = [Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)]
    curve = BSplineCurve(2, poles, [0.0, 1.0], [3, 3], False)
    flat = BSplineCurve.from_flat_knots(2, poles, [0, 0, 0, 1, 1, 1])
    assert flat == curve and curve.flat_knots() == [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    assert BSplineCurve.clamped_uniform(2, poles) == curve
    arc = curve.with_weight(1, 2**-0.5)
    assert arc.point_at(0.5).distance(Point(0, 0, 0)) == pytest.approx(1.0, abs=1e-15)
    assert curve.with_pole(-1, Point(0, 2, 0)).poles()[2] == Point(0, 2, 0)


def test_flat_knots_round_trip_rational_and_periodic_curves():
    rational = BSplineCurve.new_rational(2, _POLES, [1, 0.5, 2, 1, 0.7], [0, 1, 2, 3], [3, 1, 1, 3], False)
    periodic = BSplineCurve2D(2, [Point2D(p.x, p.y) for p in _POLES], [0, 1, 2, 3, 4, 5], [1] * 6, True)
    for curve in (rational, periodic):
        rebuilt = type(curve).from_flat_knots(
            curve.degree(), curve.poles(), curve.flat_knots(), curve.weights(), curve.is_periodic()
        )
        assert rebuilt == curve


def test_a_decreasing_flat_knot_vector_is_refused():
    with pytest.raises(ConstructionError):
        BSplineCurve.from_flat_knots(3, _POLES, [0, 0, 0, 0, 0.6, 0.4, 1, 1, 1, 1])


def test_clamped_uniform_takes_an_array_of_poles():
    numpy = pytest.importorskip("numpy")
    curve = BSplineCurve2D.clamped_uniform(3, numpy.array([[0, 0], [1, 2], [3, 2], [4, 0], [5, 1]]))
    assert curve.knots() == [0.0, 0.5, 1.0]
    assert curve.multiplicities() == [4, 1, 4]
    with pytest.raises(ConstructionError):
        BSplineCurve.clamped_uniform(5, _POLES)


def test_with_pole_and_with_weight_index_like_a_list():
    curve = BSplineCurve.clamped_uniform(3, _POLES)
    assert curve.with_pole(-1, Point(9, 9, 9)).poles()[4] == Point(9, 9, 9)
    assert curve.poles() == _POLES
    assert curve.with_weight(-5, 3.0).weights() == [3.0, 1.0, 1.0, 1.0, 1.0]
    for index in (5, -6):
        with pytest.raises(IndexError):
            curve.with_pole(index, Point(0, 0, 0))
        with pytest.raises(IndexError):
            curve.with_weight(index, 1.0)
    for weight in (0.0, -2.0, math.nan, math.inf):
        with pytest.raises(ConstructionError):
            curve.with_weight(0, weight)


def test_geometric_continuity_sees_a_smooth_shape_through_a_full_knot():
    curve = BSplineCurve.clamped_uniform(3, _POLES)
    split = curve.insert_knot(0.25, 3)
    assert split.continuity(0.25) == 0
    assert split.geometric_continuity(0.25) == 2
    assert split.geometric_continuity(0.1) == math.inf
    kinked = BSplineCurve2D(1, [Point2D(0, 0), Point2D(1, 0), Point2D(1, 1)], [0, 1, 2], [2, 1, 2], False)
    assert kinked.geometric_continuity(1.0) == 0
    # An angular tolerance past the right angle takes the kink for
    # tangent, and two straight pieces agree in curvature.
    assert kinked.geometric_continuity(1.0, Tolerance(angular=2.0)) == 2


def test_the_repr_shows_the_knots_of_a_small_curve():
    curve = BSplineCurve.clamped_uniform(2, _POLES)
    assert repr(curve) == (
        "BSplineCurve(degree=2, pole_count=5, rational=False, periodic=False, "
        "knots=[0.0, 0.3333333333333333, 0.6666666666666666, 1.0], multiplicities=[3, 1, 1, 3])"
    )
    many = BSplineCurve.clamped_uniform(1, [Point(i, 0, 0) for i in range(12)])
    assert repr(many) == "BSplineCurve(degree=1, pole_count=12, rational=False, periodic=False)"
