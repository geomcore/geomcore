"""Curve–curve extrema between 2D curves (#422)."""

import math

import pytest

import geomcore
from geomcore import (
    BSplineCurve2D,
    Circle,
    Circle2D,
    CurveCurveExtrema2D,
    CurveCurveExtremum2D,
    CurveProjection2D,
    Ellipse2D,
    Line2D,
    Point,
    Point2D,
    TrimmedCurve2D,
    Vector,
    Vector2D,
)


def test_two_circles_answer_their_gap_with_both_points():
    first = Circle2D(Point2D(0.0, 0.0), 1.0)
    second = Circle2D(Point2D(3.0, 4.0), 2.0)
    extrema = first.extrema(second)
    assert isinstance(extrema, CurveCurveExtrema2D)
    assert extrema.first is first and extrema.second is second
    assert len(extrema) == 1
    closest = extrema.closest()
    assert isinstance(closest, CurveCurveExtremum2D)
    assert closest.distance == pytest.approx(2.0)
    assert closest.first_point.distance(Point2D(0.6, 0.8)) < 1e-7
    assert closest.second_point.distance(Point2D(1.8, 2.4)) < 1e-7
    assert first.minimum_distance(second) == closest


def test_a_point_still_answers_its_projections():
    circle = Circle2D(Point2D(0.0, 0.0), 1.0)
    (projection,) = circle.extrema((3.0, 0.0))
    assert isinstance(projection, CurveProjection2D)
    assert projection.distance == pytest.approx(2.0)


def test_concentric_circles_report_their_gap_once():
    outer = Circle2D(Point2D(1.0, 2.0), 3.0)
    inner = Circle2D(Point2D(1.0, 2.0), 1.0)
    (minimum,) = outer.extrema(inner).minima
    assert minimum.distance == pytest.approx(2.0)
    assert minimum.first_parameter == 0.0


def test_an_ellipse_is_its_semi_minor_end_from_a_parallel_line():
    ellipse = Ellipse2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0), 3.0, 1.0)
    line = Line2D(Point2D(0.0, 4.0), Vector2D(1.0, 0.0))
    assert ellipse.minimum_distance(line).distance == pytest.approx(3.0)


def test_a_segment_and_an_arc_are_their_closest_ends_apart():
    segment = TrimmedCurve2D.segment(Point2D(-1.0, 3.0), Point2D(1.0, 3.0))
    arc = TrimmedCurve2D(Circle2D(Point2D(0.0, 0.0), 1.0), 0.0, 0.5 * math.pi)
    closest = segment.extrema(arc).closest()
    assert closest.distance == pytest.approx(2.0)
    assert closest.second_point.distance(Point2D(0.0, 1.0)) < 1e-7


def test_two_bsplines_are_searched_to_completeness():
    arch = BSplineCurve2D(
        2, [Point2D(0, 0), Point2D(1, 2), Point2D(2, 0)], [0.0, 1.0], [3, 3], False
    )
    rail = BSplineCurve2D(1, [Point2D(-1, 2), Point2D(3, 2)], [0.0, 1.0], [2, 2], False)
    extrema = arch.extrema(rail)
    assert extrema.is_complete()
    assert extrema[0].distance == pytest.approx(1.0)


def test_the_deprecated_function_takes_2d_pairs():
    first = Circle2D(Point2D(0.0, 0.0), 1.0)
    second = Circle2D(Point2D(4.0, 0.0), 1.0)
    with pytest.warns(DeprecationWarning):
        old = geomcore.curve_curve_extrema(first, second)  # type: ignore[deprecated]
    assert old == first.extrema(second)


def test_a_2d_curve_and_a_3d_curve_have_no_extrema():
    plane_circle = Circle2D(Point2D(0.0, 0.0), 1.0)
    space_circle = Circle(Point(0.0, 0.0, 0.0), Vector(0.0, 0.0, 1.0), 1.0)
    with pytest.raises(TypeError, match="no extrema with a 3D curve"):
        plane_circle.extrema(space_circle)  # type: ignore[call-overload]
    with pytest.raises(TypeError, match="not a 2D curve"):
        space_circle.extrema(plane_circle)  # type: ignore[call-overload]
