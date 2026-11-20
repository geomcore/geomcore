"""The offset of a 2D curve: what the binding adds, the result class and its
curve of each kind, the cusps, and the errors. Which curves offset to what,
and within which deviation, is asserted in the core tests."""

import math

import pytest

import geomcore
from geomcore import (
    BSplineCurve2D,
    Circle2D,
    CurveOffset2D,
    CurveOffset2DError,
    Ellipse2D,
    Line2D,
    Parabola2D,
    Point2D,
    TrimmedCurve2D,
    Vector2D,
)


def test_a_line_offsets_exactly_to_its_left():
    offset = Line2D(Point2D(0, 0), Vector2D(1, 0)).offset(2.0, 1e-6)
    assert isinstance(offset, CurveOffset2D)
    assert offset.is_exact() and offset.is_proven()
    assert offset.deviation == 0.0
    assert offset.cusp_parameters == []
    assert isinstance(offset.curve, Line2D)
    assert offset.curve.point_at(3.0) == Point2D(3.0, 2.0)


def test_a_circle_offsets_to_a_circle_and_refuses_a_radius_past_zero():
    circle = Circle2D(Point2D(0, 0), 3.0)
    offset = circle.offset(1.0, 1e-6, geomcore.Tolerance.default())
    assert isinstance(offset.curve, Circle2D)
    assert offset.curve.radius() == 2.0
    with pytest.raises(CurveOffset2DError, match="radius"):
        circle.offset(3.0, 1e-6)


def test_an_ellipse_offsets_to_a_bspline_within_the_deviation_with_its_cusps():
    ellipse = Ellipse2D(Point2D(0, 0), Vector2D(1, 0), 3.0, 1.0)
    offset = ellipse.offset(0.5, 1e-6)
    assert isinstance(offset.curve, BSplineCurve2D)
    assert not offset.is_exact() and not offset.is_proven()
    assert offset.deviation <= 1e-6
    assert len(offset.cusp_parameters) == 4
    for cusp in offset.cusp_parameters:
        assert 0.5 * ellipse.curvature(cusp) == pytest.approx(1.0, abs=1e-9)
    t = 1.0
    point, velocity = ellipse.point_at(t), ellipse.derivative(t, 1)
    speed = math.hypot(velocity.x, velocity.y)
    expected = Point2D(point.x - 0.5 * velocity.y / speed, point.y + 0.5 * velocity.x / speed)
    assert offset.curve.point_at(t).distance(expected) <= 1e-6
    assert repr(offset) == (
        f"CurveOffset2D(exact=False, deviation={offset.deviation!r}, proven=False, cusp_count=4)"
    )


def test_a_trimmed_curve_offsets_to_a_trimmed_curve_over_its_range():
    segment = TrimmedCurve2D.segment(Point2D(0, 0), Point2D(2, 0))
    offset = segment.offset(-1.0, 1e-6)
    assert isinstance(offset.curve, TrimmedCurve2D)
    assert offset.curve.bounds() == segment.bounds()
    assert offset.curve.start_point() == Point2D(0, -1)
    piece = TrimmedCurve2D(Parabola2D(Point2D(0, 0), Vector2D(1, 0), 1.0), -1.0, 2.0)
    offset = piece.offset(0.2, 1e-6)
    assert isinstance(offset.curve, TrimmedCurve2D)
    assert isinstance(offset.curve.basis(), BSplineCurve2D)


def test_refusals_raise_curve_offset_2d_error():
    parabola = Parabola2D(Point2D(0, 0), Vector2D(1, 0), 1.0)
    with pytest.raises(CurveOffset2DError, match="unbounded"):
        parabola.offset(0.1, 1e-6)
    line = Line2D(Point2D(0, 0), Vector2D(1, 0))
    with pytest.raises(CurveOffset2DError, match="deviation"):
        line.offset(0.1, 0.0)
    with pytest.raises(CurveOffset2DError, match="distance"):
        line.offset(math.nan, 1e-6)
    assert issubclass(CurveOffset2DError, geomcore.GeomcoreError)
