"""2D curves gain the local properties, closedness, continuity and
bounding boxes of their 3D counterparts; the plane's curvature is signed."""

import math

import pytest

from geomcore import BoundingBox2D, Point2D, Tolerance, Vector2D
from geomcore.curves import BSplineCurve2D, Circle2D, Line2D, TrimmedCurve2D

ORIGIN = Point2D.origin()
X, Y = Vector2D.unit_x(), Vector2D.unit_y()


def xy(value):
    return (value.x, value.y)


def arch():
    """An arch through (0, 0) and (2, 0) that turns right."""
    poles = [Point2D(0.0, 0.0), Point2D(1.0, 1.0), Point2D(2.0, 0.0)]
    return BSplineCurve2D(2, poles, [0.0, 1.0], [3, 3], False)


def kinked():
    poles = [Point2D(0.0, 0.0), Point2D(1.0, 0.0), Point2D(2.0, 1.0)]
    poles += [Point2D(3.0, 1.0), Point2D(4.0, 0.0)]
    return BSplineCurve2D(2, poles, [0.0, 0.5, 1.0], [3, 2, 3], False)


def test_line_is_straight():
    line = Line2D(ORIGIN, X)
    assert xy(line.tangent(5.0)) == (1.0, 0.0)
    assert line.curvature(5.0) == 0.0
    assert line.curvature_center(5.0) is None


def test_counterclockwise_circle_has_positive_curvature_about_its_center():
    circle = Circle2D(Point2D(1.0, 1.0), 2.0)
    assert xy(circle.tangent(0.0)) == (0.0, 1.0)
    assert circle.curvature(0.0) == 0.5
    assert xy(circle.curvature_center(0.0)) == pytest.approx((1.0, 1.0), abs=1e-12)


def test_reversing_a_curve_flips_the_sign_of_its_curvature():
    curve = arch()
    reversed_curve = curve.reversed()
    assert xy(reversed_curve.tangent(0.5)) == (-1.0, 0.0)
    assert reversed_curve.curvature(0.5) == pytest.approx(1.0)
    # The osculating circle is the same one, now on the left.
    assert xy(reversed_curve.curvature_center(0.5)) == pytest.approx((1.0, -0.5))


def test_bspline_arch_turns_right():
    curve = arch()
    assert curve.curvature(0.5) == pytest.approx(-1.0)
    assert xy(curve.curvature_center(0.5)) == pytest.approx((1.0, -0.5))


def test_stopped_curve_has_no_local_properties():
    poles = [Point2D(1.0, 1.0), Point2D(1.0, 1.0)]
    stopped = BSplineCurve2D(1, poles, [0.0, 1.0], [2, 2], False)
    assert stopped.tangent(0.5) is None
    assert stopped.curvature(0.5) is None
    assert stopped.curvature_center(0.5) is None


def test_tolerance_is_accepted():
    circle = Circle2D(ORIGIN, 1.0)
    tol = Tolerance.default()
    assert circle.curvature(0.3, tol) == circle.curvature(0.3)


def test_trimmed_arc_has_its_basis_curvature():
    arc = TrimmedCurve2D(Circle2D(ORIGIN, 2.0), 0.0, 1.0)
    assert arc.curvature(0.5) == 0.5
    assert xy(arc.tangent(0.0)) == (0.0, 1.0)
    assert xy(arc.curvature_center(0.5)) == pytest.approx((0.0, 0.0), abs=1e-12)


def test_closedness():
    assert not Line2D(ORIGIN, X).is_closed()
    assert Circle2D(ORIGIN, 1.0).is_closed()
    circle = Circle2D(ORIGIN, 1.0)
    assert TrimmedCurve2D(circle, 0.0, math.tau).is_closed()
    assert not TrimmedCurve2D(circle, 0.0, 1.0).is_closed()
    assert not TrimmedCurve2D.segment(ORIGIN, Point2D(1.0, 0.0)).is_closed()


def test_continuity():
    assert Circle2D(ORIGIN, 1.0).continuity(0.5) == math.inf
    assert Line2D(ORIGIN, X).min_continuity() == math.inf
    curve = kinked()
    piece = TrimmedCurve2D(curve, 0.25, 0.75)
    assert piece.continuity(0.5) == 0
    assert piece.continuity(0.25) == math.inf
    assert piece.min_continuity() == 0
    assert TrimmedCurve2D(curve, 0.0, 0.5).min_continuity() == math.inf


def test_bounding_boxes():
    assert Line2D(ORIGIN, X).bounding_box() is None
    box = Circle2D(Point2D(1.0, 2.0), 2.0).bounding_box()
    assert xy(box.min()) == (-1.0, 0.0)
    assert xy(box.max()) == (3.0, 4.0)

    quarter = TrimmedCurve2D(Circle2D(Point2D(1.0, 2.0), 2.0), 0.0, math.pi / 2)
    box = quarter.bounding_box()
    assert xy(box.min()) == pytest.approx((1.0, 2.0), abs=1e-12)
    assert xy(box.max()) == pytest.approx((3.0, 4.0), abs=1e-12)

    box = arch().bounding_box()
    assert xy(box.min()) == (0.0, 0.0)
    assert xy(box.max()) == (2.0, 1.0)


def test_bounding_box2d_value():
    a = BoundingBox2D(Point2D(0.0, 0.0), Point2D(2.0, 1.0))
    b = BoundingBox2D(Point2D(3.0, 0.0), Point2D(4.0, 1.0))
    assert not a.overlaps(b)
    assert a.enlarged(0.5).overlaps(b.enlarged(0.5))
    assert xy(a.union(b).max()) == (4.0, 1.0)
    assert a.encloses(Point2D(1.0, 0.5))
    assert xy(a.center()) == (1.0, 0.5)
    assert xy(a.diagonal()) == (2.0, 1.0)
    assert BoundingBox2D.from_points([Point2D(1.0, 0.0), Point2D(-1.0, 2.0)]) == BoundingBox2D(
        Point2D(-1.0, 0.0), Point2D(1.0, 2.0)
    )
    with pytest.raises(ValueError):
        BoundingBox2D.from_points([])
    assert repr(a) == "BoundingBox2D(min=(0.0, 0.0), max=(2.0, 1.0))"
