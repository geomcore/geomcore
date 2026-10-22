"""`Ellipse2D` mirrors the 3D ellipse in the plane."""

import math

import pytest

from geomcore import Frame2D, Point2D, Vector2D
from geomcore.curves import BSplineCurve2D, Ellipse2D, TrimmedCurve2D

ORIGIN = Point2D.origin()
X, Y = Vector2D.unit_x(), Vector2D.unit_y()


def ellipse():
    return Ellipse2D(ORIGIN, X, 3.0, 1.5)


def xy(value):
    return (value.x, value.y)


def test_construction_and_accessors():
    e = ellipse()
    assert xy(e.center()) == (0.0, 0.0)
    assert xy(e.frame().x_direction()) == (1.0, 0.0)
    assert xy(e.frame().y_direction()) == (0.0, 1.0)
    assert e.major_radius() == 3.0
    assert e.minor_radius() == 1.5
    assert e.area() == pytest.approx(math.pi * 4.5)
    assert e.bounds() == (0.0, math.tau)
    assert e.period() == math.tau
    assert e.is_periodic()
    assert e.is_closed()

    from_frame = Ellipse2D.from_frame(Frame2D(Point2D(1.0, 2.0), X, Y), 2.0, 1.0)
    assert xy(from_frame.center()) == (1.0, 2.0)

    from_points = Ellipse2D.from_center_and_points(
        ORIGIN, Point2D(3.0, 0.0), Point2D(0.0, 1.5)
    )
    assert from_points.major_radius() == pytest.approx(3.0)
    assert from_points.minor_radius() == pytest.approx(1.5)


def test_construction_validation():
    with pytest.raises(ValueError):
        Ellipse2D(ORIGIN, X, 3.0, 0.0)
    with pytest.raises(ValueError):
        Ellipse2D(ORIGIN, X, 1.0, 3.0)
    with pytest.raises(ValueError):
        Ellipse2D(ORIGIN, Vector2D(0.0, 0.0), 3.0, 1.0)
    with pytest.raises(ValueError):
        Ellipse2D(Point2D(math.nan, 0.0), X, 3.0, 1.0)
    with pytest.raises(ValueError):
        Ellipse2D.from_center_and_points(ORIGIN, ORIGIN, Point2D(0.0, 1.0))
    with pytest.raises(ValueError):
        Ellipse2D.from_center_and_points(
            ORIGIN, Point2D(3.0, 0.0), Point2D(6.0, 0.0)
        )


def test_evaluation_derivatives_and_parameter_inversion():
    e = ellipse()
    for u in [0.0, 0.7, 2.3, 5.9, -1.2]:
        p = e.eval_point(u)
        assert xy(p) == pytest.approx((3.0 * math.cos(u), 1.5 * math.sin(u)))
        assert xy(e.eval_derivative(u, 1)) == pytest.approx(
            (-3.0 * math.sin(u), 1.5 * math.cos(u))
        )
        assert xy(e.eval_derivative(u, 2)) == pytest.approx(
            (-3.0 * math.cos(u), -1.5 * math.sin(u))
        )
        d1, v1 = e.eval_d1(u)
        assert xy(d1) == xy(p)
        assert xy(v1) == xy(e.eval_derivative(u, 1))
        d2, v2, a2 = e.eval_d2(u)
        assert xy(d2) == xy(p)
        assert xy(v2) == xy(v1)
        assert xy(a2) == xy(e.eval_derivative(u, 2))
        assert e.parameter_of(p) == pytest.approx(u % math.tau)
    points = e.eval_points([0.0, math.pi / 2])
    assert xy(points[0]) == pytest.approx((3.0, 0.0))
    assert xy(points[1]) == pytest.approx((0.0, 1.5))


def test_contains_and_projection():
    e = ellipse()
    assert e.contains(Point2D(3.0, 0.0))
    assert e.contains(Point2D(0.0, 1.5))
    assert not e.contains(Point2D(4.0, 0.0))

    projection = e.project_point(Point2D(4.0, 0.0))
    assert projection.parameter == pytest.approx(0.0)
    assert xy(projection.point) == pytest.approx((3.0, 0.0))
    assert projection.distance == pytest.approx(1.0)

    batch = e.project_points([Point2D(4.0, 0.0), Point2D(0.0, 3.0)])
    assert len(batch) == 2
    assert batch[0].distance == pytest.approx(1.0)
    assert batch[1].parameter == pytest.approx(math.pi / 2)

    extrema = e.extrema(Point2D(4.0, 0.0))
    assert extrema
    assert extrema[0].distance == pytest.approx(1.0)
    distances = [item.distance for item in extrema]
    assert distances == sorted(distances)


def test_local_properties():
    e = ellipse()
    assert xy(e.tangent(0.0)) == pytest.approx((0.0, 1.0))
    assert e.curvature(0.0) == pytest.approx(4.0 / 3.0)
    assert xy(e.curvature_center(0.0)) == pytest.approx((2.25, 0.0))
    assert e.continuity(0.0) == math.inf
    assert e.min_continuity() == math.inf


def test_bounding_box_and_discretization():
    e = ellipse()
    box = e.bounding_box()
    assert xy(box.min()) == pytest.approx((-3.0, -1.5))
    assert xy(box.max()) == pytest.approx((3.0, 1.5))

    criteria = (1e-3, 0.2)
    count = e.point_count(*criteria)
    d = e.discretize(*criteria)
    assert len(d.points) == count
    assert d.points[0] == d.points[-1]
    assert max(d.chord_deviations) <= 1e-3
    assert max(d.angle_deviations) <= 0.2


def test_trimmed_curve():
    e = ellipse()
    arc = TrimmedCurve2D(e, 0.0, math.pi / 2)
    assert arc.bounds() == (0.0, math.pi / 2)
    assert xy(arc.start_point()) == pytest.approx((3.0, 0.0))
    assert xy(arc.end_point()) == pytest.approx((0.0, 1.5))
    assert xy(arc.eval_point(0.0)) == pytest.approx((3.0, 0.0))
    assert arc.contains(Point2D(3.0, 0.0))
    projection = arc.project_point(Point2D(0.0, 3.0))
    assert projection.parameter == pytest.approx(math.pi / 2)
    bspline = arc.to_bspline()
    assert isinstance(bspline, BSplineCurve2D)
    assert bspline.bounds() == (0.0, math.pi / 2)


def test_to_bspline_and_repr():
    e = ellipse()
    bspline = e.to_bspline()
    assert bspline.is_periodic()
    assert bspline.is_rational()
    assert bspline.bounds() == (0.0, math.tau)
    assert xy(bspline.eval_point(math.pi / 2)) == pytest.approx((0.0, 1.5))
    assert repr(e) == (
        "Ellipse2D(center=(0.0, 0.0), x_direction=(1.0, 0.0), "
        "major_radius=3.0, minor_radius=1.5)"
    )
