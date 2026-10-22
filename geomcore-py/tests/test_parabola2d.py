"""`Parabola2D` mirrors the 3D parabola in the plane."""

import math

import pytest

from geomcore import Frame2D, Point2D, Vector2D
from geomcore.curves import Line2D, Parabola2D, TrimmedCurve2D

ORIGIN = Point2D.origin()
X, Y = Vector2D.unit_x(), Vector2D.unit_y()


def parabola():
    return Parabola2D(ORIGIN, X, 1.0)


def xy(value):
    return (value.x, value.y)


def test_construction_and_accessors():
    p = parabola()
    assert xy(p.apex()) == (0.0, 0.0)
    assert xy(p.frame().x_direction()) == (1.0, 0.0)
    assert xy(p.frame().y_direction()) == (0.0, 1.0)
    assert p.focal() == 1.0
    assert p.bounds() == (-math.inf, math.inf)
    assert p.period() is None
    assert not p.is_periodic()
    assert not p.is_closed()

    from_frame = Parabola2D.from_frame(Frame2D(Point2D(1.0, 2.0), X, Y), 2.5)
    assert xy(from_frame.apex()) == (1.0, 2.0)
    assert from_frame.focal() == 2.5


def test_construction_validation():
    with pytest.raises(ValueError):
        Parabola2D(ORIGIN, X, 0.0)
    with pytest.raises(ValueError):
        Parabola2D(ORIGIN, X, -1.0)
    with pytest.raises(ValueError):
        Parabola2D(ORIGIN, Vector2D(0.0, 0.0), 1.0)
    with pytest.raises(ValueError):
        Parabola2D(Point2D(math.nan, 0.0), X, 1.0)
    with pytest.raises(ValueError):
        Parabola2D.from_frame(Frame2D(ORIGIN, X, Y), math.inf)


def test_evaluation_derivatives_and_parameter_inversion():
    p = parabola()
    for u in [0.0, 0.7, 2.3, -1.2]:
        point = p.eval_point(u)
        assert xy(point) == pytest.approx((u * u / 4.0, u))
        assert xy(p.eval_derivative(u, 1)) == pytest.approx((u / 2.0, 1.0))
        assert xy(p.eval_derivative(u, 2)) == pytest.approx((0.5, 0.0))
        assert xy(p.eval_derivative(u, 3)) == (0.0, 0.0)
        d1, v1 = p.eval_d1(u)
        assert xy(d1) == xy(point)
        assert xy(v1) == xy(p.eval_derivative(u, 1))
        d2, v2, a2 = p.eval_d2(u)
        assert xy(d2) == xy(point)
        assert xy(v2) == xy(v1)
        assert xy(a2) == xy(p.eval_derivative(u, 2))
        assert p.parameter_of(point) == pytest.approx(u)
    points = p.eval_points([0.0, 2.0])
    assert xy(points[0]) == (0.0, 0.0)
    assert xy(points[1]) == pytest.approx((1.0, 2.0))
    with pytest.raises(ValueError):
        p.eval_derivative(0.0, 0)


def test_contains_and_projection():
    p = parabola()
    assert p.contains(Point2D(1.0, 2.0))
    assert not p.contains(Point2D(0.0, 2.0))

    projection = p.project_point(Point2D(-1.0, 0.0))
    assert projection.parameter == pytest.approx(0.0)
    assert xy(projection.point) == pytest.approx((0.0, 0.0))
    assert projection.distance == pytest.approx(1.0)

    batch = p.project_points([Point2D(-1.0, 0.0), Point2D(1.0, 2.0)])
    assert len(batch) == 2
    assert batch[0].distance == pytest.approx(1.0)
    assert batch[1].parameter == pytest.approx(2.0)
    assert batch[1].distance == 0.0

    # Beyond the centre of curvature at the apex, three feet exist.
    extrema = p.extrema(Point2D(5.0, 0.0))
    assert len(extrema) == 3
    distances = [item.distance for item in extrema]
    assert distances == sorted(distances)
    assert extrema[-1].parameter == pytest.approx(0.0)


def test_local_properties():
    p = parabola()
    assert xy(p.tangent(0.0)) == pytest.approx((0.0, 1.0))
    # At the apex the radius of curvature is twice the focal distance, and
    # the curve turns right as it travels along +y.
    assert p.curvature(0.0) == pytest.approx(-0.5)
    assert xy(p.curvature_center(0.0)) == pytest.approx((2.0, 0.0))
    assert p.continuity(0.0) == math.inf
    assert p.min_continuity() == math.inf
    assert p.bounding_box() is None


def test_trimmed_curve():
    p = parabola()
    arc = TrimmedCurve2D(p, -2.0, 2.0)
    assert arc.bounds() == (-2.0, 2.0)
    assert xy(arc.start_point()) == pytest.approx((1.0, -2.0))
    assert xy(arc.end_point()) == pytest.approx((1.0, 2.0))
    assert arc.contains(Point2D(0.25, 1.0))
    assert not arc.contains(Point2D(4.0, 4.0))
    projection = arc.project_point(Point2D(4.0, 4.0))
    assert projection.parameter == pytest.approx(2.0)
    assert isinstance(arc.basis(), Parabola2D)

    box = arc.bounding_box()
    assert xy(box.min()) == pytest.approx((0.0, -2.0))
    assert xy(box.max()) == pytest.approx((1.0, 2.0))

    criteria = (1e-3, 0.2)
    d = arc.discretize(*criteria)
    assert len(d.points) == arc.point_count(*criteria)
    assert max(d.chord_deviations) <= 1e-3
    assert max(d.angle_deviations) <= 0.2


def test_intersect_and_repr():
    p = parabola()
    hits = p.intersect(Line2D(Point2D(1.0, 0.0), Y))
    assert sorted(hit.point.y for hit in hits.hits) == pytest.approx([-2.0, 2.0])
    assert repr(p) == "Parabola2D(apex=(0.0, 0.0), x_direction=(1.0, 0.0), focal=1.0)"
