"""`Hyperbola2D` mirrors the 3D hyperbola in the plane."""

import math

import pytest

from geomcore import Frame2D, Point2D, Vector2D
from geomcore.curves import Hyperbola2D, Line2D, TrimmedCurve2D

ORIGIN = Point2D.origin()
X, Y = Vector2D.unit_x(), Vector2D.unit_y()


def hyperbola():
    return Hyperbola2D(ORIGIN, X, 2.0, 1.0)


def xy(value):
    return (value.x, value.y)


def test_construction_and_accessors():
    h = hyperbola()
    assert xy(h.center()) == (0.0, 0.0)
    assert xy(h.frame().x_direction()) == (1.0, 0.0)
    assert xy(h.frame().y_direction()) == (0.0, 1.0)
    assert h.major_radius() == 2.0
    assert h.minor_radius() == 1.0
    assert h.bounds() == (-math.inf, math.inf)
    assert h.period() is None
    assert not h.is_periodic()
    assert not h.is_closed()

    from_frame = Hyperbola2D.from_frame(Frame2D(Point2D(1.0, 2.0), X, Y), 1.0, 2.5)
    assert xy(from_frame.center()) == (1.0, 2.0)
    assert from_frame.major_radius() == 1.0
    assert from_frame.minor_radius() == 2.5

    from_points = Hyperbola2D.from_center_and_points(
        Point2D(1.0, 1.0), Point2D(3.0, 1.0), Point2D(1.0, -2.0)
    )
    assert xy(from_points.center()) == (1.0, 1.0)
    assert from_points.major_radius() == pytest.approx(2.0)
    assert from_points.minor_radius() == pytest.approx(3.0)
    # The y direction points toward the second point.
    assert xy(from_points.frame().y_direction()) == pytest.approx((0.0, -1.0))


def test_construction_validation():
    with pytest.raises(ValueError):
        Hyperbola2D(ORIGIN, X, 0.0, 1.0)
    with pytest.raises(ValueError):
        Hyperbola2D(ORIGIN, X, 2.0, -1.0)
    with pytest.raises(ValueError):
        Hyperbola2D(ORIGIN, Vector2D(0.0, 0.0), 2.0, 1.0)
    with pytest.raises(ValueError):
        Hyperbola2D(Point2D(math.nan, 0.0), X, 2.0, 1.0)
    with pytest.raises(ValueError):
        Hyperbola2D.from_frame(Frame2D(ORIGIN, X, Y), math.inf, 1.0)
    with pytest.raises(ValueError):
        Hyperbola2D.from_center_and_points(ORIGIN, ORIGIN, Point2D(0.0, 1.0))
    with pytest.raises(ValueError):
        Hyperbola2D.from_center_and_points(
            ORIGIN, Point2D(2.0, 0.0), Point2D(5.0, 0.0)
        )


def test_evaluation_derivatives_and_parameter_inversion():
    h = hyperbola()
    for u in [0.0, 0.7, 2.3, -1.2]:
        point = h.eval_point(u)
        assert xy(point) == pytest.approx((2.0 * math.cosh(u), math.sinh(u)))
        assert xy(h.eval_derivative(u, 1)) == pytest.approx(
            (2.0 * math.sinh(u), math.cosh(u))
        )
        assert xy(h.eval_derivative(u, 2)) == pytest.approx(
            (2.0 * math.cosh(u), math.sinh(u))
        )
        assert xy(h.eval_derivative(u, 3)) == pytest.approx(
            (2.0 * math.sinh(u), math.cosh(u))
        )
        d1, v1 = h.eval_d1(u)
        assert xy(d1) == xy(point)
        assert xy(v1) == xy(h.eval_derivative(u, 1))
        d2, v2, a2 = h.eval_d2(u)
        assert xy(d2) == xy(point)
        assert xy(v2) == xy(v1)
        assert xy(a2) == xy(h.eval_derivative(u, 2))
        assert h.parameter_of(point) == pytest.approx(u)
    points = h.eval_points([0.0, 1.0])
    assert xy(points[0]) == (2.0, 0.0)
    assert xy(points[1]) == pytest.approx((2.0 * math.cosh(1.0), math.sinh(1.0)))
    with pytest.raises(ValueError):
        h.eval_derivative(0.0, 0)


def test_contains_and_projection():
    h = hyperbola()
    assert h.contains(Point2D(2.0, 0.0))
    assert h.contains(h.eval_point(1.5))
    assert not h.contains(Point2D(0.0, 0.0))
    assert not h.contains(Point2D(-2.0, 0.0))

    projection = h.project_point(Point2D(3.0, 0.0))
    assert abs(projection.parameter) == pytest.approx(0.6223625, abs=1e-6)
    assert projection.distance == pytest.approx(math.sqrt(0.8))

    batch = h.project_points([Point2D(3.0, 0.0), h.eval_point(0.7)])
    assert len(batch) == 2
    assert batch[0].distance == pytest.approx(math.sqrt(0.8))
    assert batch[1].parameter == pytest.approx(0.7)
    assert batch[1].distance == 0.0

    extrema = h.extrema(Point2D(3.0, 0.0))
    assert len(extrema) >= 2
    distances = [item.distance for item in extrema]
    assert distances == sorted(distances)


def test_local_properties():
    h = hyperbola()
    assert xy(h.tangent(0.0)) == pytest.approx((0.0, 1.0))
    # At the vertex the radius of curvature is minor^2 / major = 0.5, and the
    # curve turns right as it travels along +y.
    assert h.curvature(0.0) == pytest.approx(-2.0)
    assert xy(h.curvature_center(0.0)) == pytest.approx((2.5, 0.0))
    assert h.continuity(0.0) == math.inf
    assert h.min_continuity() == math.inf
    assert h.bounding_box() is None


def test_trimmed_curve():
    h = hyperbola()
    arc = TrimmedCurve2D(h, -1.0, 1.0)
    assert arc.bounds() == (-1.0, 1.0)
    assert xy(arc.start_point()) == pytest.approx(xy(h.eval_point(-1.0)))
    assert xy(arc.end_point()) == pytest.approx(xy(h.eval_point(1.0)))
    assert arc.contains(h.eval_point(0.25))
    assert not arc.contains(h.eval_point(3.0))
    projection = arc.project_point(h.eval_point(3.0))
    assert projection.parameter == pytest.approx(1.0)
    assert isinstance(arc.basis(), Hyperbola2D)

    box = arc.bounding_box()
    assert xy(box.min()) == pytest.approx((2.0, -math.sinh(1.0)))
    assert xy(box.max()) == pytest.approx((2.0 * math.cosh(1.0), math.sinh(1.0)))

    criteria = (1e-3, 0.2)
    d = arc.discretize(*criteria)
    assert len(d.points) == arc.point_count(*criteria)
    assert max(d.chord_deviations) <= 1e-3
    assert max(d.angle_deviations) <= 0.2


def test_intersect_and_repr():
    h = hyperbola()
    hits = h.intersect(Line2D(Point2D(3.0, 0.0), Y))
    expected = math.sqrt(5.0) / 2.0
    assert sorted(hit.point.y for hit in hits.hits) == pytest.approx(
        [-expected, expected]
    )
    assert (
        repr(h)
        == "Hyperbola2D(center=(0.0, 0.0), x_direction=(1.0, 0.0), major_radius=2.0, minor_radius=1.0)"
    )
