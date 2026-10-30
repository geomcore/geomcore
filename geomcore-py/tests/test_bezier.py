"""Tests for the Bézier curve and surface bindings."""

import math

import pytest

import geomcore
from geomcore import ConstructionError, Point2D, Point, Side, Tolerance, Vector, Vector2D
from geomcore.curves import (
    BSplineCurve,
    BSplineCurve2D,
    BezierCurve,
    BezierCurve2D,
    Line,
    Line2D,
    TrimmedCurve,
)
from geomcore.surfaces import BSplineSurface, BezierSurface, Plane

W = math.sqrt(0.5)


def cubic() -> BezierCurve:
    return BezierCurve(
        3,
        [Point(0.0, 0.0, 0.0), Point(1.0, 2.0, 0.0), Point(3.0, 2.0, 1.0), Point(4.0, 0.0, 1.0)],
    )


def saddle() -> BezierSurface:
    """The surface (2u, 2v, 4uv)."""
    poles = [[Point(float(i), float(j), float(i * j)) for j in range(3)] for i in range(3)]
    return BezierSurface(2, 2, poles)


def quarter_cylinder() -> BezierSurface:
    arc = [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    poles = [[Point(x, y, 0.0), Point(x, y, 1.0)] for x, y in arc]
    weights = [[1.0, 1.0], [W, W], [1.0, 1.0]]
    return BezierSurface.new_rational(2, 1, poles, weights)


def close(a: Point, b: Point, tol: float = 1e-12) -> bool:
    return a.distance(b) < tol


# ---- BezierCurve ----


def test_curve_accessors_and_eval():
    curve = cubic()
    assert (curve.degree(), len(curve.poles()), curve.is_rational()) == (3, 4, False)
    assert curve.weights() is None
    assert curve.bounds() == (0.0, 1.0)
    assert curve.eval_point(0.0) == Point(0.0, 0.0, 0.0)
    assert curve.eval_point(0.5) == Point(2.0, 1.5, 0.5)
    assert curve.eval_derivative(0.5, 1) == Vector(4.5, 0.0, 1.5)
    point, d1, d2 = curve.eval_d2(0.5)
    assert (point, d1) == (Point(2.0, 1.5, 0.5), Vector(4.5, 0.0, 1.5))
    assert d2 == Vector(0.0, -12.0, 0.0)
    assert curve.eval_points([0.5])[0] == point
    assert curve.eval_d1_sided(0.5, Side.BELOW) == curve.eval_d1_sided(0.5, Side.ABOVE)
    assert curve.continuity(0.5) == math.inf
    assert curve.min_continuity() == math.inf
    box = curve.bounding_box()
    assert box is not None and box.max() == Point(4.0, 2.0, 1.0)
    assert repr(curve) == "BezierCurve(degree=3, pole_count=4, rational=False)"


def test_curve_rational_quarter_circle():
    arc = BezierCurve.new_rational(
        2, [Point(1.0, 0.0, 0.0), Point(1.0, 1.0, 0.0), Point(0.0, 1.0, 0.0)], [1.0, W, 1.0]
    )
    assert arc.is_rational() and arc.weights() == [1.0, W, 1.0]
    mid = arc.eval_point(0.5)
    assert mid.x == pytest.approx(W) and mid.y == pytest.approx(W)
    assert arc.length(1e-12) == pytest.approx(math.pi / 2, abs=1e-9)


def test_curve_errors():
    with pytest.raises(ConstructionError):
        BezierCurve(0, [Point.origin()])
    with pytest.raises(ConstructionError):
        BezierCurve(3, [Point.origin()] * 3)
    with pytest.raises(ConstructionError):
        BezierCurve.new_rational(1, [Point.origin()] * 2, [1.0, 0.0])
    with pytest.raises(ValueError):
        cubic().eval_derivative(0.5, 0)


def test_curve_derivatives_of_every_order():
    # The third derivative of a cubic is the constant 6 (P3 - 3 P2 + 3 P1 - P0).
    third = cubic().eval_derivative(0.5, 3)
    assert (third.x, third.y, third.z) == pytest.approx((-12.0, 0.0, -12.0))
    fourth = cubic().eval_derivative(0.5, 4)
    assert (fourth.x, fourth.y, fourth.z) == (0.0, 0.0, 0.0)


def test_curve_editing_keeps_the_shape():
    curve = cubic()
    left, right = curve.split_at(0.25)
    assert (left.bounds(), right.bounds()) == ((0.0, 0.25), (0.25, 1.0))
    assert left.poles()[-1] == right.poles()[0]
    assert close(right.eval_point(0.6), curve.eval_point(0.6))
    assert close(curve.elevate_degree(2).eval_point(0.3), curve.eval_point(0.3))
    assert curve.elevate_degree(2).degree() == 5
    assert curve.reversed().eval_point(0.0) == curve.eval_point(1.0)
    with pytest.raises(geomcore.BSplineSplitError):
        curve.split_at(1.0)


def test_curve_round_trip_with_bspline():
    curve = cubic()
    spline = curve.to_bspline()
    assert isinstance(spline, BSplineCurve)
    assert (spline.knots(), spline.multiplicities()) == ([0.0, 1.0], [4, 4])
    assert BezierCurve.from_bspline(spline) == curve
    with pytest.raises(ConstructionError):
        BezierCurve.from_bspline(spline.insert_knot(0.5, 1))


def test_trimmed_curve_to_bspline_extracts_the_piece():
    curve = cubic()
    piece = TrimmedCurve(curve, 0.25, 0.75).to_bspline()
    assert piece == curve.to_bspline().extract(0.25, 0.75)


def test_curve_projection_and_discretization():
    segment = BezierCurve(1, [Point.origin(), Point(2.0, 0.0, 0.0)])
    projection = segment.project_point(Point(1.0, 1.0, 0.0))
    assert projection.parameter == pytest.approx(0.5)
    assert projection.distance == 1.0
    assert projection.target is segment
    assert segment.parameter_of(Point(1.0, 1.0, 0.0)) == pytest.approx(0.5)
    assert segment.contains(Point(1.0, 0.0, 0.0))
    assert not segment.contains(Point(1.0, 1.0, 0.0))
    assert segment.extrema(Point(1.0, 1.0, 0.0))[0].distance == 1.0
    assert len(segment.project_points([Point(0.0, 1.0, 0.0)])) == 1
    curve = cubic()
    points = curve.discretize(0.01)
    assert len(points.points) == curve.point_count(0.01)
    assert curve.is_closed() is False


# ---- BezierCurve2D ----


def test_curve_2d():
    curve = BezierCurve2D(2, [Point2D(0.0, 0.0), Point2D(1.0, 1.0), Point2D(2.0, 0.0)])
    assert (curve.degree(), curve.bounds()) == (2, (0.0, 1.0))
    assert curve.eval_point(0.5) == Point2D(1.0, 0.5)
    assert curve.eval_d2(0.5)[2].y == pytest.approx(-4.0)
    left, right = curve.split_at(0.5)
    assert left.poles()[-1] == right.poles()[0]
    assert curve.to_bspline().multiplicities() == [3, 3]
    assert BezierCurve2D.from_bspline(curve.to_bspline()) == curve
    assert curve.project_point(Point2D(1.0, 2.0)).target is curve
    box = curve.bounding_box()
    assert box is not None and box.max() == Point2D(2.0, 1.0)
    rational = BezierCurve2D.new_rational(
        2, [Point2D(1.0, 0.0), Point2D(1.0, 1.0), Point2D(0.0, 1.0)], [1.0, W, 1.0]
    )
    assert rational.is_rational()
    with pytest.raises(ConstructionError):
        BezierCurve2D(2, [Point2D(0.0, 0.0)] * 2)
    assert len(curve.discretize(0.01).points) == curve.point_count(0.01)


def test_bspline_2d_bezier_pieces():
    spline = BSplineCurve2D(
        3,
        [Point2D(float(i), float(i % 2)) for i in range(6)],
        [0.0, 1.0, 2.0, 3.0],
        [4, 1, 1, 4],
        False,
    )
    pieces = spline.to_bezier_pieces()
    assert len(pieces) == 3 and all(isinstance(p, BezierCurve2D) for p in pieces)
    assert pieces[1].bounds() == (1.0, 2.0)
    assert pieces[1].eval_point(1.5).distance(spline.eval_point(1.5)) < 1e-12


# ---- BezierSurface ----


def test_surface_accessors_and_eval():
    surface = saddle()
    assert (surface.u_degree(), surface.v_degree()) == (2, 2)
    assert (surface.u_pole_count(), surface.v_pole_count()) == (3, 3)
    assert surface.poles()[1][2] == Point(1.0, 2.0, 2.0)
    assert surface.weights() is None and not surface.is_rational()
    assert (surface.u_bounds(), surface.v_bounds()) == ((0.0, 1.0), (0.0, 1.0))
    assert (surface.u_period(), surface.v_period()) == (None, None)
    assert surface.eval_point(1.0, 1.0) == Point(2.0, 2.0, 4.0)
    assert close(surface.eval_point(0.25, 0.5), Point(0.5, 1.0, 0.5), 1e-14)
    point, (su, sv), (suu, suv, svv) = surface.eval_d2(0.25, 0.5)
    assert close(point, Point(0.5, 1.0, 0.5), 1e-14)
    assert (su - Vector(2.0, 0.0, 2.0)).magnitude() < 1e-12
    assert (sv - Vector(0.0, 2.0, 1.0)).magnitude() < 1e-12
    assert (suv - Vector(0.0, 0.0, 4.0)).magnitude() < 1e-12
    assert suu.magnitude() < 1e-12 and svv.magnitude() < 1e-12
    assert surface.eval_derivative(0.25, 0.5, 1, 1) == suv
    assert surface.eval_points([(0.25, 0.5)])[0] == point
    assert surface.eval_d1_sided(0.3, 0.7, Side.BELOW, Side.ABOVE) == surface.eval_d1(0.3, 0.7)
    box = surface.bounding_box()
    assert box is not None and box.max() == Point(2.0, 2.0, 4.0)
    assert repr(surface) == (
        "BezierSurface(u_degree=2, v_degree=2, pole_count=(3, 3), rational=False)"
    )
    normal = surface.normal(0.0, 0.0)
    assert normal is not None and (normal - Vector.unit_z()).magnitude() < 1e-12
    with pytest.raises(ValueError):
        surface.eval_derivative(0.5, 0.5, 0, 0)


def test_surface_rational_stays_on_the_cylinder():
    surface = quarter_cylinder()
    assert surface.is_rational()
    assert surface.weights() == [[1.0, 1.0], [W, W], [1.0, 1.0]]
    for i in range(11):
        p = surface.eval_point(i / 10, 0.5)
        assert math.hypot(p.x, p.y) == pytest.approx(1.0)
        assert p.z == pytest.approx(0.5)


def test_surface_errors():
    grid = [[Point.origin()] * 2] * 2
    with pytest.raises(ConstructionError):
        BezierSurface(0, 1, grid)
    with pytest.raises(ConstructionError):
        BezierSurface(2, 1, grid)
    with pytest.raises(ConstructionError):
        BezierSurface.new_rational(1, 1, grid, [[1.0, 1.0], [1.0]])
    with pytest.raises(ConstructionError):
        BezierSurface.new_rational(1, 1, grid, [[1.0, 1.0], [1.0, 0.0]])


def test_surface_editing_keeps_the_shape():
    surface = saddle()
    low, high = surface.split_u_at(0.25)
    assert (low.u_bounds(), high.u_bounds()) == ((0.0, 0.25), (0.25, 1.0))
    assert close(low.eval_point(0.1, 0.3), surface.eval_point(0.1, 0.3))
    assert close(high.eval_point(0.6, 0.3), surface.eval_point(0.6, 0.3))
    low, high = quarter_cylinder().split_v_at(0.4)
    assert (low.v_bounds(), high.v_bounds()) == ((0.0, 0.4), (0.4, 1.0))
    up = surface.elevate_u_degree(1)
    assert (up.u_degree(), up.u_pole_count()) == (3, 4)
    assert close(up.eval_point(0.3, 0.6), surface.eval_point(0.3, 0.6))
    assert surface.elevate_v_degree(2).v_degree() == 4
    assert close(surface.reversed_u().eval_point(0.2, 0.7), surface.eval_point(0.8, 0.7))
    assert close(surface.reversed_v().eval_point(0.2, 0.7), surface.eval_point(0.2, 0.3))
    with pytest.raises(geomcore.BSplineSplitError):
        surface.split_u_at(0.0)


def test_surface_round_trip_with_bspline():
    for surface in (saddle(), quarter_cylinder()):
        spline = surface.to_bspline()
        assert isinstance(spline, BSplineSurface)
        assert spline.u_multiplicities() == [surface.u_degree() + 1] * 2
        assert BezierSurface.from_bspline(spline) == surface
    with pytest.raises(ConstructionError):
        BezierSurface.from_bspline(saddle().to_bspline().insert_u_knot(0.5, 1))


def test_surface_projection_and_closedness():
    surface = saddle()
    on = surface.eval_point(0.4, 0.6)
    assert surface.contains(on)
    projection = surface.project_point(on)
    assert projection.distance == 0.0
    assert projection.u == pytest.approx(0.4, abs=1e-7)
    assert projection.target is surface
    assert not surface.contains(Point(0.0, 0.0, 5.0), Tolerance.default())
    assert len(surface.project_points([on, on])) == 2
    assert surface.extrema(on)[0].distance == 0.0
    assert not surface.is_u_closed() and not surface.is_v_closed()
    assert surface.curvatures(0.5, 0.5) is not None



# ---- Bézier operands of intersect ----


def test_bezier_operands_intersect_as_their_bspline():
    arch = BezierCurve(2, [Point(0.0, 0.0, 0.0), Point(1.0, 2.0, 0.0), Point(2.0, 0.0, 0.0)])
    line = Line(Point(0.0, 0.5, 0.0), Vector(1.0, 0.0, 0.0))
    result = line.intersect(arch)
    assert len(result.hits) == 2
    assert result.second is arch
    expected = line.intersect(arch.to_bspline())
    assert [h.point for h in result.hits] == [h.point for h in expected.hits]
    half = TrimmedCurve(arch, 0.0, 0.5)
    assert isinstance(half.basis(), BezierCurve)
    assert len(line.intersect(half).hits) == 1
    plane = Plane(Point(0.0, 0.5, 0.0), Vector(0.0, 1.0, 0.0))
    assert len(plane.intersect(arch).hits) == 2
    vertical = Line(Point(1.0, 1.0, 0.0), Vector(0.0, 0.0, 1.0))
    on_saddle = vertical.intersect(saddle())
    assert len(on_saddle.hits) == 1
    assert isinstance(on_saddle.surface, BezierSurface)
    level = Plane(Point(0.0, 0.0, 1.0), Vector(0.0, 0.0, 1.0))
    assert level.intersect(saddle()).curves
    arch_2d = BezierCurve2D(2, [Point2D(0.0, 0.0), Point2D(1.0, 2.0), Point2D(2.0, 0.0)])
    line_2d = Line2D(Point2D(0.0, 0.5), Vector2D(1.0, 0.0))
    assert len(line_2d.intersect(arch_2d).hits) == 2
