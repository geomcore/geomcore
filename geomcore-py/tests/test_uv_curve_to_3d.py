"""`uv_curve_to_3d` carries a UV curve, a curve of a surface's (u, v)
parameter space, onto the surface, as the exact 3D curve at the same
parameter."""

import math

import pytest

import geomcore
from geomcore import (
    BSplineCurve,
    BSplineCurve2D,
    BezierCurve2D,
    BezierSurface,
    Circle,
    Circle2D,
    Cylinder,
    Ellipse,
    Ellipse2D,
    Line,
    Line2D,
    Plane,
    Point,
    Point2D,
    Sphere,
    TrimmedCurve,
    TrimmedCurve2D,
    Vector,
    Vector2D,
)


def assert_traces(traced, curve, surface, parameters=(-1.2, 0.0, 0.4, 2.5)):
    for t in parameters:
        q = curve.point_at(t)
        assert traced.point_at(t).distance(surface.point_at(q.x, q.y)) < 1e-12


def test_lines_along_u_and_v_on_a_cylinder():
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 2.0)
    along_v = Line2D(Point2D(0.5, 1.0), Vector2D.unit_y())
    along_u = Line2D(Point2D(0.5, 1.0), Vector2D(-1.0, 0.0))
    generator = cylinder.uv_curve_to_3d(along_v)
    parallel = cylinder.uv_curve_to_3d(along_u)
    assert isinstance(generator, Line)
    assert isinstance(parallel, Circle)
    assert_traces(generator, along_v, cylinder)
    assert_traces(parallel, along_u, cylinder)


def test_curves_on_a_plane_keep_their_kind():
    plane = Plane(Point(0, 0, 1), Vector(0, 0, 1))
    circle = Circle2D(Point2D(1.0, 2.0), 0.5)
    ellipse = Ellipse2D(Point2D(1.0, 2.0), Vector2D(0.6, 0.8), 3.0, 1.0)
    assert isinstance(plane.uv_curve_to_3d(circle), Circle)
    traced = plane.uv_curve_to_3d(ellipse, geomcore.Tolerance.default())
    assert isinstance(traced, Ellipse)
    assert_traces(traced, ellipse, plane)


def test_refused_pairs_raise():
    sphere = Sphere(Point(0, 0, 0), 2.0)
    circle = Circle2D(Point2D(1.0, 0.0), 0.5)
    with pytest.raises(geomcore.UvCurveTo3DError):
        sphere.uv_curve_to_3d(circle)
    pole = Line2D(Point2D(0.0, math.pi / 2), Vector2D.unit_x())
    with pytest.raises(geomcore.UvCurveTo3DError, match="point"):
        sphere.uv_curve_to_3d(pole)


def test_a_3d_curve_is_rejected():
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 2.0)
    line = Line(Point(0, 0, 0), Vector(0, 0, 1))
    with pytest.raises(TypeError):
        cylinder.uv_curve_to_3d(line)


def test_a_straight_spline_runs_along_its_iso_curve():
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 2.0)
    poles = [Point2D(0.5, 2.0), Point2D(0.5, 0.0)]
    along_v = BSplineCurve2D(1, poles, [1.0, 3.0], [2, 2], False)
    generator = cylinder.uv_curve_to_3d(along_v)
    assert isinstance(generator, Line)
    assert_traces(generator, along_v, cylinder, (1.2, 2.0, 2.9))


def bezier_surface(degree):
    poles = [
        [
            Point(i / degree, j / degree, math.sin(i + 2 * j) / degree)
            for j in range(degree + 1)
        ]
        for i in range(degree + 1)
    ]
    return BezierSurface(degree, degree, poles)


def test_a_slanted_line_on_a_spline_surface_is_composed():
    surface = bezier_surface(3)
    line = Line2D(Point2D(0.5, 0.5), Vector2D(0.6, 0.8))
    traced = surface.uv_curve_to_3d(line)
    assert isinstance(traced, BSplineCurve)
    assert traced.degree() == 6
    first, last = traced.bounds()
    assert_traces(traced, line, surface, (first, -0.1, 0.2, last))


def test_a_composition_outside_the_domain_or_above_degree_25_raises():
    surface = bezier_surface(3)
    outside = Line2D(Point2D(5.0, 3.0), Vector2D(0.6, -0.8))
    with pytest.raises(geomcore.UvCurveTo3DError, match="domain"):
        surface.uv_curve_to_3d(outside)
    poles = [Point2D(0.1, 0.1), Point2D(0.3, 0.9), Point2D(0.7, 0.1), Point2D(0.9, 0.9)]
    cubic = BezierCurve2D(3, poles)
    with pytest.raises(geomcore.UvCurveTo3DError, match="degree 30"):
        bezier_surface(5).uv_curve_to_3d(cubic)


def test_a_piece_traces_the_piece_of_its_basis_trace():
    # #176: the stubs took any 2D curve, the runtime refused a piece.
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 2.0)
    piece = TrimmedCurve2D(Line2D(Point2D(0.5, 1.0), Vector2D(-1.0, 0.0)), 0.2, 1.5)
    arc = cylinder.uv_curve_to_3d(piece)
    assert isinstance(arc, TrimmedCurve)
    assert isinstance(arc.basis(), Circle)
    assert arc.bounds() == piece.bounds()
    assert_traces(arc, piece, cylinder, (0.2, 0.9, 1.5))
    segment = TrimmedCurve2D.segment(Point2D(0.1, 0.1), Point2D(0.8, 0.6))
    on_patch = bezier_surface(3).uv_curve_to_3d(segment)
    assert on_patch.bounds() == segment.bounds()
    assert_traces(on_patch, segment, bezier_surface(3), segment.bounds())


def test_a_spline_piece_needs_only_the_piece_in_the_domain():
    surface = bezier_surface(2)
    poles = [Point2D(-1.0, 0.5), Point2D(0.5, 0.5), Point2D(2.0, 0.5)]
    leaving = BSplineCurve2D(1, poles, [0.0, 1.0, 2.0], [2, 1, 2], False)
    with pytest.raises(geomcore.UvCurveTo3DError, match="domain"):
        surface.uv_curve_to_3d(leaving)
    inside = TrimmedCurve2D(leaving, 0.8, 1.2)
    traced = surface.uv_curve_to_3d(inside)
    assert traced.bounds() == (0.8, 1.2)
    assert_traces(traced, inside, surface, (0.8, 1.0, 1.2))


def test_a_piece_off_the_domain_or_around_twice_raises():
    line = Line2D(Point2D(0.5, 0.5), Vector2D(0.6, 0.8))
    with pytest.raises(geomcore.UvCurveTo3DError, match="domain"):
        bezier_surface(3).uv_curve_to_3d(TrimmedCurve2D(line, 0.0, 3.0))
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 2.0)
    around_twice = TrimmedCurve2D(Line2D(Point2D(0.0, 0.0), Vector2D.unit_x()), 0.0, 13.0)
    with pytest.raises(geomcore.TrimmedCurveConstructionError, match="period"):
        cylinder.uv_curve_to_3d(around_twice)


def test_a_slanted_line_on_a_cylinder_maps_within_a_deviation():
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 10.0)
    slope = Vector2D(2.0 * math.pi, 5.0)
    line = Line2D(Point2D(0.0, 0.0), slope)
    with pytest.raises(geomcore.UvCurveTo3DError):
        cylinder.uv_curve_to_3d(line)
    turn = TrimmedCurve2D(line, 0.0, slope.magnitude())
    approximation = cylinder.uv_curve_to_3d_within(turn, 1e-7)
    assert approximation.is_proven()
    assert approximation.deviation <= 1e-7
    curve = approximation.curve
    assert isinstance(curve, BSplineCurve)
    assert curve.degree() == 9
    assert curve.bounds() == (0.0, slope.magnitude())
    for t in (0.0, 1.0, 3.3, slope.magnitude()):
        uv = turn.point_at(t)
        assert curve.point_at(t).distance(cylinder.point_at(uv.x, uv.y)) <= 1e-7
    assert "degree=9" in repr(approximation)


def test_a_circle_on_a_sphere_maps_within_a_sampled_deviation():
    sphere = Sphere(Point(0, 0, 0), 2.0)
    circle = Circle2D(Point2D(1.0, 0.0), 0.5)
    approximation = sphere.uv_curve_to_3d_within(circle, 1e-6, geomcore.Tolerance.default())
    assert not approximation.is_proven()
    assert approximation.deviation <= 1e-6
    first, last = approximation.curve.bounds()
    assert (first, last) == pytest.approx(circle.bounds())


def test_uv_curve_to_3d_within_refuses_with_its_own_error():
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 10.0)
    line = Line2D(Point2D(0.0, 0.0), Vector2D(1.0, 1.0))
    with pytest.raises(geomcore.UvCurveTo3DApproximationError, match="unbounded"):
        cylinder.uv_curve_to_3d_within(line, 1e-7)
    turn = TrimmedCurve2D(line, 0.0, 1.0)
    with pytest.raises(geomcore.UvCurveTo3DApproximationError, match="deviation"):
        cylinder.uv_curve_to_3d_within(turn, -1.0)



def test_the_free_functions_warn_and_give_the_methods_results():
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 10.0)
    slope = Vector2D(2.0 * math.pi, 5.0)
    turn = TrimmedCurve2D(Line2D(Point2D(0.0, 0.0), slope), 0.0, slope.magnitude())
    with pytest.warns(DeprecationWarning, match="surface.uv_curve_to_3d_within"):
        traced = geomcore.uv_curve_to_3d_within(cylinder, turn, 1e-7)  # type: ignore[deprecated]
    assert traced.curve == cylinder.uv_curve_to_3d_within(turn, 1e-7).curve
    with pytest.warns(DeprecationWarning, match="surface.curve_3d_to_uv_within"):
        mapped = geomcore.curve_3d_to_uv_within(cylinder, traced.curve, 1e-6)  # type: ignore[deprecated]
    assert mapped.curve == cylinder.curve_3d_to_uv_within(traced.curve, 1e-6).curve


def test_the_helix_recipe_of_the_parametrization_guide():
    # #257: the recipe in docs/guides/parametrization.md, and its comments.
    from geomcore import Cylinder, Line2D, Point2D, TrimmedCurve2D, Vector2D

    radius, pitch, turns = 5.0, 2.0, 3.0
    cylinder = Cylinder(geomcore.Point.origin(), geomcore.Vector.unit_z(), radius)
    uv_line = Line2D(Point2D(0.0, 0.0), Vector2D(2 * math.pi, pitch))
    uv_piece = TrimmedCurve2D(uv_line, 0.0, turns * math.hypot(2 * math.pi, pitch))
    approximation = cylinder.uv_curve_to_3d_within(uv_piece, 1e-6)
    helix = approximation.curve
    assert (helix.degree(), len(helix.poles())) == (9, 60)
    assert approximation.is_proven()
    assert approximation.deviation < 4.2e-7
    end = helix.end_point()
    assert (end.x, end.y, end.z) == pytest.approx((5.0, 0.0, 6.0), abs=1e-12)
    for k in range(1, 30):
        point = helix.point_at(helix.bounds()[1] * k / 30)
        assert math.hypot(point.x, point.y) == pytest.approx(radius, abs=1e-6)
        assert point.z == pytest.approx(pitch * turns * k / 30, abs=1e-6)
    # Right-handed: it starts turning counterclockwise seen from above.
    assert helix.point_at(0.1).y > 0.0
