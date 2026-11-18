"""Curve–surface and surface–surface extrema, and the minimum distance of a
pair (#251): the bindings' shape, roles and errors. The geometry is
checked in Rust."""

import pytest

from geomcore import (
    BSplineSurface,
    Circle,
    CompletenessStatus,
    CurveSurfaceExtrema,
    CurveSurfaceExtremum,
    Line,
    Line2D,
    Plane,
    Point,
    Point2D,
    Sphere,
    SurfaceProjection,
    SurfaceSurfaceExtrema,
    SurfaceSurfaceExtremum,
    TrimmedCurve,
    Vector,
    Vector2D,
)


def test_a_curve_and_a_surface_give_their_minima_with_both_points():
    line = Line(Point(0.0, 0.0, 3.0), Vector(1.0, 0.0, 0.0))
    sphere = Sphere(Point.origin(), 1.0)
    extrema = line.extrema(sphere)
    assert isinstance(extrema, CurveSurfaceExtrema)
    assert extrema.is_complete()
    assert extrema.completeness.status == CompletenessStatus.COMPLETE
    (minimum,) = extrema.minima
    assert isinstance(minimum, CurveSurfaceExtremum)
    assert minimum.distance == pytest.approx(2.0)
    assert minimum.curve_point == Point(0.0, 0.0, 3.0)
    assert minimum.surface_point == sphere.point_at(*minimum.surface_parameters)
    assert line.point_at(minimum.curve_parameter) == minimum.curve_point
    assert extrema.first is line and extrema.second is sphere
    assert extrema.curve is line and extrema.surface is sphere


def test_a_surface_called_with_a_curve_keeps_the_roles():
    line = Line(Point(0.0, 0.0, 3.0), Vector(1.0, 0.0, 0.0))
    sphere = Sphere(Point.origin(), 1.0)
    extrema = sphere.extrema(line)
    assert extrema.first is sphere and extrema.second is line
    assert extrema.curve is line and extrema.surface is sphere
    assert extrema.minima == line.extrema(sphere).minima


def test_two_surfaces_give_their_minima():
    first = Sphere(Point.origin(), 1.0)
    second = Sphere(Point(5.0, 0.0, 0.0), 2.0)
    extrema = first.extrema(second)
    assert isinstance(extrema, SurfaceSurfaceExtrema)
    assert extrema.is_complete()
    (minimum,) = extrema
    assert isinstance(minimum, SurfaceSurfaceExtremum)
    assert minimum.distance == pytest.approx(2.0)
    assert minimum.first_point == Point(1.0, 0.0, 0.0)
    assert minimum.second_point == Point(3.0, 0.0, 0.0)
    assert extrema.first is first and extrema.second is second


def test_parallel_planes_report_their_distance_once():
    low = Plane(Point.origin(), Vector.unit_z())
    high = Plane(Point(0.0, 0.0, 2.0), Vector.unit_z())
    (minimum,) = low.extrema(high)
    assert minimum.distance == pytest.approx(2.0)


def test_a_curve_that_meets_a_surface_reports_each_crossing():
    line = Line(Point(-3.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0))
    sphere = Sphere(Point.origin(), 1.0)
    extrema = line.extrema(sphere)
    assert [m.distance for m in extrema] == pytest.approx([0.0, 0.0], abs=1e-12)


def test_the_minimum_distance_is_the_closest_minimum():
    circle = Circle(Point(0.0, 0.0, 4.0), Vector.unit_z(), 1.0)
    sphere = Sphere(Point.origin(), 1.0)
    plane = Plane(Point(0.0, 0.0, -2.0), Vector.unit_z())
    assert circle.minimum_distance(sphere) == circle.extrema(sphere).closest()
    assert sphere.minimum_distance(circle) == sphere.extrema(circle).closest()
    assert sphere.minimum_distance(plane) == sphere.extrema(plane).closest()
    assert sphere.minimum_distance(plane).distance == pytest.approx(1.0)
    segment = TrimmedCurve.segment(Point(-1.0, 3.0, 0.0), Point(1.0, 3.0, 0.0))
    assert circle.minimum_distance(segment) == circle.extrema(segment).closest()


def test_trimmed_surfaces_are_operands():
    plane = Plane(Point.origin(), Vector.unit_z())
    patch = plane.trimmed((-1.0, 1.0), (-1.0, 1.0))
    line = Line(Point(3.0, 0.0, 1.0), Vector(0.0, 1.0, 0.0))
    extrema = line.extrema(patch)
    assert isinstance(extrema, CurveSurfaceExtrema)
    assert extrema.surface is patch
    assert extrema.closest().distance == pytest.approx(5.0**0.5)
    assert patch.extrema(line).minima == extrema.minima
    assert isinstance(patch.extrema(Sphere(Point(0.0, 0.0, 5.0), 1.0)), SurfaceSurfaceExtrema)


def test_a_point_needs_a_bspline_or_bezier_surface():
    plane = Plane(Point.origin(), Vector.unit_z())
    with pytest.raises(TypeError, match="project_point"):
        plane.extrema(Point(0.0, 0.0, 1.0))  # type: ignore[call-overload]
    surface = BSplineSurface(
        1,
        1,
        [
            [Point(0.0, 0.0, 0.0), Point(0.0, 2.0, 0.0)],
            [Point(2.0, 0.0, 0.0), Point(2.0, 2.0, 0.0)],
        ],
        [0.0, 1.0],
        [2, 2],
        [0.0, 1.0],
        [2, 2],
        False,
        False,
    )
    (projection,) = surface.extrema(Point(1.0, 1.0, 1.0))
    assert isinstance(projection, SurfaceProjection)


def test_anything_else_raises_a_type_error():
    sphere = Sphere(Point.origin(), 1.0)
    line_2d = Line2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0))
    with pytest.raises(TypeError, match="not a 2D curve"):
        sphere.extrema(line_2d)  # type: ignore[call-overload]
    with pytest.raises(TypeError, match="a Point, a 3D curve or a surface, not str"):
        sphere.minimum_distance("sphere")  # type: ignore[call-overload]
    with pytest.raises(TypeError, match="a Point, a 3D curve or a surface, not int"):
        Line(Point.origin(), Vector.unit_x()).extrema(3)  # type: ignore[call-overload]
