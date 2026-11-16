"""TrimmedSurface: a patch of a surface, and the intersections clipped to it."""

import math
import pickle
import typing

import pytest

import geomcore
from geomcore import (
    BSplineSurface,
    Cylinder,
    Ellipse,
    Line2D,
    Plane,
    Point,
    Point2D,
    Sphere,
    SurfaceSurfaceIntersection,
    TrimmedSurface,
    TrimmedSurfaceConstructionError,
    Vector,
    Vector2D,
)


def cylinder():
    return Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 1.0)


def test_every_surface_has_trimmed():
    for cls in typing.get_args(geomcore.Surface):
        assert hasattr(cls, "trimmed"), cls.__name__


def test_trimmed_keeps_the_basis_and_ranges():
    patch = cylinder().trimmed((0.0, math.pi), (0.0, 2.0))
    assert isinstance(patch, TrimmedSurface)
    assert patch.basis() == cylinder()
    assert patch.u_range() == (0.0, math.pi)
    assert patch.v_bounds() == (0.0, 2.0)
    assert patch == TrimmedSurface(cylinder(), (0.0, math.pi), (0.0, 2.0))
    assert pickle.loads(pickle.dumps(patch)) == patch
    assert repr(patch).startswith("TrimmedSurface(basis=Cylinder(")


def test_invalid_ranges_raise():
    with pytest.raises(TrimmedSurfaceConstructionError):
        cylinder().trimmed((0.0, 7.0), (0.0, 1.0))
    with pytest.raises(TrimmedSurfaceConstructionError):
        cylinder().trimmed((1.0, 1.0), (0.0, 1.0))
    with pytest.raises(TypeError):
        TrimmedSurface(Point(0, 0, 0), (0.0, 1.0), (0.0, 1.0))


def test_a_plane_across_a_half_cylinder_gives_half_an_ellipse():
    half = cylinder().trimmed((0.0, math.pi), (-5.0, 5.0))
    plane = Plane(Point(0, 0, 0), Vector(0, 1, 1))
    for result in (plane.intersect(half), half.intersect(plane)):
        assert isinstance(result, SurfaceSurfaceIntersection)
        assert len(result.curves) == 1
        piece = result.curves[0]
        assert isinstance(piece.curve, Ellipse)
        start, end = piece.parameter_range
        assert end - start == pytest.approx(math.pi)


def test_two_patches_and_untrimmed_parameter_ranges():
    plane = Plane(Point(0, 0, 0.5), Vector(0, 0, 1))
    result = plane.intersect(cylinder())
    assert result.curves[0].parameter_range == (0.0, 2 * math.pi)
    a = plane.trimmed((-1.0, 0.0), (-2.0, 2.0))
    b = cylinder().trimmed((0.0, 2 * math.pi), (0.0, 1.0))
    clipped = a.intersect(b)
    # The circle's half with x <= 0.
    assert len(clipped.curves) == 1
    start, end = clipped.curves[0].parameter_range
    assert end - start == pytest.approx(math.pi)


def test_contact_points_outside_the_patch_are_dropped():
    sphere = Sphere(Point(0, 0, 1), 1.0)
    aside = Plane(Point(0, 0, 0), Vector(0, 0, 1)).trimmed((1.0, 2.0), (1.0, 2.0))
    assert aside.intersect(sphere).is_empty()


def test_a_curve_and_a_trimmed_surface_raise():
    patch = cylinder().trimmed((0.0, math.pi), (0.0, 2.0))
    line = geomcore.Line(Point(0, 0, 0), Vector(1, 0, 0))
    with pytest.raises(TypeError):
        patch.intersect(line)
    with pytest.raises(TypeError):
        line.intersect(patch)


def test_a_patch_answers_what_its_basis_answers():
    # #290: the patch took only 14 hand-written methods of its basis's.
    basis = cylinder()
    patch = basis.trimmed((0.0, 3.0), (-1.0, 1.0))
    assert patch.second_derivatives(1.0, 0.5) == basis.second_derivatives(1.0, 0.5)
    assert patch.derivative(1.0, 0.5, 3, 0) == basis.derivative(1.0, 0.5, 3, 0)
    assert patch.curvatures(1.0, 0.5).mean == basis.curvatures(1.0, 0.5).mean
    assert patch.continuity_u(1.0) == basis.continuity_u(1.0)
    assert patch.u_period() is None and not patch.is_u_periodic()
    points = [Point(2, 0, 0.5), Point(0, 2, 0), Point(-2, -0.1, 0)]
    projections = patch.project_points(points)
    assert [(p.u, p.v) for p in projections] == [
        (patch.project_point(point).u, patch.project_point(point).v) for point in points
    ]
    assert patch.points_at_array([[1.0, 0.5]]).tolist() == [list(patch.point_at(1.0, 0.5))]
    along_v = Line2D(Point2D(1.0, 0.0), Vector2D(0.0, 1.0))
    assert patch.uv_curve_to_3d(along_v).point_at(0.5) == basis.uv_curve_to_3d(along_v).point_at(0.5)
    assert patch.tight_bounding_box() is not None


def test_a_patch_checks_a_derivative_order_against_its_basis():
    grid = [[Point(i, j, (i - 1) * (j - 1)) for j in range(3)] for i in range(3)]
    spline = BSplineSurface.approximate(grid, 2, 2, 3, 3)
    patch = spline.trimmed((0.2, 0.8), (0.2, 0.8))
    patch.derivative(0.5, 0.5, 1, 1)
    with pytest.raises(ValueError, match=r"du \+ dv must be <= 2"):
        patch.derivative(0.5, 0.5, 2, 1)
