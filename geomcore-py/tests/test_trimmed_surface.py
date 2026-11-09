"""TrimmedSurface: a patch of a surface, and the intersections clipped to it."""

import math
import pickle

import pytest

import geomcore
from geomcore import (
    Cylinder,
    Ellipse,
    Plane,
    Point,
    Sphere,
    SurfaceSurfaceIntersection,
    TrimmedSurface,
    TrimmedSurfaceConstructionError,
    Vector,
)


def cylinder():
    return Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 1.0)


def test_every_surface_has_trimmed():
    for cls in geomcore.surfaces.__all__:
        assert hasattr(getattr(geomcore, cls), "trimmed"), cls


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
