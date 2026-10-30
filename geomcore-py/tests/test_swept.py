"""The swept surfaces: a profile curve extruded along a direction or
turned about an axis."""

import math

import pytest

from geomcore import Axis, Frame, Point, Transform, Vector
from geomcore.curves import BSplineCurve, Circle, Ellipse, Line
from geomcore.surfaces import (
    BSplineSurface,
    Cylinder,
    ExtrusionSurface,
    RevolutionSurface,
    Torus,
)

_O, _X, _Y, _Z = Point.origin(), Vector.unit_x(), Vector.unit_y(), Vector.unit_z()


def _close(p, q, tol=1e-9):
    return p.distance(q) <= tol


def _torus_profile():
    """The circle meridian of `Torus(origin, z, 3, 1)` at u = 0."""
    return Circle.from_frame(Frame(Point(3.0, 0.0, 0.0), -_Y, _X), 1.0)


def test_extruded_circle_is_a_cylinder():
    surface = ExtrusionSurface(Circle(_O, _Z, 2.0), Vector(0.0, 0.0, 5.0))
    cylinder = Cylinder(_O, _Z, 2.0)
    assert surface.direction() == _Z
    assert isinstance(surface.profile(), Circle)
    for u, v in [(0.0, 0.0), (1.0, -2.0), (4.0, 3.5)]:
        assert _close(surface.eval_point(u, v), cylinder.eval_point(u, v))
    assert surface.u_period() == pytest.approx(2.0 * math.pi)
    assert surface.v_bounds() == (-math.inf, math.inf)
    assert surface.bounding_box() is None
    projection = surface.project_point(Point(4.0, 0.0, 1.0))
    assert projection.distance == pytest.approx(2.0)
    assert (projection.u, projection.v) == pytest.approx((0.0, 1.0))


def test_extrusion_iso_curves_and_reversal():
    profile = Ellipse(_O, _Z, _X, 3.0, 1.0)
    surface = ExtrusionSurface(profile, Vector(0.0, 1.0, 1.0))
    line = surface.u_iso(0.5)
    assert isinstance(line, Line)
    assert _close(line.eval_point(2.0), surface.eval_point(0.5, 2.0))
    ellipse = surface.v_iso(2.0)
    assert isinstance(ellipse, Ellipse)
    assert _close(ellipse.eval_point(0.5), surface.eval_point(0.5, 2.0))
    reversed_u = surface.reversed_u()
    assert _close(reversed_u.eval_point(surface.u_reversed_parameter(0.5), 1.0), surface.eval_point(0.5, 1.0))
    reversed_v = surface.reversed_v()
    assert _close(reversed_v.eval_point(0.5, surface.v_reversed_parameter(1.0)), surface.eval_point(0.5, 1.0))


def test_extrusion_patch_lies_on_the_surface():
    surface = ExtrusionSurface(Ellipse(_O, _Z, _X, 3.0, 1.0), Vector(0.0, 1.0, 1.0))
    patch = surface.to_bspline_patch((0.0, 2.0), (-1.0, 1.0))
    assert isinstance(patch, BSplineSurface)
    for u, v in [(0.0, -1.0), (0.7, 0.3), (2.0, 1.0)]:
        on = patch.eval_point(u, v)
        assert surface.project_point(on).distance < 1e-9


def test_extrusion_construction_errors():
    with pytest.raises(ValueError):
        ExtrusionSurface(Circle(_O, _Z, 1.0), Vector(0.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        ExtrusionSurface(Line(_O, _Z), _Z)
    with pytest.raises(TypeError):
        ExtrusionSurface(Cylinder(_O, _Z, 1.0), _Z)


def test_revolved_circle_is_a_torus():
    surface = RevolutionSurface(_torus_profile(), Axis(_O, _Z))
    torus = Torus(_O, _Z, 3.0, 1.0)
    assert surface.axis() == Axis(_O, _Z)
    for u, v in [(0.0, 0.0), (1.0, 2.0), (4.0, 5.5)]:
        assert _close(surface.eval_point(u, v), torus.eval_point(u, v))
    box = surface.bounding_box()
    assert box is not None
    assert box.min().z == pytest.approx(-1.0) and box.max().x == pytest.approx(4.0)
    query = Point(0.0, 5.0, 0.5)
    assert surface.project_point(query).distance == pytest.approx(torus.project_point(query).distance)
    # A rational circle does not keep the angle as its parameter, so the
    # patch reaches the same points at other parameters.
    whole = surface.to_bspline()
    assert torus.project_point(whole.eval_point(1.0, 2.0)).distance < 1e-9


def test_revolution_iso_curves_and_transform():
    profile = BSplineCurve(
        2,
        [Point(1.0, 0.0, 0.0), Point(2.0, 0.0, 1.0), Point(1.5, 0.0, 2.0)],
        [0.0, 1.0],
        [3, 3],
        False,
    )
    surface = RevolutionSurface(profile, Axis(_O, _Z))
    meridian = surface.u_iso(1.0)
    assert isinstance(meridian, BSplineCurve)
    assert _close(meridian.eval_point(0.4), surface.eval_point(1.0, 0.4))
    parallel = surface.v_iso(0.4)
    assert isinstance(parallel, Circle)
    assert _close(parallel.eval_point(1.0), surface.eval_point(1.0, 0.4))
    moved = surface.transformed(Transform.translation(Vector(0.0, 0.0, 2.0)))
    assert _close(moved.eval_point(1.0, 0.4), surface.eval_point(1.0, 0.4) + Vector(0.0, 0.0, 2.0))
    assert surface.parameters_of(surface.eval_point(1.0, 0.4)) == pytest.approx((1.0, 0.4))


def test_revolution_construction_errors():
    with pytest.raises(ValueError):
        RevolutionSurface(Line(_O, _Z), Axis(_O, _Z))
    with pytest.raises(ValueError):
        RevolutionSurface(Circle(_O, _Z, 1.0), Axis(_O, _Z))


def test_swept_surfaces_intersect_numerically():
    surface = RevolutionSurface(_torus_profile(), Axis(_O, _Z))
    result = surface.intersect(Line(Point(-5.0, 0.0, 0.0), _X))
    assert sorted(h.point.x for h in result.hits) == pytest.approx([-4.0, -2.0, 2.0, 4.0], abs=1e-7)
    extrusion = ExtrusionSurface(Circle(_O, _Z, 2.0), _Z)
    hits = extrusion.intersect(Line(Point(-5.0, 0.0, 1.0), _X)).hits
    assert sorted(h.point.x for h in hits) == pytest.approx([-2.0, 2.0], abs=1e-7)
