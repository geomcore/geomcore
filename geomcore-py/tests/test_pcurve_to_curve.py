"""`pcurve_to_curve` carries a pcurve, a curve of a surface's (u, v)
parameter space, onto the surface, as the exact 3D curve at the same
parameter."""

import math

import pytest

import geomcore
from geomcore import Point, Point2D, Vector, Vector2D
from geomcore.curves import Circle, Circle2D, Ellipse, Ellipse2D, Line, Line2D
from geomcore.surfaces import Cylinder, Plane, Sphere


def assert_traces(traced, curve, surface):
    for t in (-1.2, 0.0, 0.4, 2.5):
        q = curve.eval_point(t)
        assert traced.eval_point(t).distance(surface.eval_point(q.x, q.y)) < 1e-12


def test_lines_along_u_and_v_on_a_cylinder():
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 2.0)
    along_v = Line2D(Point2D(0.5, 1.0), Vector2D.unit_y())
    along_u = Line2D(Point2D(0.5, 1.0), Vector2D(-1.0, 0.0))
    generator = cylinder.pcurve_to_curve(along_v)
    parallel = cylinder.pcurve_to_curve(along_u)
    assert isinstance(generator, Line)
    assert isinstance(parallel, Circle)
    assert_traces(generator, along_v, cylinder)
    assert_traces(parallel, along_u, cylinder)


def test_curves_on_a_plane_keep_their_kind():
    plane = Plane(Point(0, 0, 1), Vector(0, 0, 1))
    circle = Circle2D(Point2D(1.0, 2.0), 0.5)
    ellipse = Ellipse2D(Point2D(1.0, 2.0), Vector2D(0.6, 0.8), 3.0, 1.0)
    assert isinstance(plane.pcurve_to_curve(circle), Circle)
    traced = plane.pcurve_to_curve(ellipse, geomcore.Tolerance.default())
    assert isinstance(traced, Ellipse)
    assert_traces(traced, ellipse, plane)


def test_refused_pairs_raise():
    sphere = Sphere(Point(0, 0, 0), 2.0)
    circle = Circle2D(Point2D(1.0, 0.0), 0.5)
    with pytest.raises(geomcore.PcurveToCurveError):
        sphere.pcurve_to_curve(circle)
    pole = Line2D(Point2D(0.0, math.pi / 2), Vector2D.unit_x())
    with pytest.raises(geomcore.PcurveToCurveError, match="point"):
        sphere.pcurve_to_curve(pole)


def test_a_3d_curve_is_rejected():
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 2.0)
    line = Line(Point(0, 0, 0), Vector(0, 0, 1))
    with pytest.raises(TypeError):
        cylinder.pcurve_to_curve(line)
