"""Methods say what they compute in plain words (#202)."""

import math

import pytest

from geomcore import Axis, Point, Point2D, Transform, Vector, Vector2D
from geomcore.curves import Circle, Ellipse, Ellipse2D, Parabola, Parabola2D


def test_a_quarter_turn_of_a_planar_vector_is_its_perpendicular():
    assert Vector2D(3.0, 4.0).perpendicular() == Vector2D(-4.0, 3.0)
    assert not hasattr(Vector2D(3.0, 4.0), "perp")


def test_a_parabola_names_its_focal_distance():
    parabola = Parabola(Point.origin(), Vector.unit_z(), Vector.unit_x(), focal_distance=1.5)
    assert parabola.focal_distance() == 1.5
    planar = Parabola2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0), focal_distance=2.0)
    assert planar.focal_distance() == 2.0
    assert repr(planar).endswith("focal_distance=2.0)")


def test_composition_takes_the_other_transform():
    shift = Transform.translation(Vector.unit_x())
    turn = Transform.rotation(Axis(Point.origin(), Vector.unit_z()), math.pi / 2.0)
    point = shift.then(other=turn).apply_point(Point.origin())
    assert (point.x, point.y) == pytest.approx((0.0, 1.0))


def test_closed_conics_name_their_length_and_enclosed_area():
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    assert circle.length() == pytest.approx(math.tau * 2.0)
    assert circle.enclosed_area() == pytest.approx(math.pi * 4.0)
    ellipse = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    assert ellipse.enclosed_area() == pytest.approx(math.pi * 4.5)
    planar = Ellipse2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0), 3.0, 1.5)
    assert planar.enclosed_area() == pytest.approx(math.pi * 4.5)
    for old in ("circumference", "disk_area"):
        assert not hasattr(circle, old)
    assert not hasattr(ellipse, "area")
