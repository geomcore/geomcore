"""Conics name their foci, directrices and eccentricity (#189)."""

import math

import pytest

import geomcore
from geomcore import (
    Cone,
    Ellipse,
    Ellipse2D,
    Hyperbola,
    Hyperbola2D,
    Line,
    Line2D,
    Parabola,
    Parabola2D,
    Point,
    Point2D,
    Vector,
    Vector2D,
)


def test_ellipse_points_sum_their_focal_distances_to_the_major_axis():
    ellipse = Ellipse(Point(1.0, 2.0, 0.0), Vector.unit_z(), Vector.unit_x(), 3.0, 1.2)
    near, far = ellipse.foci()
    assert ellipse.eccentricity() == pytest.approx(math.sqrt(1.0 - (1.2 / 3.0) ** 2))
    directrices = ellipse.directrices()
    assert directrices is not None
    for u in (0.0, 0.7, 2.5):
        point = ellipse.point_at(u)
        assert point.distance(near) + point.distance(far) == pytest.approx(6.0)
        to_directrix = directrices[0].project_point(point).distance
        assert point.distance(near) == pytest.approx(ellipse.eccentricity() * to_directrix)


def test_a_round_ellipse_has_no_directrices():
    round_ = Ellipse2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0), 2.0, 2.0)
    assert round_.eccentricity() == 0.0
    assert round_.directrices() is None
    assert round_.foci() == (Point2D(0.0, 0.0), Point2D(0.0, 0.0))


def test_hyperbola_other_branch_mirrors_through_the_center():
    hyperbola = Hyperbola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 2.0, 1.5)
    other = hyperbola.other_branch()
    assert other.point_at(0.0).x == pytest.approx(-2.0)
    near, far = hyperbola.foci()
    point = hyperbola.point_at(0.9)
    assert point.distance(far) - point.distance(near) == pytest.approx(4.0)
    assert hyperbola.eccentricity() == pytest.approx(2.5 / 2.0)
    assert len(hyperbola.directrices()) == 2
    planar = Hyperbola2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0), 2.0, 1.5)
    assert planar.other_branch().point_at(0.0).x == pytest.approx(-2.0)
    assert planar.foci()[0].x == pytest.approx(2.5)


def test_parabola_from_focus_and_directrix_keeps_equal_distances():
    directrix = Line(Point(0.0, -1.0, 0.0), Vector.unit_x())
    parabola = Parabola.from_focus_and_directrix(Point(0.0, 1.0, 0.0), directrix)
    assert parabola.eccentricity() == 1.0
    assert parabola.focal_distance() == pytest.approx(1.0)
    assert parabola.focus().distance(Point(0.0, 1.0, 0.0)) == pytest.approx(0.0)
    point = parabola.point_at(1.3)
    assert point.distance(parabola.focus()) == pytest.approx(
        parabola.directrix().project_point(point).distance
    )
    planar = Parabola2D.from_focus_and_directrix(
        Point2D(0.0, 1.0), Line2D(Point2D(0.0, -1.0), Vector2D(1.0, 0.0))
    )
    assert planar.apex() == Point2D(0.0, 0.0)
    assert planar.directrix().origin().y == pytest.approx(-1.0)
    with pytest.raises(geomcore.ConstructionError):
        Parabola.from_focus_and_directrix(Point(3.0, -1.0, 0.0), directrix)


def test_cone_from_apex_starts_at_its_apex():
    cone = Cone.from_apex(Point(1.0, 2.0, 3.0), Vector.unit_z(), 0.3)
    assert cone.apex().distance(Point(1.0, 2.0, 3.0)) == pytest.approx(0.0)
    assert cone.reference_radius() == 0.0


def test_the_package_reports_its_version():
    assert geomcore.__version__.count(".") == 2
