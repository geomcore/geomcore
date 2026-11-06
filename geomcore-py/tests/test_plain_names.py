"""Methods say what they compute in plain words (#202)."""

import enum
import math

import pytest

from geomcore import (
    Axis,
    Circle,
    Ellipse,
    Ellipse2D,
    Line,
    Parabola,
    Parabola2D,
    Plane,
    Point,
    Point2D,
    Transform,
    Vector,
    Vector2D,
)


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


def test_enumerations_are_real_enums_that_list_their_members():
    from geomcore import (
        CompletenessStatus,
        IncompleteReason,
        IntersectionMethod,
        PointParametrization,
        Side,
    )

    for enumeration in (Side, IntersectionMethod, CompletenessStatus, IncompleteReason):
        assert issubclass(enumeration, enum.IntEnum)
    assert [member.name for member in IncompleteReason] == [
        "BUDGET_EXHAUSTED",
        "UNSEPARATED",
        "UNBOUNDED",
    ]
    assert list(PointParametrization) == ["chordal", "centripetal"]
    assert CompletenessStatus.__module__ == "geomcore"
    assert "numeric search proved" in (CompletenessStatus.__doc__ or "")
    line = Line(Point.origin(), Vector.unit_z())
    result = line.intersect(Plane(Point(0.0, 0.0, 2.0), Vector.unit_z()))
    assert result.method is IntersectionMethod.EXACT
    assert result.completeness.status is CompletenessStatus.COMPLETE


def test_a_parametrization_is_an_enum_member_or_its_string():
    from geomcore import BSplineCurve, PointParametrization

    points = [Point(0.0, 0.0, 0.0), Point(1.0, 2.0, 0.0), Point(3.0, 2.0, 0.0), Point(4.0, 0.0, 0.0)]
    by_member = BSplineCurve.interpolate(points, 3, PointParametrization.CHORDAL)
    by_string = BSplineCurve.interpolate(points, 3, "chordal")
    # Four points at degree 3 have no interior knot: the spacing shows in
    # the curve between the points instead.
    assert by_member.point_at(0.5) == by_string.point_at(0.5)
    assert by_member.point_at(0.5) != BSplineCurve.interpolate(points, 3).point_at(0.5)
    with pytest.raises(ValueError, match="chordale"):
        BSplineCurve.interpolate(points, 3, "chordale")
