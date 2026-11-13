"""Error messages name the cause and the next step (#335): one test per
message the UX study found naming the wrong cause or stopping short."""

import re

import numpy
import pytest

from geomcore import (
    BSplineCurve,
    BSplineSurface,
    Circle,
    Circle2D,
    Cone,
    Ellipse,
    Line,
    Plane,
    Point,
    Point2D,
    PointParametrization,
    TrimmedCurve,
    Vector,
)

ORIGIN = Point.origin()
Z = Vector.unit_z()
FOUR_POINTS = [Point(0, 0, 0), Point(1, 1, 0), Point(2, 0, 0), Point(3, 1, 0)]


def raises(exception, text):
    return pytest.raises(exception, match=re.escape(text))


def test_an_x_direction_along_the_normal_is_named():
    with raises(ValueError, "x direction is parallel to the normal"):
        Ellipse(ORIGIN, Z, Z, 3.0, 2.0)


def test_two_dimensional_extrema_say_what_they_take():
    circle = Circle2D(Point2D.origin(), 1.0)
    with raises(
        TypeError,
        "extrema of a 2D curve takes a Point2D, not Circle2D; "
        "the extrema between two 2D curves are not supported",
    ):
        circle.extrema(Circle2D(Point2D(3.0, 0.0), 1.0))


def test_a_cone_half_angle_gives_its_unit_and_range():
    message = "half-angle must be in radians, with 0 < |half-angle| < pi/2"
    with raises(ValueError, message):
        Cone(ORIGIN, Z, 45.0, 1.0)
    with raises(ValueError, message):
        Cone.from_apex(ORIGIN, Z, 30.0)


def test_an_empty_trimmed_range_uses_the_parameter_names():
    with raises(ValueError, "the start parameter is not less than the end parameter"):
        TrimmedCurve(Circle(ORIGIN, Z, 1.0), 2.0, 1.0)


def test_a_split_outside_the_bounds_gives_the_bounds():
    spline = BSplineCurve.interpolate(FOUR_POINTS, 3)
    with raises(ValueError, "a parameter is outside the bounds: the bounds are [0.0, 1.0]"):
        spline.split_at(5.0)
    with raises(ValueError, "the bounds are [0.0, 1.0]"):
        spline.extract(0.5, 5.0)


def test_interpolation_counts_the_points_a_degree_needs():
    with raises(ValueError, "degree 3 needs at least 4 points, got 2"):
        BSplineCurve.interpolate([Point(0, 0, 0), Point(1, 0, 0)], 3)


def test_an_unknown_parametrization_lists_the_valid_ones():
    with raises(
        ValueError,
        "'uniform' is not a valid PointParametrization; "
        "use a member or one of 'chordal', 'centripetal'",
    ):
        BSplineCurve.interpolate(FOUR_POINTS, 3, "uniform")
    assert PointParametrization("chordal") is PointParametrization.CHORDAL


def test_a_clamped_uniform_spline_counts_the_poles_a_degree_needs():
    with raises(ValueError, "degree 3 needs at least 4 poles, got 3"):
        BSplineCurve.clamped_uniform(3, FOUR_POINTS[:3])


def test_normalized_warns_that_it_may_give_none():
    assert Vector.zero().normalized() is None
    assert "check the result" in (Vector.normalized.__doc__ or "")
    with pytest.raises(TypeError):
        Line(ORIGIN, Vector.zero().normalized())


def test_a_curve_off_the_surface_names_its_largest_distance():
    # The segment climbs from the plane to 0.995 above it; the refusal
    # names that, not the first sample past half the deviation.
    segment = TrimmedCurve(Line(ORIGIN, Vector(1.0, 0.0, 0.1)), 0.0, 10.0)
    with pytest.raises(ValueError, match=r"the curve is 0\.995\d* from the surface"):
        Plane(ORIGIN, Z).curve_3d_to_uv_within(segment, 0.2)


def test_a_grid_shape_names_two_free_axes():
    with raises(ValueError, "points must have shape (N, M, 3), not (4, 4)"):
        BSplineSurface.approximate(numpy.zeros((4, 4)), 3, 3, 4, 4)
