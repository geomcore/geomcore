"""Hits, overlaps and boxes name the same things the same way (#200)."""

import math

import pytest

from geomcore import (
    BezierSurface,
    BoundingBox,
    Circle,
    CurveSurfaceOverlap,
    Line,
    Plane,
    Point,
    Sphere,
    Vector,
)


def test_a_curve_surface_overlap_is_typed_and_names_its_range():
    line = Line(Point.origin(), Vector.unit_x())
    (overlap,) = line.intersect(Plane(Point.origin(), Vector.unit_z())).overlaps
    assert isinstance(overlap, CurveSurfaceOverlap)
    assert overlap.curve_range == (-math.inf, math.inf)
    assert repr(overlap) == "CurveSurfaceOverlap(curve_range=(-inf, inf))"

    equator = Circle(Point.origin(), Vector.unit_z(), 1.0)
    (on_sphere,) = equator.intersect(Sphere(Point.origin(), 1.0)).overlaps
    assert on_sphere.curve_range == pytest.approx((0.0, math.tau))


def test_a_curve_surface_hit_names_its_surface_parameters_u_and_v():
    # The unit square in z = 0: a numeric search, which reports (u, v).
    corners = [[Point(float(i), float(j), 0.0) for j in range(2)] for i in range(2)]
    square = BezierSurface(1, 1, corners)
    (hit,) = Line(Point(0.25, 0.5, -1.0), Vector.unit_z()).intersect(square).hits
    assert (hit.u, hit.v) == pytest.approx((0.25, 0.5))
    assert repr(hit) == (
        "CurveSurfaceHit(curve_parameter=1.0, u=0.25, v=0.5, "
        "point=Point(0.25, 0.5, 0.0), multiplicity=1)"
    )

    # A closed-form hit carries its surface parameters too (#177).
    plane = Plane(Point(0.0, 0.0, 2.0), Vector.unit_z())
    (on_plane,) = Line(Point.origin(), Vector.unit_z()).intersect(plane).hits
    assert (on_plane.u, on_plane.v) == (0.0, 0.0)


def test_a_bounding_box_names_its_corners_like_a_region():
    box = BoundingBox(Point(1.0, 2.0, 3.0), Point.origin())
    assert box.min_corner() == Point.origin()
    assert box.max_corner() == Point(1.0, 2.0, 3.0)
    assert not hasattr(box, "min") and not hasattr(box, "max")
