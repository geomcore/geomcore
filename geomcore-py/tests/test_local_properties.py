"""Local properties, closedness, continuity, surface reversal and one-sided
evaluation: the queries added in 0.36."""

import math

import pytest

import geomcore
from geomcore import (
    Axis3D,
    Frame3D,
    FrameConstructionError,
    Point3D,
    Side,
    SurfaceCurvatures,
    Tolerance,
    Transform,
    Vector3D,
)
from geomcore.curves import (
    BSplineCurve3D,
    Circle3D,
    Ellipse3D,
    Hyperbola3D,
    Line3D,
    Parabola3D,
    TrimmedCurve3D,
)
from geomcore.surfaces import BSplineSurface, Cone, Cylinder, Plane, Sphere, Torus

ORIGIN = Point3D.origin()
X, Y, Z = Vector3D.unit_x(), Vector3D.unit_y(), Vector3D.unit_z()


def coordinates(value):
    return (value.x, value.y, value.z)


def close(a, b, abs=1e-12):
    return coordinates(a) == pytest.approx(coordinates(b), abs=abs)


def kinked_polyline():
    """Two straight pieces meeting at a right angle at u = 1."""
    poles = [Point3D(0.0, 0.0, 0.0), Point3D(1.0, 0.0, 0.0), Point3D(1.0, 1.0, 0.0)]
    return BSplineCurve3D(1, poles, [0.0, 1.0, 2.0], [2, 1, 2], False)


def periodic_quadratic():
    """A closed quadratic through four poles, C1 at each of its knots."""
    poles = [
        Point3D(1.0, 0.0, 0.0),
        Point3D(0.0, 1.0, 0.0),
        Point3D(-1.0, 0.0, 0.5),
        Point3D(0.0, -1.0, 0.0),
    ]
    return BSplineCurve3D(2, poles, [0.0, 1.0, 2.0, 3.0, 4.0], [1, 1, 1, 1, 1], True)


def creased_surface():
    """Two bilinear patches meeting at a crease along u = 1, flat then rising."""

    def flat(x):
        return [Point3D(x, 0.0, 0.0), Point3D(x, 1.0, 0.0)]

    def tilted(x):
        return [Point3D(x, 0.0, 0.0), Point3D(x, 1.0, 1.0)]

    return BSplineSurface(
        1, 1, [flat(0.0), flat(1.0), tilted(2.0)],
        [0.0, 1.0, 2.0], [2, 1, 2], [0.0, 1.0], [2, 2], False, False,
    )


# ---------------------------------------------------------------------------
# Curve tangent, curvature, centre of curvature
# ---------------------------------------------------------------------------


def test_circle_local_properties_match_closed_form():
    circle = Circle3D(ORIGIN, Z, 2.0)
    assert close(circle.tangent(0.0), Y)
    assert circle.curvature(0.0) == pytest.approx(0.5)
    assert close(circle.curvature_center(1.3), ORIGIN)


def test_line_is_straight():
    line = Line3D(ORIGIN, X)
    assert close(line.tangent(5.0, Tolerance.default()), X)
    assert line.curvature(5.0) == 0.0
    assert line.curvature_center(5.0) is None


def test_conic_curvatures_match_closed_forms():
    ellipse = Ellipse3D(ORIGIN, Z, X, 3.0, 2.0)
    # At the end of the major axis the curvature is a / b².
    assert ellipse.curvature(0.0) == pytest.approx(3.0 / 4.0)
    assert close(ellipse.curvature_center(0.0), Point3D(3.0 - 4.0 / 3.0, 0.0, 0.0))
    parabola = Parabola3D(ORIGIN, Z, X, 0.5)
    # At the apex the radius of curvature is twice the focal length.
    assert parabola.curvature(0.0) == pytest.approx(1.0)
    hyperbola = Hyperbola3D(ORIGIN, Z, X, 3.0, 2.0)
    # At the vertex the curvature is a / b².
    assert hyperbola.curvature(0.0) == pytest.approx(3.0 / 4.0)
    for curve in [ellipse, parabola, hyperbola]:
        tangent = curve.tangent(0.7)
        assert tangent.magnitude() == pytest.approx(1.0)
        assert tangent.cross(curve.eval_d1(0.7)[1]).magnitude() < 1e-12


def test_bspline_cusp_has_no_tangent_or_curvature():
    # A quadratic whose two middle poles coincide comes to a cusp at u = 1.
    poles = [
        Point3D(0.0, 0.0, 0.0),
        Point3D(1.0, 1.0, 0.0),
        Point3D(1.0, 1.0, 0.0),
        Point3D(2.0, 0.0, 0.0),
    ]
    cusp = BSplineCurve3D(2, poles, [0.0, 1.0, 2.0], [3, 1, 3], False)
    _, speed = cusp.eval_d1(1.0)
    assert speed.magnitude() == 0.0
    assert cusp.tangent(1.0) is None
    assert cusp.curvature(1.0) is None
    assert cusp.curvature_center(1.0) is None
    assert close(cusp.tangent(0.5), Vector3D(1.0, 1.0, 0.0).normalized())


def test_rescaling_the_parameter_keeps_the_answers():
    poles = [Point3D(0.0, 0.0, 0.0), Point3D(1.0, 2.0, 0.0), Point3D(3.0, 0.0, 1.0)]
    unit = BSplineCurve3D(2, poles, [0.0, 1.0], [3, 3], False)
    wide = BSplineCurve3D(2, poles, [0.0, 1000.0], [3, 3], False)
    assert close(unit.tangent(0.3), wide.tangent(300.0))
    assert unit.curvature(0.3) == pytest.approx(wide.curvature(300.0))
    assert close(unit.curvature_center(0.3), wide.curvature_center(300.0), abs=1e-9)


def test_trimmed_curve_local_properties_follow_the_basis():
    circle = Circle3D(ORIGIN, Z, 2.0)
    arc = TrimmedCurve3D(circle, 0.5, 2.0)
    assert close(arc.tangent(1.0), circle.tangent(1.0))
    assert arc.curvature(1.0) == pytest.approx(0.5)
    assert close(arc.curvature_center(1.0), ORIGIN)


def test_curve_local_properties_at_nan_are_none():
    circle = Circle3D(ORIGIN, Z, 2.0)
    assert circle.tangent(math.nan) is None
    assert circle.curvature(math.nan) is None
    assert circle.curvature_center(math.nan) is None


# ---------------------------------------------------------------------------
# Closedness and continuity
# ---------------------------------------------------------------------------


def test_curve_closedness():
    assert Circle3D(ORIGIN, Z, 1.0).is_closed()
    assert not Line3D(ORIGIN, X).is_closed(Tolerance.default())
    assert not kinked_polyline().is_closed()
    assert periodic_quadratic().is_closed()
    # A clamped B-spline whose ends coincide is closed without a period.
    loop = BSplineCurve3D(
        1,
        [Point3D(0.0, 0.0, 0.0), Point3D(1.0, 0.0, 0.0), Point3D(0.0, 1.0, 0.0),
         Point3D(0.0, 0.0, 0.0)],
        [0.0, 1.0, 2.0, 3.0], [2, 1, 1, 2], False,
    )
    assert loop.is_closed() and not loop.is_periodic()
    circle = Circle3D(ORIGIN, Z, 1.0)
    assert TrimmedCurve3D(circle, 1.0, 1.0 + 2.0 * math.pi).is_closed()
    assert not TrimmedCurve3D(circle, 1.0, 2.0).is_closed()


def test_curve_continuity():
    assert not hasattr(geomcore, "math.inf")
    smooth = Circle3D(ORIGIN, Z, 1.0).continuity(0.0)
    assert smooth == math.inf and isinstance(smooth, float)
    assert Line3D(ORIGIN, X).min_continuity() == math.inf
    kinked = kinked_polyline()
    kink = kinked.continuity(1.0)
    assert kink == 0 and type(kink) is int
    # A finite order compares below smooth, so plain comparisons and
    # min() read naturally.
    assert kinked.continuity(0.5) >= 2 and not kink >= 1
    assert min(kinked.continuity(u) for u in (0.5, 1.0, 1.5)) == 0
    assert kinked.continuity(1.0 + 1e-12) == 0  # within tol.parametric
    assert kinked.continuity(0.5) == math.inf
    # The ends of a non-periodic curve report smooth.
    assert kinked.continuity(0.0) == math.inf
    assert kinked.continuity(2.0) == math.inf
    assert kinked.min_continuity() == 0
    assert kinked.continuity(math.nan) == math.inf
    # A periodic curve's seam is a knot, at either bound.
    periodic = periodic_quadratic()
    assert periodic.continuity(0.0) == 1
    assert periodic.continuity(4.0) == 1
    assert periodic.min_continuity() == 1
    # A trimmed curve reports its basis inside the range, and only the
    # knots inside the range in its minimum.
    assert TrimmedCurve3D(kinked, 0.5, 1.5).continuity(1.0) == 0
    assert TrimmedCurve3D(kinked, 0.2, 0.8).min_continuity() == math.inf


def test_surface_closedness_and_continuity():
    cylinder = Cylinder(ORIGIN, Z, 1.0)
    assert cylinder.is_u_closed() and not cylinder.is_v_closed()
    torus = Torus(ORIGIN, Z, 3.0, 1.0)
    assert torus.is_u_closed(Tolerance.default()) and torus.is_v_closed()
    plane = Plane(ORIGIN, Z)
    assert not plane.is_u_closed() and not plane.is_v_closed()
    for surface in [plane, cylinder, Cone(ORIGIN, Z, 0.4, 2.0), Sphere(ORIGIN, 1.0), torus]:
        assert surface.continuity_u(0.3) == math.inf
        assert surface.continuity_v(0.3) == math.inf
        assert surface.min_continuity_u() == math.inf
        assert surface.min_continuity_v() == math.inf
    creased = creased_surface()
    assert creased.continuity_u(1.0) == 0
    assert creased.continuity_u(0.5, Tolerance.default()) == math.inf
    assert creased.continuity_v(0.5) == math.inf
    assert creased.min_continuity_u() == 0
    assert creased.min_continuity_v() == math.inf
    assert not creased.is_u_closed() and not creased.is_v_closed()
    # Coinciding first and last pole rows close a clamped surface in u.
    def row(x):
        return [Point3D(x, 0.0, 0.0), Point3D(x, 1.0, 0.0)]

    folded = BSplineSurface(
        1, 1, [row(0.0), row(1.0), row(0.0)],
        [0.0, 1.0, 2.0], [2, 1, 2], [0.0, 1.0], [2, 2], False, False,
    )
    assert folded.is_u_closed() and not folded.is_v_closed()
    assert not folded.is_u_periodic()


# ---------------------------------------------------------------------------
# Surface curvatures
# ---------------------------------------------------------------------------


def test_cylinder_curvatures_match_closed_form():
    cylinder = Cylinder(ORIGIN, Z, 2.0)
    c = cylinder.curvatures(0.0, 3.0)
    assert isinstance(c, SurfaceCurvatures)
    assert (c.min_curvature, c.max_curvature) == (-0.5, 0.0)
    assert (c.mean, c.gaussian) == (-0.25, 0.0)
    assert c.min_direction == Y
    assert c.max_direction == Z
    assert c == cylinder.curvatures(0.0, 3.0, Tolerance.default())
    assert repr(c).startswith("SurfaceCurvatures(min_curvature=-0.5, max_curvature=0.0,")


@pytest.mark.parametrize(
    "surface",
    [
        Sphere(ORIGIN, 2.0),
        Torus(ORIGIN, Z, 3.0, 1.0),
        Cone(ORIGIN, Z, 0.4, 2.0),
        Plane(ORIGIN, Z),
    ],
)
def test_principal_directions_are_a_perpendicular_unit_tangent_pair(surface):
    c = surface.curvatures(0.3, 0.2)
    _, (su, sv) = surface.eval_d1(0.3, 0.2)
    normal = su.cross(sv).normalized()
    for direction in (c.min_direction, c.max_direction):
        assert direction.magnitude() == pytest.approx(1.0)
        assert abs(direction.dot(normal)) < 1e-12
    assert abs(c.min_direction.dot(c.max_direction)) < 1e-12
    assert c.min_curvature <= c.max_curvature
    assert c.mean == pytest.approx((c.min_curvature + c.max_curvature) / 2)
    assert c.gaussian == pytest.approx(c.min_curvature * c.max_curvature)


def test_sphere_bends_away_from_its_outward_normal():
    c = Sphere(ORIGIN, 2.0).curvatures(0.3, 0.2)
    assert c.min_curvature == pytest.approx(-0.5)
    assert c.max_curvature == pytest.approx(-0.5)


def test_curvatures_are_none_where_the_normal_degenerates():
    sphere = Sphere(ORIGIN, 2.0)
    assert sphere.curvatures(0.3, math.pi / 2) is None
    cone = Cone(ORIGIN, Z, 0.4, 2.0)
    assert cone.curvatures(0.0, -2.0 / math.sin(0.4)) is None
    assert sphere.curvatures(math.nan, 0.2) is None
    # A B-spline whose first pole row collapses to one point.
    apex = Point3D(0.0, 0.0, 1.0)
    collapsed = BSplineSurface(
        1, 1,
        [[apex, apex], [Point3D(1.0, 0.0, 0.0), Point3D(0.0, 1.0, 0.0)]],
        [0.0, 1.0], [2, 2], [0.0, 1.0], [2, 2], False, False,
    )
    assert collapsed.curvatures(0.0, 0.5) is None
    assert collapsed.curvatures(0.5, 0.5) is not None


# ---------------------------------------------------------------------------
# Surface reversal, frames and orientation
# ---------------------------------------------------------------------------

ELEMENTARY = [
    Plane(ORIGIN, Z),
    Cylinder(ORIGIN, Z, 1.0),
    Cone(ORIGIN, Z, 0.4, 2.0),
    Sphere(ORIGIN, 2.0),
    Torus(ORIGIN, Z, 3.0, 1.0),
]


def normal_at(surface, u, v):
    return surface.normal() if isinstance(surface, Plane) else surface.normal(u, v)


@pytest.mark.parametrize("surface", ELEMENTARY, ids=lambda s: type(s).__name__)
def test_reversal_keeps_the_points_and_flips_the_normal(surface):
    u, v = 0.3, 0.2
    point = surface.eval_point(u, v)
    reversed_u = surface.reversed_u()
    reversed_v = surface.reversed_v()
    assert type(reversed_u) is type(surface)
    u_there = surface.u_reversed_parameter(u)
    v_there = surface.v_reversed_parameter(v)
    assert reversed_u.eval_point(u_there, v).distance(point) < 1e-12
    assert reversed_v.eval_point(u, v_there).distance(point) < 1e-12
    assert not reversed_u.frame().is_direct()
    assert close(normal_at(reversed_u, u_there, v), -normal_at(surface, u, v))
    assert close(normal_at(reversed_v, u, v_there), -normal_at(surface, u, v))
    # Reversing twice gives back the surface.
    assert reversed_u.reversed_u() == surface


def test_bspline_surface_reversed_parameters():
    surface = creased_surface()
    u_there = surface.u_reversed_parameter(0.5)
    assert u_there == 1.5
    assert surface.v_reversed_parameter(0.25) == 0.75
    point = surface.eval_point(0.5, 0.25)
    assert surface.reversed_u().eval_point(u_there, 0.25).distance(point) < 1e-12
    assert surface.reversed_v().eval_point(0.5, 0.75).distance(point) < 1e-12


def test_reversed_surface_curvatures_change_sign():
    sphere = Sphere(ORIGIN, 2.0)
    reversed_sphere = sphere.reversed_u()
    c = reversed_sphere.curvatures(sphere.u_reversed_parameter(0.3), 0.2)
    assert c.min_curvature == pytest.approx(0.5)


def test_frame_from_axes_accepts_either_handedness():
    right = Frame3D.from_axes(ORIGIN, X, Y, Z)
    assert right.is_direct()
    left = Frame3D.from_axes(ORIGIN, X, -Y, Z)
    assert not left.is_direct()
    assert close(left.y_direction(), -Y)
    world = Frame3D.world()
    assert world.is_direct()
    for flipped in (world.flipped_x(), world.flipped_y(), world.flipped_z()):
        assert not flipped.is_direct()
    assert close(world.flipped_x().x_direction(), -X)
    assert close(world.flipped_y().y_direction(), -Y)
    assert close(world.flipped_z().z_direction(), -Z)
    assert world.flipped_z().flipped_z() == world


@pytest.mark.parametrize(
    "axes",
    [
        (Vector3D(math.nan, 0.0, 0.0), Y, Z),
        (X, Vector3D(0.0, math.inf, 0.0), Z),
        (Vector3D(0.0, 0.0, 0.0), Y, Z),
        (X, Vector3D(1.0, 1.0, 0.0), Z),
    ],
)
def test_frame_from_axes_rejects_bad_directions(axes):
    with pytest.raises(FrameConstructionError):
        Frame3D.from_axes(ORIGIN, *axes)
    with pytest.raises(FrameConstructionError):
        Frame3D.from_axes(Point3D(math.nan, 0.0, 0.0), X, Y, Z)


def test_plane_and_circle_normals_follow_orientation():
    left = Frame3D.from_axes(ORIGIN, X, -Y, Z)
    assert close(Plane.from_frame(left).normal(), -Z)
    assert close(Circle3D.from_frame(left, 1.0).normal(), -Z)
    assert close(Plane(ORIGIN, Z).normal(), Z)


def test_mirrored_surface_maps_its_parametrization_and_flips_its_normal():
    sphere = Sphere(Point3D(1.0, 2.0, 0.0), 2.0)
    mirror = Transform.mirror_plane(Frame3D.world())
    image = sphere.transformed(mirror)
    u, v = 0.3, 0.2
    assert image.eval_point(u, v).distance(mirror.apply_point(sphere.eval_point(u, v))) < 1e-12
    assert not image.frame().is_direct()
    # The image's normal is the opposite of the mirrored normal: inward.
    assert close(image.normal(u, v), -mirror.apply_vector(sphere.normal(u, v)))
    assert close(image.reversed_u().normal(image.u_reversed_parameter(u), v),
                 mirror.apply_vector(sphere.normal(u, v)))
    circle = Circle3D(ORIGIN, Z, 1.0).transformed(Transform.mirror_axis(Axis3D(ORIGIN, X)))
    assert close(circle.normal(), -Z)


# ---------------------------------------------------------------------------
# One-sided evaluation
# ---------------------------------------------------------------------------


def test_side_is_a_hashable_named_enum():
    assert Side.BELOW != Side.ABOVE
    assert Side.BELOW.name == "BELOW" and Side.ABOVE.name == "ABOVE"
    assert repr(Side.ABOVE) == "Side.ABOVE"
    assert {Side.BELOW: 1, Side.ABOVE: 2}[Side.ABOVE] == 2
    assert geomcore.Side is Side


def test_curve_sided_evaluation_at_a_kink():
    curve = kinked_polyline()
    below = curve.eval_d1_sided(1.0, Side.BELOW)
    above = curve.eval_d1_sided(1.0, Side.ABOVE)
    assert below[0] == above[0] == Point3D(1.0, 0.0, 0.0)
    assert below[1] == X and above[1] == Y
    # Plain evaluation takes the above side.
    assert curve.eval_d1(1.0)[1] == above[1]
    assert curve.eval_derivative_sided(1.0, 1, Side.BELOW) == X
    point, first, second = curve.eval_d2_sided(1.0, Side.BELOW)
    assert (point, first, second) == (Point3D(1.0, 0.0, 0.0), X, Vector3D(0.0, 0.0, 0.0))
    # Strictly inside a span both sides agree.
    assert curve.eval_d1_sided(0.5, Side.BELOW) == curve.eval_d1_sided(0.5, Side.ABOVE)


def test_curve_sided_evaluation_wraps_a_periodic_seam():
    curve = periodic_quadratic()
    last_span = curve.eval_derivative(3.5, 2)  # constant on a quadratic span
    first_span = curve.eval_derivative(0.5, 2)
    assert first_span != last_span
    for seam in (0.0, 4.0):
        assert close(curve.eval_derivative_sided(seam, 2, Side.BELOW), last_span)
        assert close(curve.eval_derivative_sided(seam, 2, Side.ABOVE), first_span)
        assert close(curve.eval_d2_sided(seam, Side.BELOW)[2], last_span)


def test_curve_sided_evaluation_rejects_unsupported_orders():
    curve = kinked_polyline()
    with pytest.raises(ValueError, match="order must be >= 1"):
        curve.eval_derivative_sided(1.0, 0, Side.BELOW)
    with pytest.raises(ValueError, match="at most 2"):
        curve.eval_derivative_sided(1.0, 3, Side.BELOW)


def test_sided_evaluation_at_nan_gives_nan():
    point, first = kinked_polyline().eval_d1_sided(math.nan, Side.BELOW)
    assert all(math.isnan(c) for c in coordinates(point) + coordinates(first))
    point, (du, dv) = creased_surface().eval_d1_sided(math.nan, 0.5, Side.BELOW, Side.ABOVE)
    assert all(math.isnan(c) for c in coordinates(point))


def test_surface_sided_evaluation_across_a_crease():
    surface = creased_surface()
    below = surface.eval_d1_sided(1.0, 0.5, Side.BELOW, Side.ABOVE)
    above = surface.eval_d1_sided(1.0, 0.5, Side.ABOVE, Side.ABOVE)
    assert below[1][0] == Vector3D(1.0, 0.0, 0.0)
    assert above[1][0] == Vector3D(1.0, 0.0, 0.5)
    assert surface.eval_d1(1.0, 0.5)[1][0] == above[1][0]
    assert surface.eval_derivative_sided(1.0, 0.5, 1, 1, Side.BELOW, Side.ABOVE) == Vector3D(
        0.0, 0.0, 0.0
    )
    assert surface.eval_derivative_sided(1.0, 0.5, 1, 1, Side.ABOVE, Side.ABOVE) == Vector3D(
        0.0, 0.0, 1.0
    )
    point, (du, dv), (duu, duv, dvv) = surface.eval_d2_sided(1.0, 0.5, Side.BELOW, Side.BELOW)
    assert du == Vector3D(1.0, 0.0, 0.0) and duv == Vector3D(0.0, 0.0, 0.0)
    with pytest.raises(ValueError, match="du \\+ dv must be >= 1"):
        surface.eval_derivative_sided(1.0, 0.5, 0, 0, Side.BELOW, Side.BELOW)
    with pytest.raises(ValueError, match="du \\+ dv must be <= 2"):
        surface.eval_derivative_sided(1.0, 0.5, 2, 1, Side.BELOW, Side.BELOW)
