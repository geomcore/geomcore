"""Behavioral tests for the geomcore Python bindings."""

import math
from typing import get_args

import pytest

import geomcore
from geomcore import (
    Axis,
    BSplineCurve,
    BSplineCurve2D,
    BSplineSurface,
    BezierCurve,
    BezierSurface,
    BoundingBox,
    Circle,
    CurveApproximation,
    Cylinder,
    Ellipse,
    Ellipse2D,
    Frame,
    Hyperbola,
    Hyperbola2D,
    Line,
    Line2D,
    Parabola,
    Parabola2D,
    Plane,
    Point,
    Point2D,
    Sphere,
    SurfaceApproximation,
    SurfaceSurfaceIntersection,
    Tolerance,
    Transform,
    Vector,
    Vector2D,
)


def test_readme_snippet():
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    point = circle.point_at(math.pi / 4)
    assert point.x == pytest.approx(2.0 * math.cos(math.pi / 4))
    assert point.y == pytest.approx(2.0 * math.sin(math.pi / 4))
    assert point.z == pytest.approx(0.0)


def test_eval_points_bulk():
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    params = [i / 100 * math.tau for i in range(100)]
    points = circle.points_at(params)
    assert len(points) == 100
    for u, p in zip(params, points):
        assert p.x == pytest.approx(2.0 * math.cos(u))
        assert p.y == pytest.approx(2.0 * math.sin(u))


def test_constructor_error_raises_value_error():
    with pytest.raises(ValueError) as excinfo:
        Circle(Point.origin(), Vector.unit_z(), -1.0)
    assert str(excinfo.value)


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_non_finite_inputs_raise_value_error(bad):
    from geomcore import Parabola, Plane

    with pytest.raises(ValueError, match="NaN or infinite"):
        Circle(Point.origin(), Vector.unit_z(), bad)
    with pytest.raises(ValueError, match="NaN or infinite"):
        Sphere(Point(bad, 0.0, 0.0), 1.0)
    with pytest.raises(ValueError, match="NaN or infinite"):
        Plane(Point.origin(), Vector(bad, 0.0, 1.0))
    with pytest.raises(ValueError, match="NaN or infinite"):
        Parabola(Point.origin(), Vector.unit_z(), Vector.unit_x(), bad)


def test_zero_size_conics_raise_value_error():
    from geomcore import Ellipse, Hyperbola, Parabola

    o, z, x = Point.origin(), Vector.unit_z(), Vector.unit_x()
    with pytest.raises(ValueError, match="not positive"):
        Parabola(o, z, x, 0.0)
    with pytest.raises(ValueError, match="not positive"):
        Ellipse(o, z, x, 3.0, 0.0)
    with pytest.raises(ValueError, match="not positive"):
        Hyperbola(o, z, x, 2.0, 0.0)


def test_init_constructs():
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    assert circle.radius() == 2.0
    # `new` is no longer a constructor: only __init__ builds objects.
    assert not hasattr(Circle, "new")


def test_curve_3d_to_uv_cylinder():
    # Coaxial circle on a cylinder of the same radius: the UV curve is the
    # horizontal line v = 0 in (u, v) space.
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    uv_curve = cylinder.curve_3d_to_uv(circle)
    assert isinstance(uv_curve, Line2D)
    origin = uv_curve.origin()
    direction = uv_curve.direction()
    assert origin.x == pytest.approx(0.0)
    assert origin.y == pytest.approx(0.0)
    assert direction.x == pytest.approx(1.0)
    assert direction.y == pytest.approx(0.0)


def test_curve_3d_to_uv_takes_a_tolerance():
    # A radius 1e-4 too wide lies on the cylinder at confusion 1e-3 only.
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 1.0)
    circle = Circle(Point(0.0, 0.0, 2.0), Vector.unit_z(), 1.0 + 1e-4)
    loose = Tolerance(confusion=1e-3)
    assert isinstance(cylinder.curve_3d_to_uv(circle, loose), Line2D)
    assert isinstance(cylinder.curve_3d_to_uv(circle, tol=loose), Line2D)
    with pytest.raises(ValueError):
        cylinder.curve_3d_to_uv(circle)


def test_curve_3d_to_uv_not_analytic_raises():
    sphere = Sphere(Point.origin(), 2.0)
    line = geomcore.Line(Point(2.0, 0.0, 0.0), Vector.unit_z())
    with pytest.raises(ValueError):
        sphere.curve_3d_to_uv(line)


def test_curve_3d_to_uv_plane_returns_the_new_2d_types():
    plane = Plane(Point(0.0, 0.0, 1.0), Vector.unit_z())
    center = Point(2.0, 1.0, 1.0)

    ellipse = Ellipse(center, Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    uv_curve = plane.curve_3d_to_uv(ellipse)
    assert isinstance(uv_curve, Ellipse2D)
    assert uv_curve.major_radius() == pytest.approx(3.0)
    assert uv_curve.minor_radius() == pytest.approx(1.5)
    assert uv_curve.center().x == pytest.approx(2.0)
    assert uv_curve.center().y == pytest.approx(1.0)

    parabola = Parabola(center, Vector.unit_z(), Vector.unit_x(), 0.5)
    uv_curve = plane.curve_3d_to_uv(parabola)
    assert isinstance(uv_curve, Parabola2D)
    assert uv_curve.focal_distance() == pytest.approx(0.5)

    hyperbola = Hyperbola(center, Vector.unit_z(), Vector.unit_x(), 2.0, 1.0)
    uv_curve = plane.curve_3d_to_uv(hyperbola)
    assert isinstance(uv_curve, Hyperbola2D)
    assert uv_curve.major_radius() == pytest.approx(2.0)

    # Each uv_curve maps back onto the 3D curve: surface(q(t)) == curve(t).
    for curve in (ellipse, parabola, hyperbola):
        uv_curve = plane.curve_3d_to_uv(curve)
        for t in (0.0, 0.7, 1.9):
            q = uv_curve.point_at(t)
            expected = curve.point_at(t)
            assert q.x == pytest.approx(expected.x)
            assert q.y == pytest.approx(expected.y)


def test_curve_3d_to_uv_plane_bspline_returns_bspline_2d():
    plane = Plane(Point(0.0, 0.0, 1.0), Vector.unit_z())
    curve = BSplineCurve(
        2,
        [Point(0.0, 0.0, 1.0), Point(1.0, 2.0, 1.0), Point(3.0, 0.0, 1.0)],
        [0.0, 1.0],
        [3, 3],
        False,
    )
    uv_curve = plane.curve_3d_to_uv(curve)
    assert isinstance(uv_curve, BSplineCurve2D)
    assert [(p.x, p.y) for p in uv_curve.poles()] == [(0.0, 0.0), (1.0, 2.0), (3.0, 0.0)]
    assert uv_curve.knots() == [0.0, 1.0]


def test_curve_3d_to_uv_off_plane_and_curved_surfaces_raise():
    plane = Plane(Point.origin(), Vector.unit_z())
    raised = Ellipse(Point(0.0, 0.0, 1.0), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    with pytest.raises(ValueError):
        plane.curve_3d_to_uv(raised)
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    flat = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    with pytest.raises(ValueError):
        cylinder.curve_3d_to_uv(flat)


def test_curve_3d_to_uv_maps_a_bezier_curve_on_a_plane_and_back():
    plane = Plane(Point(0.0, 0.0, 1.0), Vector.unit_z())
    curve = geomcore.BezierCurve(
        2, [Point(0.0, 0.0, 1.0), Point(1.0, 2.0, 1.0), Point(3.0, 0.0, 1.0)]
    )
    uv_curve = plane.curve_3d_to_uv(curve)
    assert isinstance(uv_curve, geomcore.BezierCurve2D)
    assert [(p.x, p.y) for p in uv_curve.poles()] == [(0.0, 0.0), (1.0, 2.0), (3.0, 0.0)]
    traced = plane.uv_curve_to_3d(uv_curve)
    for t in (0.0, 0.4, 1.0):
        assert traced.point_at(t).distance(curve.point_at(t)) < 1e-12


def test_curve_3d_to_uv_takes_a_3d_curve_only():
    plane = Plane(Point.origin(), Vector.unit_z())
    with pytest.raises(TypeError):
        plane.curve_3d_to_uv(Line2D(geomcore.Point2D(0.0, 0.0), geomcore.Vector2D(1.0, 0.0)))


def test_bspline_curve_degree_one_line():
    curve = BSplineCurve(
        1,
        [Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0)],
        [0.0, 1.0],
        [2, 2],
        False,
    )
    mid = curve.point_at(0.5)
    assert mid.x == pytest.approx(1.0)
    assert mid.y == pytest.approx(0.0)
    assert curve.degree() == 1
    assert not curve.is_periodic()


def test_transform_rotation():
    axis = Axis(Point.origin(), Vector.unit_z())
    rot = Transform.rotation(axis, math.pi / 2)
    p = rot.apply_point(Point(1.0, 0.0, 0.0))
    assert p.x == pytest.approx(0.0)
    assert p.y == pytest.approx(1.0)


def test_transform_composition():
    t = Transform.translation(Vector(1.0, 0.0, 0.0))
    rot = Transform.rotation(Axis(Point.origin(), Vector.unit_z()), math.pi / 2)
    composed = t.then(rot)
    p = composed.apply_point(Point.origin())
    assert p.x == pytest.approx(0.0)
    assert p.y == pytest.approx(1.0)


def test_sphere_parameters_round_trip():
    sphere = Sphere(Point.origin(), 3.0)
    p = sphere.point_at(0.7, 0.4)
    u, v = sphere.parameters_of(p)
    assert u == pytest.approx(0.7)
    assert v == pytest.approx(0.4)


def test_surface_normal_is_unit_and_recovers_at_singular_points():
    sphere = Sphere(Point.origin(), 3.0)
    assert sphere.normal(0.0, 0.0).components() == pytest.approx((1.0, 0.0, 0.0))
    pole = sphere.normal(0.3, math.pi / 2, Tolerance.default())
    assert pole.components() == pytest.approx((0.0, 0.0, 1.0), abs=1e-12)
    cone = geomcore.Cone(Point.origin(), Vector.unit_z(), 0.4, 2.0)
    v_apex = -2.0 / math.sin(0.4)
    # Each meridian has its own limit at the apex, so there is no normal.
    assert cone.normal(0.0, v_apex) is None
    beside = cone.normal(0.0, v_apex + 1e-9)
    assert beside.components() == pytest.approx(
        (math.cos(0.4), 0.0, -math.sin(0.4)), abs=1e-9
    )
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 1.0)
    assert cylinder.normal(0.0, 5.0).components() == pytest.approx((1.0, 0.0, 0.0))
    assert sphere.normal(float("nan"), 0.0) is None



def test_plane_normal_takes_the_parameters_every_surface_takes():
    # #190: the plane's own zero-argument `normal()` hid the shared
    # `normal(u, v, tol)`, so a loop over surfaces raised on the plane.
    plane = Plane(Point(0.0, 0.0, 1.0), Vector(0.0, 0.0, 2.0))
    surfaces = [plane, Cylinder(Point.origin(), Vector.unit_z(), 1.0), Sphere(Point.origin(), 3.0)]
    for surface in surfaces:
        assert surface.normal(0.3, 0.2) is not None
        assert surface.normal(0.3, 0.2, Tolerance.default()) is not None
    # The normal is the same at every point; the plane's axis carries it.
    assert plane.normal(5.0, -7.0) == Vector.unit_z()
    assert plane.axis().origin() == Point(0.0, 0.0, 1.0)
    assert plane.axis().direction() == Vector.unit_z()
    reversed_plane = plane.reversed_u()
    assert reversed_plane.axis().direction() == -Vector.unit_z()
    assert reversed_plane.normal(0.3, 0.2) == -Vector.unit_z()


def test_frame3_accessors():
    frame = Frame(Point.origin(), Vector.unit_z(), Vector.unit_x())
    assert frame.z_direction().components() == pytest.approx((0.0, 0.0, 1.0))
    assert frame.x_direction().components() == pytest.approx((1.0, 0.0, 0.0))


def test_docstrings_present():
    assert Circle.__doc__
    assert Circle.point_at.__doc__


def test_tolerance_defaults():
    tol = Tolerance.default()
    assert tol.confusion == 1e-7
    assert tol.angular == 1e-12
    assert tol.parametric == 1e-9
    assert Tolerance().confusion == tol.confusion
    custom = Tolerance(confusion=1e-3)
    assert custom.confusion == 1e-3
    assert custom.angular == tol.angular
    # Models fit in a cube of MODEL_SIZE; angles derive from distances.
    assert Tolerance.MODEL_SIZE == 1e5
    assert tol.angular_for_extent(10.0) == 1e-8
    assert tol.angular_for_extent(1e9) == tol.angular


def test_points_coincide_within_confusion():
    p = Point(1.0, 2.0, 3.0)
    assert p.coincides(Point(1.0, 2.0, 3.0 + 5e-8))
    assert not p.coincides(Point(1.0, 2.0, 3.0 + 5e-7))
    near = Point(1.0, 2.0005, 3.0)
    assert not p.coincides(near)
    assert p.coincides(near, Tolerance(confusion=1e-3))
    assert p.coincides(near, tol=Tolerance(confusion=1e-3))
    q = Point2D(1.0, 2.0)
    assert q.coincides(Point2D(1.0 + 5e-8, 2.0))
    assert not q.coincides(Point2D(1.0 + 5e-7, 2.0))


def test_a_number_as_tol_is_its_confusion():
    p = Point(1.0, 2.0, 3.0)
    near = Point(1.0, 2.0005, 3.0)
    assert p.coincides(near, 1e-3)
    assert p.coincides(near, tol=1)
    assert not p.coincides(near, 1e-4)
    with pytest.raises(TypeError, match="tol must be a Tolerance or a number"):
        p.coincides(near, "1e-3")  # type: ignore[arg-type]


def test_confusion_steers_the_size_of_a_fitted_intersection_curve():
    def wave(phase: float, height: float) -> BSplineSurface:
        rows = [
            [Point(x / 2, y / 2, height + math.sin(x / 2 + phase) * math.cos(y / 2)) for y in range(9)]
            for x in range(9)
        ]
        return BSplineSurface.approximate(rows, 3, 3, 6, 6)

    first, second = wave(0.0, 0.0), wave(1.0, 0.2)
    loose = first.intersect(second, 1e-3)
    assert loose == first.intersect(second, Tolerance(confusion=1e-3))
    fine = first.intersect(second)

    def poles(result: SurfaceSurfaceIntersection) -> int:
        curves = [piece.curve for piece in result.curves]
        assert all(isinstance(curve, BSplineCurve) for curve in curves)
        return sum(len(curve.poles()) for curve in curves if isinstance(curve, BSplineCurve))

    assert poles(loose) < poles(fine) / 5


def _bspline_curve(end):
    return BSplineCurve(1, [Point.origin(), end], [0.0, 1.0], [2, 2], False)


def _bspline_surface(corner_z):
    poles = [
        [Point(0.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)],
        [Point(1.0, 0.0, 0.0), Point(1.0, 1.0, corner_z)],
    ]
    return geomcore.BSplineSurface(
        1, 1, poles, [0.0, 1.0], [2, 2], [0.0, 1.0], [2, 2], False, False
    )


_O, _X, _Z = Point.origin(), Vector.unit_x(), Vector.unit_z()

# For each class: a factory building one value, and a value that differs from it.
VALUE_TYPES = [
    (lambda: Point(1.0, 2.0, 3.0), Point(1.0, 2.0, 3.0 + 1e-12)),
    (lambda: Point2D(1.0, 2.0), Point2D(1.0, 2.5)),
    (lambda: Vector(1.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0)),
    (lambda: Vector2D(1.0, 0.0), Vector2D(0.0, 1.0)),
    (lambda: Axis(_O, _Z), Axis(_O, _X)),
    (lambda: Frame(_O, _Z, _X), Frame(Point(0.0, 0.0, 1.0), _Z, _X)),
    (lambda: Tolerance(), Tolerance(confusion=1e-3)),
    (lambda: Transform.translation(_X), Transform.identity()),
    (lambda: BoundingBox(_O, Point(1.0, 1.0, 1.0)), BoundingBox(_O, Point(2.0, 1.0, 1.0))),
    (lambda: geomcore.Plane(_O, _Z), geomcore.Plane(_O, _X)),
    (lambda: Cylinder(_O, _Z, 1.0), Cylinder(_O, _Z, 2.0)),
    (lambda: geomcore.Cone(_O, _Z, 0.3, 1.0), geomcore.Cone(_O, _Z, 0.4, 1.0)),
    (lambda: Sphere(_O, 1.0), Sphere(_O, 2.0)),
    (lambda: geomcore.Torus(_O, _Z, 3.0, 1.0), geomcore.Torus(_O, _Z, 3.0, 0.5)),
    (lambda: _bspline_surface(0.0), _bspline_surface(1.0)),
    (
        lambda: geomcore.ExtrusionSurface(Circle(_O, _Z, 2.0), _Z),
        geomcore.ExtrusionSurface(Circle(_O, _Z, 2.0), Vector(0.0, 1.0, 1.0)),
    ),
    (
        lambda: geomcore.RevolutionSurface(geomcore.Line(Point(1.0, 0.0, 0.0), _Z), Axis(_O, _Z)),
        geomcore.RevolutionSurface(geomcore.Line(Point(2.0, 0.0, 0.0), _Z), Axis(_O, _Z)),
    ),
    (lambda: geomcore.Line(_O, _X), geomcore.Line(_O, _Z)),
    (lambda: Circle(_O, _Z, 2.0), Circle(_O, _Z, 2.5)),
    (
        lambda: geomcore.Ellipse(_O, _Z, _X, 3.0, 1.0),
        geomcore.Ellipse(_O, _Z, _X, 3.0, 2.0),
    ),
    (lambda: geomcore.Parabola(_O, _Z, _X, 1.0), geomcore.Parabola(_O, _Z, _X, 2.0)),
    (
        lambda: geomcore.Hyperbola(_O, _Z, _X, 3.0, 1.0),
        geomcore.Hyperbola(_O, _Z, _X, 3.0, 2.0),
    ),
    (lambda: _bspline_curve(Point(1.0, 0.0, 0.0)), _bspline_curve(Point(2.0, 0.0, 0.0))),
    (lambda: Line2D(Point2D.origin(), Vector2D(1.0, 0.0)), Line2D(Point2D.origin(), Vector2D(0.0, 1.0))),
    (lambda: geomcore.Circle2D(Point2D.origin(), 1.0), geomcore.Circle2D(Point2D.origin(), 2.0)),
    (
        lambda: geomcore.TrimmedCurve.segment(_O, Point(1.0, 0.0, 0.0)),
        geomcore.TrimmedCurve.segment(_O, Point(2.0, 0.0, 0.0)),
    ),
    (
        lambda: geomcore.TrimmedCurve2D.segment(Point2D.origin(), Point2D(1.0, 0.0)),
        geomcore.TrimmedCurve2D.segment(Point2D.origin(), Point2D(2.0, 0.0)),
    ),
    (lambda: Circle(_O, _Z, 2.0).discretize(0.1), Circle(_O, _Z, 2.0).discretize(0.01)),
    (
        lambda: geomcore.Circle2D(Point2D.origin(), 1.0).discretize(0.1),
        geomcore.Circle2D(Point2D.origin(), 1.0).discretize(0.01),
    ),
]


@pytest.mark.parametrize(
    "make, different", VALUE_TYPES, ids=[type(d).__name__ for _, d in VALUE_TYPES]
)
def test_equality_compares_values_exactly(make, different):
    value = make()
    assert value == make()
    assert not value != make()
    assert value != different
    assert value != "not a geometry value"


@pytest.mark.parametrize(
    "make, different", VALUE_TYPES, ids=[type(d).__name__ for _, d in VALUE_TYPES]
)
def test_value_types_are_unhashable(make, different):
    with pytest.raises(TypeError, match="unhashable"):
        hash(make())


def test_equal_points_from_different_computations():
    curve = _bspline_curve(Point(2.0, 0.0, 0.0))
    assert curve.point_at(0.5) == Point(1.0, 0.0, 0.0)
    assert curve.point_at(0.0) == Point.origin()


def test_contains_curve():
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    assert circle.contains(Point(2.0, 0.0, 0.0))
    assert circle.contains(circle.point_at(1.0))
    assert not circle.contains(Point(3.0, 0.0, 0.0))
    assert not circle.contains(Point.origin())
    # Explicit tolerance widens the acceptance band.
    near = Point(2.0 + 1e-6, 0.0, 0.0)
    assert not circle.contains(near)
    assert circle.contains(near, Tolerance(confusion=1e-5))


def test_contains_surface():
    from geomcore import Plane, Torus

    plane = Plane(Point.origin(), Vector.unit_z())
    assert plane.contains(Point.origin())
    assert plane.contains(Point(1.0, 2.0, 0.0))
    assert not plane.contains(Point(0.0, 0.0, 1.0))

    sphere = Sphere(Point.origin(), 3.0)
    assert sphere.contains(Point(3.0, 0.0, 0.0))
    assert not sphere.contains(Point.origin())

    torus = Torus(Point.origin(), Vector.unit_z(), 4.0, 1.0)
    assert torus.contains(Point(5.0, 0.0, 0.0))
    assert not torus.contains(Point.origin())


def test_project_point_curve():
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    projection = circle.project_point(Point(3.0, 0.0, 0.0))
    assert projection.parameter == pytest.approx(0.0)
    assert projection.distance == pytest.approx(1.0)
    projection = circle.project_point(Point(2.0, 0.0, 0.0))
    assert projection.distance == 0.0
    # Batch variant agrees with scalar calls.
    batch = circle.project_points([Point(3.0, 0.0, 0.0), Point(2.0, 0.0, 0.0)])
    assert (batch[0].parameter, batch[0].distance) == pytest.approx((0.0, 1.0))
    assert batch[1].distance == 0.0


def test_project_point_surface():
    from geomcore import Plane

    plane = Plane(Point.origin(), Vector.unit_z())
    projection = plane.project_point(Point(1.0, 2.0, 3.0))
    assert (projection.u, projection.v) == pytest.approx((1.0, 2.0))
    assert projection.distance == pytest.approx(3.0)
    batch = plane.project_points([Point.origin(), Point(0.0, 0.0, 2.0)])
    assert batch[0].distance == 0.0
    assert batch[1].distance == pytest.approx(2.0)


def _curves(result):
    """The curves of a surface-surface result."""
    return [piece.curve for piece in result.curves]


def test_intersect_plane_plane():
    from geomcore import Line, Plane

    xy = Plane(Point.origin(), Vector.unit_z())
    zy = Plane(Point.origin(), Vector.unit_x())
    (line,) = _curves(xy.intersect(zy))
    assert isinstance(line, Line)
    assert line.origin().x == pytest.approx(0.0)
    # Direction follows the ordered pair (n1 x n2 = Z x X = Y).
    assert line.direction().components() == pytest.approx((0.0, 1.0, 0.0))

    parallel = xy.intersect(Plane(Point(0.0, 0.0, 1.0), Vector.unit_z()))
    assert parallel.is_empty() and parallel.is_complete()

    coincident = xy.intersect(xy)
    assert coincident.is_coincident()
    assert coincident.curves == []


def test_intersect_plane_sphere():
    from geomcore import Plane

    plane = Plane(Point.origin(), Vector.unit_z())
    (circle,) = _curves(plane.intersect(Sphere(Point.origin(), 2.0)))
    assert circle.radius() == pytest.approx(2.0)

    tangent = plane.intersect(Sphere(Point(0.0, 0.0, 2.0), 2.0))
    assert tangent.curves == []
    (point,) = tangent.contact_points
    assert (point.x, point.y, point.z) == pytest.approx((0.0, 0.0, 0.0))

    assert plane.intersect(Sphere(Point(0.0, 0.0, 5.0), 2.0)).is_empty()


def test_intersect_sphere_sphere():
    s1 = Sphere(Point.origin(), 2.0)
    (circle,) = _curves(s1.intersect(Sphere(Point(3.0, 0.0, 0.0), 2.0)))
    assert circle.center().x == pytest.approx(1.5)

    (point,) = s1.intersect(Sphere(Point(4.0, 0.0, 0.0), 2.0)).contact_points
    assert point.x == pytest.approx(2.0)

    assert s1.intersect(Sphere(Point(5.0, 0.0, 0.0), 2.0)).is_empty()
    assert s1.intersect(s1).is_coincident()


def test_intersect_plane_cylinder():
    from geomcore import Cylinder, Plane

    plane = Plane(Point.origin(), Vector.unit_z())
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    (circle,) = _curves(plane.intersect(cylinder))
    assert circle.radius() == pytest.approx(2.0)

    side = Plane(Point(0.0, 1.0, 0.0), Vector.unit_y())
    lines = _curves(side.intersect(cylinder))
    assert len(lines) == 2
    assert lines[0].direction().components() == pytest.approx((0.0, 0.0, 1.0))

    assert side.intersect(Cylinder(Point.origin(), Vector.unit_z(), 0.5)).is_empty()
    miss = Plane(Point(0.0, 3.0, 0.0), Vector.unit_y())
    assert miss.intersect(cylinder).is_empty()


def test_intersect_plane_cone():
    from geomcore import Cone, Frame, Hyperbola, Plane

    plane = Plane(Point.origin(), Vector.unit_z())
    cone = Cone.from_frame(Frame.world(), 0.4, 2.0)
    (circle,) = _curves(plane.intersect(cone))
    assert circle.radius() == pytest.approx(2.0)

    # The cone is double: a plane beyond the apex cuts the other nappe.
    below = Plane(Point(0.0, 0.0, -6.0), Vector.unit_z())
    (circle,) = _curves(below.intersect(cone))
    assert circle.radius() == pytest.approx(0.5368, abs=1e-3)

    h1, h2 = _curves(Plane(Point(3.0, 0.0, 0.0), Vector.unit_x()).intersect(cone))
    assert isinstance(h1, Hyperbola) and isinstance(h2, Hyperbola)
    assert h1.center().distance(h2.center()) < 1e-9

    apex_plane = Plane(cone.apex(), Vector.unit_z())
    (apex,) = apex_plane.intersect(cone).contact_points
    assert apex.distance(cone.apex()) < 1e-9


def test_intersect_line_surface():
    from geomcore import Cylinder, Line, Plane

    line = Line(Point.origin(), Vector.unit_x())
    hits = line.intersect(Sphere(Point.origin(), 2.0)).hits
    assert [h.curve_parameter for h in hits] == pytest.approx([-2.0, 2.0])
    assert [h.point.x for h in hits] == pytest.approx([-2.0, 2.0])

    (hit,) = Line(Point.origin(), Vector.unit_z()).intersect(
        Plane(Point.origin(), Vector.unit_z())
    ).hits
    assert hit.curve_parameter == pytest.approx(0.0)
    on_plane = line.intersect(Plane(Point.origin(), Vector.unit_z()))
    assert on_plane.hits == [] and [overlap.curve_range for overlap in on_plane.overlaps] == [(-math.inf, math.inf)]

    generator = Line(Point(2.0, 0.0, 0.0), Vector.unit_z())
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    assert [overlap.curve_range for overlap in generator.intersect(cylinder).overlaps] == [(-math.inf, math.inf)]


def test_intersect_symmetric_quadrics():
    from geomcore import Cone, Cylinder, Ellipse, Frame, IntersectionMethod

    sphere = Sphere(Point.origin(), 3.0)
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)

    c1, c2 = _curves(sphere.intersect(cylinder))
    assert c1.radius() == pytest.approx(2.0)
    assert abs(c1.center().z) == pytest.approx(math.sqrt(5.0))

    assert len(sphere.intersect(Cylinder(Point.origin(), Vector.unit_z(), 3.0)).curves) == 1
    # Off-axis, a space quartic fitted from its closed form: here one loop
    # over an arc of the cylinder's generators.
    off = Cylinder(Point(0.0, 1.0, 0.0), Vector.unit_z(), 2.0)
    result = sphere.intersect(off)
    assert result.method == IntersectionMethod.EXACT and result.is_complete()
    assert len(result.curves) == 1
    # Viviani's curve crosses itself where the sphere touches the cylinder.
    viviani = Sphere(Point.origin(), 2.0).intersect(
        Cylinder(Point(1.0, 0.0, 0.0), Vector.unit_z(), 1.0)
    )
    assert viviani.method == IntersectionMethod.EXACT and len(viviani.curves) == 1
    assert viviani.contact_points[0].distance(Point(2.0, 0.0, 0.0)) < 1e-12

    cone = Cone.from_frame(Frame.world(), 0.4, 2.0)
    r1, r2 = _curves(sphere.intersect(cone))
    assert sorted([r1.radius(), r2.radius()]) == pytest.approx([0.7746, 2.6188], abs=1e-3)
    # Sphere around the apex: one ring on each nappe.
    small = Cone.from_frame(Frame.world(), 0.4, 0.5)
    rings = sorted(_curves(sphere.intersect(small)), key=lambda c: c.center().z)
    lower, upper = rings
    assert upper.radius() == pytest.approx(1.5786, abs=1e-3)
    assert lower.radius() == pytest.approx(0.7302, abs=1e-3)
    assert lower.center().z < small.apex().z < upper.center().z

    c1 = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    c2 = Cylinder(Point(3.0, 0.0, 0.0), Vector.unit_z(), 2.0)
    l1, _ = _curves(c1.intersect(c2))
    assert l1.origin().x == pytest.approx(1.5)
    assert c1.intersect(c1).is_coincident()
    # A thinner pipe through it: two loops of the quartic, fitted from its
    # closed form, each with its curves in both cylinders' parameters.
    crossed = Cylinder(Point.origin(), Vector.unit_x(), 1.5)
    result = c1.intersect(crossed)
    assert result.method == IntersectionMethod.EXACT and result.is_complete()
    assert len(result.curves) == 2
    for piece in result.curves:
        assert isinstance(piece.curve, BSplineCurve)
        p = piece.curve.point_at(0.3)
        uv = piece.second_uv_curve.point_at(0.3)
        assert crossed.point_at(uv.x, uv.y).distance(p) < 1e-6
        assert math.hypot(p.x, p.y) == pytest.approx(2.0, abs=1e-6)
        assert math.hypot(p.y, p.z) == pytest.approx(1.5, abs=1e-6)
    # The tee of two equal pipes: two ellipses, crossing where they touch.

    tee = c1.intersect(Cylinder(Point.origin(), Vector.unit_x(), 2.0))
    assert tee.method == IntersectionMethod.EXACT
    ellipses = _curves(tee)
    assert len(ellipses) == 2 and all(isinstance(e, Ellipse) for e in ellipses)
    assert [e.major_radius() for e in ellipses] == pytest.approx([2.0 * math.sqrt(2.0)] * 2)
    touching = sorted(tee.contact_points, key=lambda p: p.y)
    assert touching[0].distance(Point(0.0, -2.0, 0.0)) < 1e-12
    assert touching[1].distance(Point(0.0, 2.0, 0.0)) < 1e-12


def test_intersect_cone_cylinder():
    from geomcore import Cone, Cylinder, Frame, IntersectionMethod

    cone = Cone.from_frame(Frame.world(), 0.4, 2.0)
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    lower, upper = sorted(_curves(cone.intersect(cylinder)), key=lambda c: c.center().z)
    assert (lower.radius(), upper.radius()) == pytest.approx((2.0, 2.0))
    assert upper.center().distance(Point.origin()) < 1e-9
    assert lower.center().z == pytest.approx(2.0 * cone.apex().z)

    off = Cylinder(Point(1.0, 0.0, 0.0), Vector.unit_z(), 2.0)
    assert cone.intersect(off).method == IntersectionMethod.NUMERIC


def test_intersect_cone_cone():
    from geomcore import Cone, Frame, IntersectionMethod

    c1 = Cone.from_frame(Frame.world(), 0.4, 2.0)
    frame2 = Frame(Point(0.0, 0.0, 1.0), Vector.unit_z(), Vector.unit_x())
    c2 = Cone.from_frame(frame2, 0.6, 1.0)
    r1, r2 = _curves(c1.intersect(c2))
    assert sorted([r1.radius(), r2.radius()]) == pytest.approx([1.1155, 4.7245], abs=1e-3)

    assert c1.intersect(c1).is_coincident()
    tilted = Cone(Point.origin(), Vector.unit_x(), 0.4, 1.5)
    assert c1.intersect(tilted).method == IntersectionMethod.NUMERIC

    # Axes meeting at the origin, both cones as far from it: two conics
    # through the two points where the cones touch the sphere of radius 2·cos 0.4.
    around = Cone(Point.origin(), Vector.unit_x(), 0.4, 2.0)
    conics = c1.intersect(around)
    assert conics.method == IntersectionMethod.EXACT and conics.is_complete()
    assert len(conics.contact_points) == 2
    for point in conics.contact_points:
        assert point.distance(Point.origin()) == pytest.approx(2.0 * math.cos(0.4))

    # One apex, axes apart: four lines through it, each a generator of both.
    narrow = Cone(Point.origin(), Vector.unit_z(), math.pi / 4.0, 0.0)
    wide = Cone(Point.origin(), Vector.unit_x(), math.pi / 3.0, 0.0)
    shared = narrow.intersect(wide)
    assert shared.method == IntersectionMethod.EXACT and shared.is_complete()
    assert shared.contact_points == [Point.origin()]
    assert len(shared.curves) == 4
    for piece in shared.curves:
        assert isinstance(piece.curve, Line) and not piece.tangent
        d = piece.curve.direction()
        assert (abs(d.x), abs(d.y), abs(d.z)) == pytest.approx((0.5, 0.5, math.sqrt(0.5)))
        uv = piece.second_uv_curve.point_at(2.0)
        assert wide.point_at(uv.x, uv.y).distance(piece.curve.point_at(2.0)) < 1e-12


def test_intersection_curves_carry_their_uv_curves():
    from geomcore import BSplineCurve2D, Circle2D, Cylinder, Line2D, Plane, Sphere

    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 1.0)

    # Exact: a circle in the plane's parameters, the iso-line v = 2 in the
    # cylinder's, both in the circle's parameter.
    (piece,) = Plane(Point(0.0, 0.0, 2.0), Vector.unit_z()).intersect(cylinder).curves
    assert isinstance(piece.first_uv_curve, Circle2D)
    assert isinstance(piece.second_uv_curve, Line2D)
    uv = piece.second_uv_curve.point_at(0.7)
    assert cylinder.point_at(uv.x, uv.y).distance(piece.curve.point_at(0.7)) < 1e-12
    assert piece.deviation == 0.0
    assert "deviation=0.0" in repr(piece)

    # Exact with no closed form on the cylinder: an ellipse.
    oblique = Plane(Point.origin(), Vector(0.0, 1.0, 1.0))
    (piece,) = oblique.intersect(cylinder).curves
    assert piece.second_uv_curve is None

    # Marched: fitted on both surfaces, in the curve's parameter.
    sphere = Sphere(Point.origin(), 2.0)
    around = Cylinder(Point(0.3, 0.0, 0.0), Vector.unit_z(), 1.0)
    for piece in sphere.intersect(around).curves:
        first, second = piece.first_uv_curve, piece.second_uv_curve
        assert isinstance(first, BSplineCurve2D) and isinstance(second, BSplineCurve2D)
        assert 0.0 < piece.deviation < 1e-6
        start, end = piece.curve.bounds()
        for k in range(11):
            t = start + (end - start) * k / 10
            p = piece.curve.point_at(t)
            u, v = first.point_at(t).x, first.point_at(t).y
            assert sphere.point_at(u, v).distance(p) <= 1.1 * piece.deviation
            u, v = second.point_at(t).x, second.point_at(t).y
            assert around.point_at(u, v).distance(p) <= 1.1 * piece.deviation


def test_intersect_torus_symmetric():
    from geomcore import Cone, Cylinder, IntersectionMethod, Plane, Torus

    torus = Torus(Point.origin(), Vector.unit_z(), 4.0, 1.0)

    c1, c2 = _curves(torus.intersect(Plane(Point.origin(), Vector.unit_z())))
    assert sorted([c1.radius(), c2.radius()]) == pytest.approx([3.0, 5.0])

    assert len(torus.intersect(Plane(Point.origin(), Vector.unit_x())).curves) == 2

    # Bitangent plane through the center: Villarceau circles.
    tilt = math.asin(1.0 / 4.0)
    villarceau = Plane(Point.origin(), Vector(0.0, math.sin(tilt), math.cos(tilt)))
    c1, c2 = _curves(torus.intersect(villarceau))
    assert [c1.radius(), c2.radius()] == pytest.approx([4.0, 4.0])
    assert sorted([c1.center().x, c2.center().x]) == pytest.approx([-1.0, 1.0])

    (piece,) = torus.intersect(Sphere(Point.origin(), 5.0)).curves
    assert piece.tangent is True
    assert piece.curve.radius() == pytest.approx(5.0)

    assert len(torus.intersect(Cylinder(Point.origin(), Vector.unit_z(), 4.5)).curves) == 2

    rings = torus.intersect(Cone.from_frame(Frame.world(), 1.1, 2.0)).curves
    assert [ring.tangent for ring in rings] == [False, False]
    rings = torus.intersect(Cone.from_frame(Frame.world(), 1.4, 0.0)).curves
    assert len(rings) == 4

    other = Torus(Point(0.0, 0.0, 1.5), Vector.unit_z(), 4.0, 1.0)
    assert len(torus.intersect(other).curves) == 2
    assert torus.intersect(torus).is_coincident()
    tilted = Torus(Point.origin(), Vector.unit_x(), 4.0, 1.0)
    assert torus.intersect(tilted).method == IntersectionMethod.NUMERIC


def test_conic_surface_intersections():
    from geomcore import Cone, Cylinder, Ellipse, Hyperbola, Parabola, Plane

    ellipse = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    assert len(ellipse.intersect(Plane(Point.origin(), Vector.unit_x())).hits) == 2

    parabola = Parabola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 1.0)
    on_plane = parabola.intersect(Plane(Point.origin(), Vector.unit_z()))
    assert [overlap.curve_range for overlap in on_plane.overlaps] == [(-math.inf, math.inf)]

    hyperbola = Hyperbola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 2.0, 1.0)
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    assert len(hyperbola.intersect(cylinder).hits) >= 1

    cone = Cone.from_frame(Frame.world(), 0.4, 2.0)
    assert ellipse.intersect(cone).hits


def test_circle_surface_and_measures():
    import math

    from geomcore import Cone, Cylinder, Ellipse, Plane

    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    assert circle.length() == pytest.approx(math.tau * 2.0)
    assert circle.enclosed_area() == pytest.approx(math.pi * 4.0)

    hits = circle.intersect(Plane(Point.origin(), Vector.unit_x())).hits
    assert len(hits) == 2

    assert circle.intersect(Sphere(Point.origin(), 3.0)).is_empty()

    on_cylinder = circle.intersect(Cylinder(Point.origin(), Vector.unit_z(), 2.0))
    assert [overlap.curve_range for overlap in on_cylinder.overlaps] == [pytest.approx((0.0, math.tau))]

    cone = Cone.from_frame(Frame.world(), 0.4, 2.0)
    assert [overlap.curve_range for overlap in circle.intersect(cone).overlaps] == [pytest.approx((0.0, math.tau))]


    ellipse = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    assert ellipse.enclosed_area() == pytest.approx(math.pi * 4.5)

    sphere = Sphere(Point.origin(), 3.0)
    assert sphere.area() == pytest.approx(4.0 * math.pi * 9.0)
    assert sphere.volume() == pytest.approx(4.0 / 3.0 * math.pi * 27.0)


def test_interpolate_curve():
    curve = BSplineCurve.interpolate(
        [Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 0.0), Point(2.0, 0.0, 0.0)], 2
    )
    assert curve.contains(Point(0.0, 0.0, 0.0))
    assert curve.contains(Point(2.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        BSplineCurve.interpolate([Point.origin()], 1)


def test_approximate_curve():
    import math

    pts = [Point(2.0 * math.cos(i / 12 * math.tau), 2.0 * math.sin(i / 12 * math.tau), 0.0) for i in range(12)]
    curve = BSplineCurve.approximate(pts, degree=3, pole_count=7)
    assert curve.degree() == 3
    with pytest.raises(ValueError):
        BSplineCurve.approximate(pts, 3, 20)


def test_curve_surface_intersect_any_pair():
    from geomcore import IntersectionMethod, Line

    line = Line(Point.origin(), Vector.unit_x())
    result = line.intersect(Sphere(Point.origin(), 2.0))
    assert [h.curve_parameter for h in result.hits] == pytest.approx([-2.0, 2.0])
    assert line.intersect(Sphere(Point(0.0, 0.0, 5.0), 2.0)).hits == []

    patch = _bspline_surface(0.0)
    vertical = Line(Point(0.5, 0.5, -1.0), Vector.unit_z())
    result = vertical.intersect(patch)
    assert result.method == IntersectionMethod.NUMERIC
    (hit,) = result.hits
    assert hit.point.z == pytest.approx(0.0, abs=1e-9)


def test_intersect_surfaces_march_without_closed_form():
    from geomcore import CompletenessStatus, IntersectionMethod, Plane

    patch = _bspline_surface(0.0)
    result = patch.intersect(Plane(Point(0.5, 0.0, 0.0), Vector.unit_x()))
    assert result.method == IntersectionMethod.NUMERIC
    assert result.completeness.status == CompletenessStatus.SAMPLED
    (piece,) = result.curves
    assert isinstance(piece.curve, BSplineCurve)

    miss = Plane(Point(0.0, 0.0, 5.0), Vector.unit_z())
    assert patch.intersect(miss).is_empty()


@pytest.mark.parametrize("height", [0.3, 0.8, -0.45])
def test_intersect_plane_sphere_off_grid(height: float) -> None:
    # Planes between the seed grid lines used to find no curve (#6).
    from geomcore import Plane

    plane = Plane(Point(0.0, 0.0, height), Vector.unit_z())
    (piece,) = plane.intersect(Sphere(Point.origin(), 1.0)).curves
    p = piece.curve.point_at(0.5)
    assert math.hypot(p.x, p.y) == pytest.approx(math.sqrt(1.0 - height**2), abs=1e-3)


def test_curve_curve_generic():
    from geomcore import Circle, Line

    l1 = Line(Point.origin(), Vector.unit_x())
    l2 = Line(Point(0.0, 0.0, 1.0), Vector.unit_y())
    extrema = l1.extrema(l2)
    assert extrema.first is l1 and extrema.second is l2
    assert len(extrema.minima) == 1
    assert extrema.minima[0].distance == pytest.approx(1.0)

    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    hits = l1.intersect(circle).hits
    assert [h.first_parameter for h in hits] == pytest.approx([-2.0, 2.0])


def test_curve_3d_to_uv_within():
    import math

    from geomcore import (
        BSplineCurve2D,
        Curve3DToUvApproximationError,
        CurveApproximation2D,
        Line,
        TrimmedCurve,
    )

    # A turn of a circle around a cylinder, starting at u = pi / 2, so it
    # crosses the seam.
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    turn = TrimmedCurve(circle, math.pi / 2, math.pi / 2 + math.tau)
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    approximation = cylinder.curve_3d_to_uv_within(turn, 1e-7)
    assert isinstance(approximation, CurveApproximation2D)
    assert approximation.deviation <= 1e-7
    assert not approximation.is_proven()
    uv_curve = approximation.curve
    assert isinstance(uv_curve, BSplineCurve2D)
    assert uv_curve.bounds() == turn.bounds()
    previous = uv_curve.point_at(math.pi / 2).x
    for k in range(1, 101):
        t = math.pi / 2 + math.tau * k / 100
        u, v = uv_curve.point_at(t).x, uv_curve.point_at(t).y
        assert previous < u < previous + 0.1
        assert v == pytest.approx(0.0, abs=1e-9)
        assert cylinder.point_at(u, v).distance(turn.point_at(t)) <= 1e-7
        previous = u
    assert repr(approximation).startswith("CurveApproximation2D(degree=5, pole_count=")

    with pytest.raises(Curve3DToUvApproximationError, match="unbounded"):
        cylinder.curve_3d_to_uv_within(Line(Point.origin(), Vector.unit_z()), 1e-7)
    with pytest.raises(Curve3DToUvApproximationError, match="from the surface"):
        wider = Circle(Point.origin(), Vector.unit_z(), 3.0)
        cylinder.curve_3d_to_uv_within(wider, 1e-7)

    # A quarter of a meridian, ending at the sphere's south pole: the UV
    # curve keeps the meridian's u up to the pole.
    sphere = Sphere(Point.origin(), 1.0)
    meridian = Circle(Point.origin(), Vector.unit_y(), 1.0)
    to_pole = TrimmedCurve(meridian, 0.0, math.pi / 2)
    uv_curve = sphere.curve_3d_to_uv_within(to_pole, 1e-6).curve
    start, end = uv_curve.point_at(0.0), uv_curve.point_at(math.pi / 2)
    assert end.x == pytest.approx(start.x, abs=1e-9)
    assert end.y == pytest.approx(-math.pi / 2, abs=1e-12)


def test_intersect_curves_2d_and_lines():
    from geomcore import Circle2D, Line, Line2D

    l1 = Line2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0))
    l2 = Line2D(Point2D(0.0, 0.0), Vector2D(0.0, 1.0))
    (hit,) = l1.intersect(l2).hits
    assert (hit.first_parameter, hit.second_parameter) == pytest.approx((0.0, 0.0))

    circle = Circle2D(Point2D(0.0, 0.0), 2.0)
    hits = l1.intersect(circle).hits
    assert [h.first_parameter for h in hits] == pytest.approx([-2.0, 2.0])

    p1, _ = [h.point for h in circle.intersect(Circle2D(Point2D(3.0, 0.0), 2.0)).hits]
    assert p1.x == pytest.approx(1.5)

    a = Line(Point.origin(), Vector.unit_x())
    b = Line(Point.origin(), Vector.unit_y())
    (hit,) = a.intersect(b).hits
    assert (hit.point.x, hit.point.y, hit.point.z) == pytest.approx((0.0, 0.0, 0.0))
    skew = a.intersect(Line(Point(0.0, 0.0, 1.0), Vector.unit_y()))
    assert skew.is_empty() and skew.is_complete()


def test_intersect_curves_3d():
    from geomcore import Circle, Line

    line = Line(Point.origin(), Vector.unit_x())
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    hits = line.intersect(circle).hits
    assert [h.first_parameter for h in hits] == pytest.approx([-2.0, 2.0])

    c1 = Circle(Point.origin(), Vector.unit_z(), 2.0)
    c2 = Circle(Point(3.0, 0.0, 0.0), Vector.unit_z(), 2.0)
    p1, _ = [h.point for h in c1.intersect(c2).hits]
    assert p1.x == pytest.approx(1.5)
    tilted = Circle(Point.origin(), Vector.unit_x(), 2.0)
    assert sorted(h.point.y for h in c1.intersect(tilted).hits) == pytest.approx([-2.0, 2.0])
    (overlap,) = c1.intersect(c1).overlaps
    assert overlap.first_range == pytest.approx((0.0, math.tau))
    lifted = Circle(Point(0.0, 0.0, 1.0), Vector.unit_z(), 2.0)
    assert c1.intersect(lifted).is_empty()


def _bilinear_patch():
    from geomcore import BSplineSurface

    return BSplineSurface(
        1,
        1,
        [
            [Point(0.0, 0.0, 0.0), Point(0.0, 2.0, 0.0)],
            [Point(2.0, 0.0, 0.0), Point(2.0, 2.0, 0.0)],
        ],
        [0.0, 1.0],
        [2, 2],
        [0.0, 1.0],
        [2, 2],
        False,
        False,
    )


def test_bspline_surface_projection():
    surface = _bilinear_patch()
    projection = surface.project_point(Point(1.0, 1.0, 1.0))
    assert (projection.u, projection.v) == pytest.approx((0.5, 0.5))
    assert projection.distance == pytest.approx(1.0)
    assert surface.contains(Point(1.0, 1.0, 0.0))
    assert not surface.contains(Point(1.0, 1.0, 1.0))
    ext = surface.extrema(Point(1.0, 1.0, 1.0))
    assert len(ext) >= 1
    assert ext[0].distance == pytest.approx(1.0)


def test_bspline_curve_projection():
    curve = BSplineCurve(
        1,
        [Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0)],
        [0.0, 1.0],
        [2, 2],
        False,
    )
    projection = curve.project_point(Point(1.0, 1.0, 0.0))
    assert projection.parameter == pytest.approx(0.5)
    assert projection.distance == pytest.approx(1.0)
    assert curve.contains(Point(1.0, 0.0, 0.0))
    assert not curve.contains(Point(1.0, 1.0, 0.0))
    # Beyond the end: the endpoint wins.
    projection = curve.project_point(Point(5.0, 0.0, 0.0))
    assert projection.parameter == pytest.approx(1.0)
    assert projection.distance == pytest.approx(3.0)
    ext = curve.extrema(Point(1.0, 1.0, 0.0))
    assert len(ext) >= 1
    assert ext[0].distance == pytest.approx(1.0)


def test_conic_extrema_and_projection():
    from geomcore import Ellipse, Hyperbola, Parabola

    ellipse = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    ext = ellipse.extrema(Point(4.0, 0.0, 0.0))
    assert len(ext) >= 2
    assert ext[0].distance == pytest.approx(1.0)
    projection = ellipse.project_point(Point(0.0, 3.0, 0.0))
    assert projection.distance == pytest.approx(1.5)

    parabola = Parabola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 1.0)
    projection = parabola.project_point(Point(1.0, 2.0, 0.0))
    assert projection.parameter == pytest.approx(2.0)
    assert projection.distance == 0.0

    hyperbola = Hyperbola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 2.0, 1.0)
    projection = hyperbola.project_point(Point(3.0, 0.0, 0.0))
    assert projection.distance == pytest.approx(math.sqrt(0.8))


def _hit_tuples(result):
    """Curve-curve hits as `(first_parameter, point, second_parameter,
    multiplicity)` tuples."""
    return [(h.first_parameter, h.point, h.second_parameter, h.multiplicity) for h in result.hits]


def test_line_conic_intersections():
    from geomcore import Ellipse, Hyperbola, Line, Parabola

    tol = Tolerance.default()
    ellipse = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    (t1, p1, u1, _), (t2, p2, u2, _) = _hit_tuples(
        Line(Point.origin(), Vector.unit_x()).intersect(ellipse)
    )
    assert (t1, t2) == pytest.approx((-3.0, 3.0))
    assert (u1, u2) == pytest.approx(
        (math.pi, 0.0)
    )

    parabola = Parabola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 1.0)
    # The axis crosses the parabola once, at its vertex.
    [(t, p, u, mult)] = _hit_tuples(Line(Point.origin(), Vector.unit_x()).intersect(parabola))
    assert mult == 1
    assert t == pytest.approx(0.0) and u == pytest.approx(0.0)

    hyperbola = Hyperbola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 2.0, 1.0)
    [(t, p, u, mult)] = _hit_tuples(
        Line(Point(2.0, -5.0, 0.0), Vector.unit_y()).intersect(hyperbola, tol)
    )
    assert mult == 2
    assert (p.x, p.y, p.z) == pytest.approx((2.0, 0.0, 0.0))

    miss = Line(Point(0.0, 0.0, 1.0), Vector.unit_x())
    assert miss.intersect(ellipse).is_empty()
    assert miss.intersect(parabola).is_empty()
    assert miss.intersect(hyperbola).is_empty()


def test_circle_conic_intersections():
    from geomcore import Circle, Ellipse, Hyperbola, Parabola

    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    ellipse = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    hits = _hit_tuples(circle.intersect(ellipse))
    assert len(hits) == 4
    for s, p, u, mult in hits:
        assert mult == 1
        assert circle.contains(p)
        assert ellipse.contains(p)
        _assert_same_point(circle.point_at(s), p)
        _assert_same_point(ellipse.point_at(u), p)

    same = circle.intersect(Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 2.0, 2.0))
    assert len(same.overlaps) == 1 and same.hits == []
    # Crossing planes: the tilted ellipse meets the shared line at
    # y = ±1.5, inside the circle.
    tilted = Ellipse(Point.origin(), Vector.unit_x(), Vector.unit_z(), 3.0, 1.5)
    assert circle.intersect(tilted).is_empty()

    unit = Circle(Point.origin(), Vector.unit_z(), 1.0)
    parabola = Parabola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 1.0)
    assert len(unit.intersect(parabola).hits) == 2

    big = Circle(Point.origin(), Vector.unit_z(), 3.0)
    hyperbola = Hyperbola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 2.0, 1.0)
    hits = big.intersect(hyperbola).hits
    assert len(hits) == 2
    assert all(h.point.x > 0.0 for h in hits)


def test_line_torus_intersections():
    from geomcore import Line, Torus

    torus = Torus(Point.origin(), Vector.unit_z(), 4.0, 1.0)
    hits = Line(Point.origin(), Vector.unit_x()).intersect(torus).hits
    assert [h.curve_parameter for h in hits] == pytest.approx([-5.0, -3.0, 3.0, 5.0])
    for hit in hits:
        assert hit.multiplicity == 1
        assert torus.contains(hit.point)
        projection = torus.project_point(hit.point)
        assert projection.distance == pytest.approx(0.0)
        assert (hit.u, hit.v) == pytest.approx((projection.u, projection.v))

    (hit,) = Line(Point(0.0, 5.0, 0.0), Vector.unit_x()).intersect(torus).hits
    assert hit.curve_parameter == pytest.approx(0.0)
    assert hit.multiplicity == 2

    assert Line(Point(0.0, 6.0, 0.0), Vector.unit_x()).intersect(torus).is_empty()


def test_ellipse_ellipse_intersections():
    from geomcore import Ellipse

    e1 = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    e2 = Ellipse(Point(0.5, 0.0, 0.0), Vector.unit_z(), Vector.unit_x(), 2.0, 2.0)
    hits = _hit_tuples(e1.intersect(e2))
    assert len(hits) == 4
    for s, p, t, _ in hits:
        assert e1.contains(p)
        assert e2.contains(p)

    grazing = Ellipse(Point(5.0, 0.0, 0.0), Vector.unit_z(), Vector.unit_x(), 2.0, 2.0)
    [(s, p, t, mult)] = _hit_tuples(e1.intersect(grazing))
    assert mult == 2
    assert (p.x, p.y, p.z) == pytest.approx((3.0, 0.0, 0.0))

    far = Ellipse(Point(10.0, 0.0, 0.0), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    assert e1.intersect(far).is_empty()
    assert len(e1.intersect(e1).overlaps) == 1
    # Crossing planes: both pass through (0, ±1.5, 0).
    tilted = Ellipse(Point.origin(), Vector.unit_x(), Vector.unit_z(), 3.0, 1.5)
    hits = _hit_tuples(e1.intersect(tilted))
    assert sorted(p.y for _, p, _, _ in hits) == pytest.approx([-1.5, 1.5])
    lifted = Ellipse(Point(0.0, 0.0, 1.0), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5)
    assert e1.intersect(lifted).is_empty()


def test_conic_conic_intersections():
    from geomcore import Ellipse, Hyperbola, Parabola

    o, z, x = Point.origin(), Vector.unit_z(), Vector.unit_x()
    ellipse = Ellipse(o, z, x, 3.0, 2.0)
    hyperbola = Hyperbola(o, z, x, 1.0, 0.5)
    parabola = Parabola(o, z, x, 1.0)
    # The branch x > 0 meets the ellipse at (√153/5, ±√128/10).
    hits = _hit_tuples(ellipse.intersect(hyperbola))
    hx, hy = math.sqrt(153) / 5, math.sqrt(128) / 10
    assert sorted(p.y for _, p, _, _ in hits) == pytest.approx([-hy, hy])
    assert [p.x for _, p, _, _ in hits] == pytest.approx([hx, hx])
    pairs = [
        (ellipse, hyperbola),
        (ellipse, parabola),
        (parabola, hyperbola),
        (parabola, Parabola(Point(2.0, 0.0, 0.0), z, Vector(-1.0, 0.0, 0.0), 1.0)),
        (hyperbola, Hyperbola(Point(0.0, -2.0, 0.0), z, Vector.unit_y(), 1.0, 1.0)),
    ]
    for a, b in pairs:
        hits = _hit_tuples(a.intersect(b))
        assert hits, (a, b)
        for s, p, t, _ in hits:
            # One conic is solved for; a parameter not computed is None.
            assert s is not None or t is not None
            _assert_same_point(a.point_at(s), p)
            _assert_same_point(b.point_at(t), p)
        # The receiver's parameter comes first.
        back = _hit_tuples(b.intersect(a))
        on_b = sorted(t for _, p, t, _ in hits)
        assert on_b == pytest.approx(sorted(t for t, p, _, _ in back))
        on_a = sorted(s for s, p, _, _ in hits)
        assert on_a == pytest.approx(sorted(s for _, p, s, _ in back))
    assert len(parabola.intersect(parabola).overlaps) == 1
    lifted = Parabola(Point(0.0, 0.0, 1.0), z, x, 1.0)
    assert parabola.intersect(lifted).is_empty()


def _rigid_motion():
    axis = Axis(Point(0.5, -1.0, 0.0), Vector(1.0, 2.0, 2.0))
    return Transform.translation(Vector(3.0, -1.0, 2.0)).then(Transform.rotation(axis, 0.7))


def _assert_same_point(a, b):
    assert (a.x, a.y, a.z) == pytest.approx((b.x, b.y, b.z), abs=1e-9)


def test_curves_transformed_commute_with_eval():
    from geomcore import Ellipse, Hyperbola, Line, Parabola

    t = _rigid_motion()
    curves = [
        Line(Point(1.0, -2.0, 0.5), Vector(1.0, 2.0, 2.0)),
        Circle(Point(1.0, 0.0, 0.0), Vector.unit_z(), 2.0),
        Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.5),
        Parabola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 0.5),
        Hyperbola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 2.0, 1.0),
        BSplineCurve(
            2,
            [Point(0.0, 0.0, 0.0), Point(1.0, 2.0, 0.0), Point(3.0, 1.0, 1.0)],
            [0.0, 1.0],
            [3, 3],
            False,
        ),
    ]
    for curve in curves:
        moved = curve.transformed(t)
        assert type(moved) is type(curve)
        for u in (0.0, 0.3, 0.9):
            _assert_same_point(moved.point_at(u), t.apply_point(curve.point_at(u)))


def test_surfaces_transformed_commute_with_eval():
    from geomcore import Cone, Plane, Torus

    t = _rigid_motion()
    surfaces = [
        Plane(Point(0.0, 0.0, 1.0), Vector(0.0, 1.0, 1.0)),
        Cylinder(Point.origin(), Vector.unit_z(), 2.0),
        Cone(Point.origin(), Vector.unit_z(), 0.4, 1.0),
        Sphere(Point(1.0, 2.0, 3.0), 2.0),
        Torus(Point.origin(), Vector.unit_z(), 3.0, 1.0),
        _bilinear_patch(),
    ]
    for surface in surfaces:
        moved = surface.transformed(t)
        assert type(moved) is type(surface)
        for u, v in ((0.0, 0.0), (0.3, 0.6), (0.9, 0.2)):
            _assert_same_point(moved.point_at(u, v), t.apply_point(surface.point_at(u, v)))


def test_transformed_scales_radius_and_rejects_collapse():
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    assert cylinder.transformed(Transform.scaling(Point.origin(), 3.0)).radius() == pytest.approx(6.0)
    with pytest.raises(ValueError):
        cylinder.transformed(Transform.scaling(Point.origin(), 0.0))


def test_eval_derivative_rejects_unsupported_orders():
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    with pytest.raises(ValueError):
        circle.derivative(0.0, 0)
    assert circle.derivative(0.0, 5).components()[1] == pytest.approx(2.0)
    curve = BSplineCurve.interpolate(
        [Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 0.0), Point(2.0, 0.0, 0.0)], 2
    )
    with pytest.raises(ValueError):
        curve.derivative(0.5, 0)
    assert curve.derivative(0.5, 3).components() == (0.0, 0.0, 0.0)
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 1.0)
    with pytest.raises(ValueError):
        cylinder.derivative(0.0, 0.0, 0, 0)
    third = cylinder.derivative(0.0, 0.0, 3, 0).components()
    assert third == pytest.approx((0.0, -1.0, 0.0))
    assert cylinder.derivative(0.0, 0.0, 2, 1).components() == (0.0, 0.0, 0.0)
    saddle = BezierSurface(1, 1, [[Point(0, 0, 0), Point(0, 1, 0)], [Point(1, 0, 0), Point(1, 1, 1)]])
    for du, dv in [(3, 0), (2, 1)]:
        with pytest.raises(ValueError, match="<= 2"):
            saddle.derivative(0.5, 0.5, du, dv)


def test_conic_projection_of_non_finite_point_is_nan():
    from geomcore import Ellipse, Hyperbola, Parabola

    o, z, x = Point.origin(), Vector.unit_z(), Vector.unit_x()
    for curve in [Ellipse(o, z, x, 3.0, 2.0), Hyperbola(o, z, x, 3.0, 2.0), Parabola(o, z, x, 1.0)]:
        for p in [Point(math.nan, 0.0, 0.0), Point(math.inf, 0.0, 0.0)]:
            projection = curve.project_point(p)
            assert math.isnan(projection.parameter) and math.isnan(projection.distance)
            assert curve.extrema(p) == []


def test_trimmed_curve_segment_and_arc():
    from geomcore import Line, TrimmedCurve

    segment = TrimmedCurve.segment(Point.origin(), Point(3.0, 4.0, 0.0))
    assert segment.bounds() == (0.0, 5.0)
    assert isinstance(segment.basis(), Line)
    assert segment.end_point().y == pytest.approx(4.0)
    # Beyond the end: the projection clamps to it.
    projection = segment.project_point(Point(6.0, 8.0, 0.0))
    assert projection.parameter == pytest.approx(5.0)
    assert projection.distance == pytest.approx(5.0)
    assert not segment.contains(Point(6.0, 8.0, 0.0))

    # The right half of the unit circle, across the seam at angle 0.
    circle = Circle(Point.origin(), Vector.unit_z(), 1.0)
    arc = TrimmedCurve(circle, 1.5 * math.pi, 2.5 * math.pi)
    assert isinstance(arc.basis(), Circle)
    projection = arc.project_point(Point(3.0, 0.0, 0.0))
    assert projection.parameter == pytest.approx(2.0 * math.pi)
    assert projection.distance == pytest.approx(2.0)
    projection = arc.project_point(Point(-3.0, 0.1, 0.0))
    assert projection.parameter == pytest.approx(2.5 * math.pi)
    assert arc.contains(Point(0.0, 1.0, 0.0))
    assert not arc.contains(Point(-1.0, 0.0, 0.0))
    assert arc.extrema(Point(3.0, 0.0, 0.0))[0].parameter == pytest.approx(2.0 * math.pi)


def test_trimmed_curves_name_their_range_by_its_ends():
    from geomcore import Line, TrimmedCurve

    line = Line(Point.origin(), Vector.unit_x())
    piece = TrimmedCurve(line, start_parameter=1.0, end_parameter=3.0)
    assert (piece.start_parameter(), piece.end_parameter()) == piece.bounds() == (1.0, 3.0)
    assert piece.point_at(piece.start_parameter()) == piece.start_point()
    assert repr(piece) == (
        "TrimmedCurve(start_parameter=1.0, end_parameter=3.0, "
        "start_point=(1.0, 0.0, 0.0), end_point=(3.0, 0.0, 0.0))"
    )


def test_trimmed_curve_rejects_bad_ranges():
    from geomcore import Line, TrimmedCurve

    circle = Circle(Point.origin(), Vector.unit_z(), 1.0)
    for first, last in [(1.0, 1.0), (2.0, 1.0), (0.0, 7.0), (math.nan, 1.0)]:
        with pytest.raises(ValueError):
            TrimmedCurve(circle, first, last)
    with pytest.raises(ValueError):
        TrimmedCurve.segment(Point.origin(), Point.origin())
    with pytest.raises(TypeError):
        TrimmedCurve(TrimmedCurve(Line(Point.origin(), Vector.unit_x()), 0.0, 1.0), 0.0, 1.0)
    curve = BSplineCurve(1, [Point.origin(), Point(2.0, 0.0, 0.0)], [0.0, 1.0], [2, 2], False)
    with pytest.raises(ValueError):
        TrimmedCurve(curve, -0.5, 0.5)
    with pytest.raises(ValueError):
        TrimmedCurve(curve, 0.0, 1.0).derivative(0.5, 0)


def test_intersect_searches_only_the_range():
    from geomcore import TrimmedCurve

    segment = TrimmedCurve.segment(Point(-2.0, 0.5, 0.0), Point(2.0, 0.5, 0.0))
    circle = Circle(Point.origin(), Vector.unit_z(), 1.0)
    assert len(segment.intersect(circle).hits) == 2
    quarter = TrimmedCurve(circle, 0.0, math.pi / 2)
    (hit,) = segment.intersect(quarter).hits
    assert hit.point.x > 0.0

    a = TrimmedCurve.segment(Point.origin(), Point(1.0, 0.0, 0.0))
    b = TrimmedCurve.segment(Point(3.0, -1.0, 1.0), Point(3.0, 1.0, 1.0))
    assert a.extrema(b).minima[0].distance == pytest.approx(math.sqrt(5.0))

    sphere = Sphere(Point.origin(), 1.0)
    half = TrimmedCurve.segment(Point.origin(), Point(2.0, 0.0, 0.0))
    (hit,) = half.intersect(sphere).hits
    assert hit.curve_parameter == pytest.approx(1.0)


def test_trimmed_curve_2d():
    from geomcore import Circle2D, TrimmedCurve2D

    arc = TrimmedCurve2D(Circle2D(Point2D.origin(), 2.0), 0.0, math.pi / 2)
    assert isinstance(arc.basis(), Circle2D)
    projection = arc.project_point(Point2D(1.0, 1.0))
    assert projection.parameter == pytest.approx(math.pi / 4)
    projection = arc.project_point(Point2D(-1.0, -5.0))
    assert projection.parameter == 0.0
    assert arc.contains(Point2D(0.0, 2.0))
    assert not arc.contains(Point2D(-2.0, 0.0))

    segment = TrimmedCurve2D.segment(Point2D.origin(), Point2D(0.0, 2.0))
    assert isinstance(segment.basis(), Line2D)
    assert segment.bounds() == (0.0, 2.0)
    with pytest.raises(TypeError):
        TrimmedCurve2D(Circle(Point.origin(), Vector.unit_z(), 1.0), 0.0, 1.0)


def test_trimmed_curve_arcs():
    from geomcore import ConstructionError, TrimmedCurve

    # Three quarters of the unit circle, clockwise seen from +z.
    start, middle, end = Point(1.0, 0.0, 0.0), Point(-1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)
    arc = TrimmedCurve.arc_through(start, middle, end)
    assert arc.sweep_angle() == pytest.approx(1.5 * math.pi)
    assert arc.basis().normal().z == pytest.approx(-1.0)
    assert arc.end_point().distance(end) == pytest.approx(0.0, abs=1e-12)
    assert arc.length(1e-9) == pytest.approx(1.5 * math.pi)

    quarter = TrimmedCurve.arc(Point.origin(), Point(2.0, 0.0, 0.0), Point(0.0, 2.0, 0.0), Vector.unit_z())
    assert quarter.sweep_angle() == pytest.approx(math.pi / 2)
    rest = TrimmedCurve.arc(Point.origin(), Point(2.0, 0.0, 0.0), Point(0.0, 2.0, 0.0), -Vector.unit_z())
    assert rest.sweep_angle() == pytest.approx(1.5 * math.pi)
    full = TrimmedCurve.arc(Point.origin(), start, start, Vector.unit_z())
    assert full.sweep_angle() == pytest.approx(math.tau)

    with pytest.raises(ConstructionError):
        TrimmedCurve.arc_through(start, start, end)
    with pytest.raises(ConstructionError):
        TrimmedCurve.arc(Point.origin(), start, Point(0.0, 1.5, 0.0), Vector.unit_z())
    with pytest.raises(TypeError):
        TrimmedCurve.segment(Point.origin(), start).sweep_angle()


def test_trimmed_curve_2d_arcs():
    from geomcore import Circle2D, TrimmedCurve2D

    start, middle, end = Point2D(1.0, 0.0), Point2D(0.0, 1.0), Point2D(-1.0, 0.0)
    left = TrimmedCurve2D.arc_through(start, middle, end)
    right = TrimmedCurve2D.arc_through(end, middle, start)
    assert left.basis().is_counterclockwise()
    assert not right.basis().is_counterclockwise()
    assert left.sweep_angle() == pytest.approx(math.pi)

    clockwise = TrimmedCurve2D.arc(Point2D.origin(), Point2D(2.0, 0.0), Point2D(0.0, 2.0), False)
    assert clockwise.sweep_angle() == pytest.approx(1.5 * math.pi)
    assert clockwise.length(1e-9) == pytest.approx(3.0 * math.pi)
    assert Circle2D(Point2D.origin(), 2.0).length() == pytest.approx(4.0 * math.pi)


def test_lengths_of_bounded_curves():
    from geomcore import TrimmedCurve

    assert TrimmedCurve.segment(Point.origin(), Point(1.0, 2.0, 2.0)).length(1e-9) == 3.0
    ellipse = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 2.0, 1.0)
    assert ellipse.length(1e-12) == pytest.approx(9.688448220547675, abs=1e-11)
    quarter = TrimmedCurve(ellipse, 0.0, math.pi / 2)
    assert quarter.length(1e-12) == pytest.approx(ellipse.length(1e-12) / 4.0, abs=1e-11)
    ellipse_2d = Ellipse2D(Point2D.origin(), Vector2D.unit_x(), 2.0, 1.0)
    assert ellipse_2d.length(1e-12) == pytest.approx(ellipse.length(1e-12), abs=1e-11)


def test_curve_3d_to_uv_maps_a_trimmed_curve_to_a_trimmed_uv_curve():
    from geomcore import Line2D, TrimmedCurve, TrimmedCurve2D

    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    circle = Circle(Point(0.0, 0.0, 3.0), Vector.unit_z(), 2.0)
    edge = TrimmedCurve(circle, 0.5, 2.0)
    uv_edge = cylinder.curve_3d_to_uv(edge)
    assert isinstance(uv_edge, TrimmedCurve2D)
    assert isinstance(uv_edge.basis(), Line2D)
    assert (uv_edge.start_parameter(), uv_edge.end_parameter()) == (0.5, 2.0)
    end = uv_edge.end_point()
    assert (end.x, end.y) == (pytest.approx(2.0), pytest.approx(3.0))


def _intersect_shapes():
    from geomcore import Circle2D, Cone, Ellipse, Hyperbola, Line, Parabola, Plane, Torus

    o, z, x = Point.origin(), Vector.unit_z(), Vector.unit_x()
    return [
        Plane(Point(0.0, 0.0, 0.5), z),
        Cylinder(o, z, 2.0),
        Cone(o, z, 0.4, 2.0),
        Sphere(o, 3.0),
        Torus(o, z, 4.0, 1.0),
        Line(Point(0.0, -0.5, 0.25), x),
        Circle(o, z, 2.5),
        Ellipse(o, z, x, 3.0, 1.5),
        Parabola(o, z, x, 1.0),
        Hyperbola(o, z, x, 2.0, 1.0),
        Line2D(Point2D(0.0, 0.5), Vector2D(1.0, 0.0)),
        Circle2D(Point2D.origin(), 1.0),
    ]


def _is_surface(shape):
    return isinstance(shape, get_args(geomcore.Surface))


def _rounded(p):
    return tuple(round(c, 9) for c in ((p.x, p.y, p.z) if hasattr(p, "z") else (p.x, p.y)))


def _hit_points(result):
    """The points of a result's hits, rounded and sorted."""
    return sorted(_rounded(h.point) for h in result.hits)


def test_intersect_gives_the_same_result_in_both_orders():
    shapes = _intersect_shapes()
    pairs = 0
    for a in shapes:
        for b in shapes:
            if ("2D" in type(a).__name__) != ("2D" in type(b).__name__):
                continue
            forward, backward = a.intersect(b), b.intersect(a)
            assert forward.first is a and forward.second is b
            assert backward.first is b and backward.second is a
            assert type(forward) is type(backward)
            assert forward.method == backward.method, (a, b)
            if _is_surface(a) and _is_surface(b):
                assert len(forward.curves) == len(backward.curves), (a, b)
                assert forward.is_coincident() == backward.is_coincident()
            elif _is_surface(a) != _is_surface(b):
                # The curve's parameter keeps its role in either order.
                assert forward.hits == backward.hits, (a, b)
                assert forward.overlaps == backward.overlaps, (a, b)
            else:
                assert _hit_points(forward) == _hit_points(backward), (a, b)
                assert len(forward.overlaps) == len(backward.overlaps), (a, b)
            pairs += 1
    assert pairs == 104


def test_intersect_rejects_mixed_dimensions_and_foreign_objects():
    from geomcore import Circle2D

    o, z = Point.origin(), Vector.unit_z()
    circle = Circle(o, z, 2.5)
    for a, b in [
        (circle, Circle2D(Point2D.origin(), 1.0)),
        (Sphere(o, 1.0), Circle2D(Point2D.origin(), 1.0)),
        (circle, 5),
    ]:
        with pytest.raises(TypeError):
            a.intersect(b)
        if hasattr(b, "intersect"):
            with pytest.raises(TypeError):
                b.intersect(a)


def _xyz(p):
    return (p.x, p.y, p.z)


def test_discretize_circle_and_bspline():
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    # 2(1 - cos(pi/8)) = 0.152...: eight segments meet 0.16.
    d = circle.discretize(0.16)
    assert len(d.points) == len(d.parameters) == 9
    assert len(d.chord_deviations) == len(d.angle_deviations) == 8
    assert d.parameters[0] == 0.0 and d.parameters[-1] == pytest.approx(2.0 * math.pi)
    assert _xyz(d.points[0]) == _xyz(d.points[-1])
    assert all(dev == pytest.approx(2.0 * (1.0 - math.cos(math.pi / 8))) for dev in d.chord_deviations)
    assert all(angle == pytest.approx(math.pi / 8) for angle in d.angle_deviations)
    # A tighter angle limit asks for more points.
    assert len(circle.discretize(0.16, max_angle_deviation=0.1).points) == 33

    curve = BSplineCurve(
        2,
        [Point(0.0, 0.0, 0.0), Point(1.0, 2.0, 0.0), Point(2.0, 0.0, 0.0)],
        [0.0, 1.0],
        [3, 3],
        False,
    )
    d = curve.discretize(1e-3)
    assert _xyz(d.points[0]) == (0.0, 0.0, 0.0)
    assert _xyz(d.points[-1]) == (2.0, 0.0, 0.0)
    assert all(dev <= 1e-3 for dev in d.chord_deviations)
    for u, p in zip(d.parameters, d.points):
        assert curve.point_at(u).distance(p) < 1e-12

    with pytest.raises(ValueError, match="chord deviation"):
        circle.discretize(0.0)
    with pytest.raises(ValueError, match="angle deviation"):
        curve.discretize(1e-3, max_angle_deviation=-1.0)


def test_discretize_ellipse_and_trimmed_curves():
    from geomcore import Ellipse, Parabola, TrimmedCurve

    # A round ellipse is a circle: 2(1 - cos(pi/8)) = 0.152... gives eight segments.
    round_ellipse = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 2.0, 2.0)
    d = round_ellipse.discretize(0.16)
    assert len(d.points) == 9
    assert _xyz(d.points[0]) == _xyz(d.points[-1])
    flat = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 5.0, 0.5)
    d = flat.discretize(1e-3, max_angle_deviation=0.1)
    assert all(dev <= 1e-3 for dev in d.chord_deviations)
    assert all(angle <= 0.1 for angle in d.angle_deviations)
    # A segment is its own chord.
    segment = TrimmedCurve.segment(Point(0.0, 0.0, 0.0), Point(1.0, 2.0, 2.0))
    d = segment.discretize(1e-7, max_angle_deviation=1e-9)
    assert [_xyz(p) for p in d.points] == [(0.0, 0.0, 0.0), (1.0, 2.0, 2.0)]
    assert d.chord_deviations == [0.0]
    arc = TrimmedCurve(Parabola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 0.5), -2.0, 3.0)
    d = arc.discretize(1e-3)
    assert d.parameters[0] == -2.0 and d.parameters[-1] == 3.0
    assert all(dev <= 1e-3 for dev in d.chord_deviations)


def test_discretize_2d_curves():
    from geomcore import Circle2D, TrimmedCurve2D

    circle = Circle2D(Point2D(0.5, 0.25), 2.0)
    # 4 sin^2(pi/16) = 0.152...: eight segments meet 0.16.
    d = circle.discretize(0.16)
    assert isinstance(d, geomcore.Discretization2D)
    assert len(d.points) == 9 and len(d.chord_deviations) == 8
    assert (d.points[0].x, d.points[0].y) == (d.points[-1].x, d.points[-1].y)
    assert all(angle == pytest.approx(math.pi / 8) for angle in d.angle_deviations)
    arc = TrimmedCurve2D(circle, 0.5, 2.0)
    d = arc.discretize(1e-3, max_angle_deviation=0.05)
    assert d.parameters[0] == 0.5 and d.parameters[-1] == 2.0
    assert all(dev <= 1e-3 for dev in d.chord_deviations)
    assert all(angle <= 0.05 for angle in d.angle_deviations)
    segment = TrimmedCurve2D.segment(Point2D(0.0, 0.0), Point2D(3.0, 4.0))
    d = segment.discretize(1e-7)
    assert [(p.x, p.y) for p in d.points] == [(0.0, 0.0), (3.0, 4.0)]
    assert d.chord_deviations == [0.0] and d.angle_deviations == [0.0]
    with pytest.raises(ValueError):
        segment.discretize(0.0)
    # Any positive deviation is accepted, however small; a segment stays two points.
    assert len(segment.discretize(1e-9).points) == 2
    arc = geomcore.TrimmedCurve2D(geomcore.Circle2D(Point2D(0.0, 0.0), 1.0), 0.0, 1.0)
    fine = arc.discretize(1e-9)
    assert all(dev <= 1e-9 for dev in fine.chord_deviations)
    assert len(fine.points) > len(arc.discretize(1e-3).points)


def _cubic_curve():
    return BSplineCurve(
        3,
        [Point(0.0, 0.0, 0.0), Point(1.0, 2.0, 0.0), Point(3.0, 2.0, 1.0), Point(4.0, 0.0, 1.0)],
        [0.0, 1.0],
        [4, 4],
        False,
    )


def test_bspline_curve_getters():
    curve = _cubic_curve()
    assert [ (p.x, p.y, p.z) for p in curve.poles() ] == [
        (0.0, 0.0, 0.0), (1.0, 2.0, 0.0), (3.0, 2.0, 1.0), (4.0, 0.0, 1.0)]
    assert curve.weights() is None
    assert curve.knots() == [0.0, 1.0]
    assert curve.multiplicities() == [4, 4]
    rational = BSplineCurve.new_rational(
        2,
        [Point(1.0, 0.0, 0.0), Point(1.0, 1.0, 0.0), Point(0.0, 1.0, 0.0)],
        [1.0, math.sqrt(2) / 2, 1.0],
        [0.0, 1.0],
        [3, 3],
        False,
    )
    assert rational.is_rational()
    assert rational.weights() == pytest.approx([1.0, math.sqrt(2) / 2, 1.0])


def test_bspline_curve_editing_round_trip():
    curve = _cubic_curve()
    refined = curve.refine_knots([0.25, 0.5, 0.5, 0.75])
    assert len(refined.poles()) == 8
    for u in (0.0, 0.3, 0.7, 1.0):
        assert refined.point_at(u).distance(curve.point_at(u)) < 1e-12
    simpler, removed = refined.remove_knots([0.25, 0.5, 0.5, 0.75], 1e-9)
    assert removed == 4
    assert simpler.knots() == curve.knots()
    single = curve.insert_knot(0.5, 2)
    assert len(single.poles()) == 6
    back, removed = single.remove_knot(0.5, 2, 1e-9)
    assert removed == 2
    assert back.knots() == curve.knots()
    raised = curve.elevate_degree(1)
    assert raised.degree() == 4
    assert raised.point_at(0.4).distance(curve.point_at(0.4)) < 1e-12
    with pytest.raises(ValueError):
        curve.insert_knot(2.0, 1)
    with pytest.raises(ValueError):
        curve.remove_knots([0.5], 0.0)


def test_bspline_curve_split_extract_reverse():
    curve = _cubic_curve()
    first, second = curve.split_at(0.4)
    assert first.bounds() == pytest.approx((0.0, 0.4))
    assert second.bounds() == pytest.approx((0.4, 1.0))
    assert first.point_at(0.4).distance(curve.point_at(0.4)) < 1e-12
    assert second.point_at(0.4).distance(curve.point_at(0.4)) < 1e-12
    piece = curve.extract(0.2, 0.8)
    assert piece.bounds() == pytest.approx((0.2, 0.8))
    assert piece.point_at(0.5).distance(curve.point_at(0.5)) < 1e-12
    rev = curve.reversed()
    assert rev.point_at(0.3).distance(curve.point_at(0.7)) < 1e-12
    with pytest.raises(ValueError):
        curve.split_at(0.0)
    with pytest.raises(ValueError):
        curve.extract(0.8, 0.2)


def test_bspline_curve_bezier_pieces_and_length():
    curve = BSplineCurve(
        3,
        [Point(float(i), float(i % 2), 0.0) for i in range(6)],
        [0.0, 1.0, 2.0, 3.0],
        [4, 1, 1, 4],
        False,
    )
    pieces = curve.to_bezier_pieces()
    assert len(pieces) == 3
    assert all(isinstance(piece, BezierCurve) for piece in pieces)
    assert pieces[1].bounds() == pytest.approx((1.0, 2.0))
    assert len(pieces[1].poles()) == 4
    assert pieces[1].point_at(1.5).distance(curve.point_at(1.5)) < 1e-12
    assert pieces[1].to_bspline().multiplicities() == [4, 4]
    segment = BSplineCurve(
        1, [Point.origin(), Point(2.0, 0.0, 0.0)], [0.0, 1.0], [2, 2], False
    )
    assert segment.length(max_length_error=1e-9) == pytest.approx(2.0)
    assert _cubic_curve().length(1e-9) == pytest.approx(5.39065, abs=1e-4)


def test_bspline_curve_approximate_within():
    points = [Point(float(i), 0.1 * float(i), 0.0) for i in range(6)]
    fit = BSplineCurve.approximate_within(points, 2, 1e-6)
    assert isinstance(fit, CurveApproximation)
    assert fit.deviation <= 1e-6 and fit.is_proven()
    for p in points:
        assert fit.curve.contains(p, Tolerance(1e-3, 1e-9, 1e-9))
    assert repr(fit).startswith("CurveApproximation(degree=2, pole_count=")


def test_conic_to_bspline():
    from geomcore import Ellipse, TrimmedCurve

    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    spline = circle.to_bspline()
    assert spline.is_periodic()
    assert spline.point_at(math.pi / 2).distance(circle.point_at(math.pi / 2)) < 1e-12
    ellipse = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.0)
    e_spline = ellipse.to_bspline()
    assert e_spline.point_at(0.0).distance(ellipse.point_at(0.0)) < 1e-12
    arc = TrimmedCurve(circle, 0.0, math.pi)
    arc_spline = arc.to_bspline()
    assert arc_spline.bounds() == pytest.approx((0.0, math.pi))
    # The B-spline parameter matches the angle at the knots and each
    # piece's middle, not in between — but every point lies on the circle.
    assert arc_spline.point_at(math.pi / 2).distance(circle.point_at(math.pi / 2)) < 1e-12
    for u in (0.0, 0.3, 1.0, 2.0, math.pi):
        assert circle.contains(arc_spline.point_at(u))
    segment = TrimmedCurve.segment(Point.origin(), Point(1.0, 0.0, 0.0))
    seg_spline = segment.to_bspline()
    assert seg_spline.point_at(0.5).distance(Point(0.5, 0.0, 0.0)) < 1e-12


def test_curve_bounds_and_boxes():
    from geomcore import Ellipse, Hyperbola, Line, Parabola, TrimmedCurve

    line = Line(Point.origin(), Vector.unit_x())
    assert line.bounds() == (float("-inf"), float("inf"))
    assert line.bounding_box() is None
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    assert circle.bounds() == pytest.approx((0.0, math.tau))
    box = circle.bounding_box()
    assert (box.min_corner().x, box.max_corner().x) == pytest.approx((-2.0, 2.0))
    assert (box.min_corner().y, box.max_corner().y) == pytest.approx((-2.0, 2.0))
    assert box.encloses(Point(1.0, 1.0, 0.0))
    assert not box.encloses(Point(3.0, 0.0, 0.0))
    assert Parabola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 1.0).bounding_box() is None
    assert Hyperbola(
        Point.origin(), Vector.unit_z(), Vector.unit_x(), 1.0, 1.0
    ).bounding_box() is None
    assert Ellipse(
        Point.origin(), Vector.unit_z(), Vector.unit_x(), 2.0, 1.0
    ).bounding_box() is not None
    curve = _cubic_curve()
    cbox = curve.bounding_box()
    assert (cbox.min_corner().x, cbox.max_corner().x) == pytest.approx((0.0, 4.0))
    trimmed = TrimmedCurve(circle, 0.0, math.pi / 2)
    tbox = trimmed.bounding_box()
    assert tbox is not None and tbox.encloses(Point(2.0, 0.0, 0.0))


def test_bounding_box_ops():
    box = BoundingBox(Point(1.0, 0.0, 5.0), Point(0.0, 2.0, 3.0))
    assert (box.min_corner().x, box.min_corner().y, box.min_corner().z) == (0.0, 0.0, 3.0)
    assert (box.max_corner().x, box.max_corner().y, box.max_corner().z) == (1.0, 2.0, 5.0)
    assert (box.center().x, box.center().y, box.center().z) == (0.5, 1.0, 4.0)
    assert box.diagonal().components() == (1.0, 2.0, 2.0)
    other = BoundingBox(Point(2.0, -1.0, 0.0), Point(3.0, 0.0, 0.5))
    both = box.union(other)
    assert (both.min_corner().x, both.max_corner().x) == (0.0, 3.0)
    grown = box.enlarged(0.5)
    assert grown.min_corner().x == pytest.approx(-0.5)
    assert box.overlaps(other) or not box.overlaps(other)
    assert not box.overlaps(BoundingBox(Point(5.0, 5.0, 5.0), Point(6.0, 6.0, 6.0)))
    assert BoundingBox.from_points([Point(1.0, 2.0, 3.0)]).min_corner().x == 1.0
    with pytest.raises(ValueError):
        BoundingBox.from_points([])


def test_surface_getters_and_uv_bounds():
    surface = _bilinear_patch()
    assert (surface.u_degree(), surface.v_degree()) == (1, 1)
    assert surface.u_pole_count() == 2 and surface.v_pole_count() == 2
    assert [(p.x, p.y) for p in surface.poles()[1]] == [(2.0, 0.0), (2.0, 2.0)]
    assert surface.weights() is None
    assert surface.u_knots() == [0.0, 1.0] and surface.v_knots() == [0.0, 1.0]
    assert surface.u_multiplicities() == [2, 2]
    assert surface.u_bounds() == (0.0, 1.0) and surface.v_bounds() == (0.0, 1.0)
    sphere = Sphere(Point.origin(), 2.0)
    assert sphere.u_bounds() == pytest.approx((0.0, math.tau))
    assert sphere.bounding_box().max_corner().x == pytest.approx(2.0)
    assert sphere.bounding_box().min_corner().x == pytest.approx(-2.0)
    from geomcore import Plane

    plane = Plane(Point.origin(), Vector.unit_z())
    assert plane.bounding_box() is None
    assert plane.u_bounds() == (float("-inf"), float("inf"))


def test_surface_editing_round_trip():
    surface = _bilinear_patch()
    refined = surface.refine_u_knots([0.25, 0.5])
    assert refined.u_pole_count() == 4
    assert refined.point_at(0.3, 0.7).distance(surface.point_at(0.3, 0.7)) < 1e-12
    simpler, removed = refined.remove_u_knots([0.25, 0.5], 1e-9)
    assert removed == 2
    assert simpler.u_knots() == surface.u_knots()
    one, n = surface.insert_v_knot(0.5, 1).remove_v_knot(0.5, 1, 1e-9)
    assert n == 1 and one.v_knots() == surface.v_knots()
    raised = surface.elevate_u_degree(1)
    assert raised.u_degree() == 2
    assert raised.point_at(0.4, 0.6).distance(surface.point_at(0.4, 0.6)) < 1e-12
    raised_v = surface.elevate_v_degree(2)
    assert raised_v.v_degree() == 3
    with pytest.raises(ValueError):
        surface.insert_u_knot(2.0, 1)


def test_surface_split_extract_reverse_and_patches():
    surface = _bilinear_patch()
    left, right = surface.split_u_at(0.5)
    assert left.u_bounds() == pytest.approx((0.0, 0.5))
    assert right.point_at(0.75, 0.25).distance(surface.point_at(0.75, 0.25)) < 1e-12
    bottom, _top = surface.split_v_at(0.5)
    assert bottom.v_bounds() == pytest.approx((0.0, 0.5))
    piece = surface.extract((0.25, 0.75), (0.25, 0.75))
    assert piece.point_at(0.5, 0.5).distance(surface.point_at(0.5, 0.5)) < 1e-12
    assert surface.extract_u(0.25, 0.75).u_knots() == pytest.approx([0.25, 0.75])
    assert surface.extract_v(0.25, 0.75).v_knots() == pytest.approx([0.25, 0.75])
    assert surface.reversed_u().point_at(0.25, 0.5).distance(
        surface.point_at(0.75, 0.5)) < 1e-12
    assert surface.reversed_v().point_at(0.5, 0.25).distance(
        surface.point_at(0.5, 0.75)) < 1e-12
    from geomcore import BSplineSurface

    net = BSplineSurface(
        2, 2,
        [[Point(float(i), float(j), float(i * j)) for j in range(3)] for i in range(4)],
        [0.0, 1.0, 2.0], [3, 1, 3], [0.0, 1.0], [3, 3], False, False,
    )
    patches = net.to_bezier_patches()
    assert (len(patches), len(patches[0])) == (2, 1)
    assert isinstance(patches[1][0], BezierSurface)
    assert patches[1][0].u_bounds() == pytest.approx((1.0, 2.0))
    assert patches[1][0].point_at(1.5, 0.5).distance(net.point_at(1.5, 0.5)) < 1e-12
    assert patches[1][0].to_bspline().u_multiplicities() == [3, 3]
    with pytest.raises(ValueError):
        surface.split_u_at(0.0)


def test_surface_approximate():
    from geomcore import BSplineSurface

    grid = [
        [Point(x, y, math.sin(3 * x) * math.cos(2 * y)) for j in range(8) for y in [j / 7]]
        for i in range(8) for x in [i / 7]
    ]
    fit = BSplineSurface.approximate(grid, 3, 3, 5, 5)
    assert (fit.u_pole_count(), fit.v_pole_count()) == (5, 5)
    loose_fit = BSplineSurface.approximate_within(grid, 3, 3, 0.05)
    assert isinstance(loose_fit, SurfaceApproximation)
    assert loose_fit.deviation <= 0.05 and loose_fit.is_proven()
    assert repr(loose_fit).startswith("SurfaceApproximation(u_degree=3, v_degree=3, pole_count=")
    loose = loose_fit.surface
    assert loose.u_pole_count() <= 8 and loose.v_pole_count() <= 8
    assert loose.point_at(0.5, 0.5).distance(fit.point_at(0.5, 0.5)) < 0.2


def test_surface_to_bspline_conversions():
    from geomcore import BSplineSurface, Cone, Cylinder, Plane, Torus

    plane = Plane(Point.origin(), Vector.unit_z())
    patch = plane.to_bspline_patch((0.0, 1.0), (0.0, 2.0))
    assert isinstance(patch, BSplineSurface)
    assert patch.point_at(0.5, 1.0).distance(plane.point_at(0.5, 1.0)) < 1e-12
    cylinder = Cylinder(Point.origin(), Vector.unit_z(), 2.0)
    tube = cylinder.to_bspline_patch((0.0, math.tau), (-1.0, 4.0))
    assert tube.is_u_periodic() and not tube.is_v_periodic()
    assert tube.point_at(math.pi / 2, 3.0).distance(
        cylinder.point_at(math.pi / 2, 3.0)) < 1e-12
    sphere = Sphere(Point.origin(), 2.0)
    ball = sphere.to_bspline()
    assert ball.point_at(0.0, 0.0).distance(sphere.point_at(0.0, 0.0)) < 1e-9
    torus = Torus(Point.origin(), Vector.unit_z(), 3.0, 1.0)
    donut = torus.to_bspline()
    assert donut.point_at(0.0, 0.0).distance(torus.point_at(0.0, 0.0)) < 1e-9
    cone = Cone(Point.origin(), Vector.unit_z(), 0.5, 1.0)
    assert cone.to_bspline_patch((0.0, 1.0), (0.0, 2.0)).point_at(
        0.5, 1.0).distance(cone.point_at(0.5, 1.0)) < 1e-9
    with pytest.raises(ValueError):
        plane.to_bspline_patch((0.0, 0.0), (0.0, 1.0))


def _iso_surfaces():
    from geomcore import Cone, Plane, Torus

    return [
        Plane(Point(1.0, 2.0, 3.0), Vector(0.0, 0.0, 1.0)),
        Cylinder(Point.origin(), Vector.unit_z(), 2.0),
        Cone(Point.origin(), Vector.unit_z(), 0.4, 1.0),
        Sphere(Point(1.0, 0.0, 0.0), 3.0),
        Torus(Point.origin(), Vector.unit_z(), 5.0, 1.5),
        Sphere(Point.origin(), 2.0).to_bspline_patch((0.0, 3.0), (-1.0, 1.0)),
    ]


@pytest.mark.parametrize("index", range(6))
def test_iso_curves_are_the_surface_with_one_parameter_fixed(index):
    s = _iso_surfaces()[index]

    def finite(lo, hi):
        return (lo, hi) if math.isfinite(lo) and math.isfinite(hi) else (-2.0, 2.0)

    u0, u1 = finite(*s.u_bounds())
    v0, v1 = finite(*s.v_bounds())
    for k in (0.2, 0.5, 0.9):
        u, v = u0 + k * (u1 - u0), v0 + k * (v1 - v0)
        for t in (0.1, 0.6):
            tv, tu = v0 + t * (v1 - v0), u0 + t * (u1 - u0)
            for got, want in (
                (s.u_iso(u).point_at(tv), s.point_at(u, tv)),
                (s.v_iso(v).point_at(tu), s.point_at(tu, v)),
            ):
                assert (got.x, got.y, got.z) == pytest.approx((want.x, want.y, want.z), abs=1e-9)


def test_iso_curves_have_the_kernel_types():
    from geomcore import Line

    plane, cylinder, cone, sphere, torus, patch = _iso_surfaces()
    assert isinstance(plane.u_iso(0.0), Line) and isinstance(plane.v_iso(0.0), Line)
    assert isinstance(cylinder.u_iso(0.0), Line) and isinstance(cylinder.v_iso(0.0), Circle)
    assert isinstance(cone.u_iso(0.0), Line) and isinstance(cone.v_iso(0.0), Circle)
    assert isinstance(sphere.u_iso(0.0), Circle) and isinstance(sphere.v_iso(0.0), Circle)
    assert isinstance(torus.u_iso(0.0), Circle) and isinstance(torus.v_iso(0.0), Circle)
    assert isinstance(patch.u_iso(1.0), BSplineCurve) and isinstance(patch.v_iso(0.0), BSplineCurve)


def test_iso_curve_errors_are_value_errors():
    with pytest.raises(ValueError):
        Sphere(Point.origin(), 1.0).u_iso(math.nan)
    with pytest.raises(ValueError):
        Sphere(Point.origin(), 2.0).to_bspline_patch((0.0, 3.0), (-1.0, 1.0)).u_iso(10.0)


def test_a_pole_iso_is_a_circle_of_radius_zero():
    pole = Sphere(Point.origin(), 1.0).v_iso(math.pi / 2)
    assert pole.radius() == pytest.approx(0.0, abs=1e-12)
    assert pole.center().z == pytest.approx(1.0)


def _cubic_curve_2d():
    from geomcore import BSplineCurve2D

    return BSplineCurve2D(
        3,
        [Point2D(0.0, 0.0), Point2D(1.0, 2.0), Point2D(3.0, 2.0), Point2D(4.0, 0.0)],
        [0.0, 1.0],
        [4, 4],
        False,
    )


def test_bspline_curve_2d_matches_its_3d_twin():
    planar = _cubic_curve_2d()
    spatial = BSplineCurve(
        3,
        [Point(p.x, p.y, 0.0) for p in planar.poles()],
        [0.0, 1.0],
        [4, 4],
        False,
    )
    assert planar.degree() == 3 and not planar.is_rational()
    assert planar.bounds() == (0.0, 1.0) and planar.period() is None
    for u in [0.0, 0.25, 0.5, 0.9, 1.0]:
        p, q = planar.point_at(u), spatial.point_at(u)
        assert (p.x, p.y) == pytest.approx((q.x, q.y), abs=1e-14)
        _, d1 = planar.first_derivatives(u)
        _, e1 = spatial.first_derivatives(u)
        assert (d1.x, d1.y) == pytest.approx((e1.x, e1.y), abs=1e-14)
    assert planar.length(1e-9) == pytest.approx(spatial.length(1e-9), rel=1e-12)
    left, right = planar.split_at(0.5)
    assert left.bounds() == (0.0, 0.5) and right.bounds() == (0.5, 1.0)
    assert planar.reversed().point_at(0.0).x == pytest.approx(4.0)
    assert len(planar.insert_knot(0.5, 1).poles()) == 5


def test_bspline_curve_2d_projection_and_discretization():
    curve = _cubic_curve_2d()
    on = curve.point_at(0.3)
    assert curve.contains(on)
    assert curve.parameter_of(on) == pytest.approx(0.3)
    assert curve.parameter_of(on, Tolerance.default()) == pytest.approx(0.3)
    projection = curve.project_point(Point2D(2.0, 3.0))
    assert projection.point.x == pytest.approx(2.0)
    assert projection.distance == pytest.approx(3.0 - curve.point_at(0.5).y)
    nan = curve.project_point(Point2D(math.nan, 0.0))
    assert math.isnan(nan.parameter) and math.isnan(nan.distance)
    d = curve.discretize(1e-3)
    assert isinstance(d, geomcore.Discretization2D)
    assert d.parameters[0] == 0.0 and d.parameters[-1] == 1.0
    assert all(dev <= 1e-3 for dev in d.chord_deviations)
    assert curve.point_count(1e-3) == len(d.points)


def test_bspline_curve_2d_intersections_are_numeric_and_complete():
    from geomcore import Circle2D, IntersectionMethod, TrimmedCurve2D

    curve = _cubic_curve_2d()
    line = Line2D(Point2D(0.0, 1.0), Vector2D(1.0, 0.0))
    result = curve.intersect(line)
    assert result.method == IntersectionMethod.NUMERIC and result.is_complete()
    assert len(result.hits) == 2
    for hit in result.hits:
        assert hit.point.y == pytest.approx(1.0)
        assert curve.point_at(hit.first_parameter).x == pytest.approx(hit.point.x)
    # The same crossings, seen from the line.
    assert len(line.intersect(curve).hits) == 2

    circle = Circle2D(Point2D(0.0, 0.0), 2.0)
    spline = circle.to_bspline()
    assert spline.is_rational() and spline.is_closed()
    other = Circle2D(Point2D(3.0, 0.0), 2.0).to_bspline()
    hits = spline.intersect(other).hits
    assert sorted(h.point.y for h in hits) == pytest.approx([-math.sqrt(1.75), math.sqrt(1.75)])

    arc = TrimmedCurve2D(circle, 0.0, math.pi)
    half = arc.to_bspline()
    assert half.point_at(half.bounds()[1]).x == pytest.approx(-2.0)
