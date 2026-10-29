"""Every curve and surface class answers its tier's shared method set and
pairs with every other class.

The geometry is checked elsewhere; this file checks the bindings' shape:
that each class carries the methods its core trait gives it, that
`intersect` accepts every class on either side, and that each class
survives the conversions through the core enums.
"""

import inspect
from concurrent.futures import ThreadPoolExecutor

import pytest

import geomcore
from geomcore import (
    CurveCurveIntersection,
    CurveCurveIntersection2D,
    CurveSurfaceIntersection,
    Point2D,
    Point,
    SurfaceSurfaceIntersection,
    Vector2D,
    Vector,
)
from geomcore import curves, surfaces

CURVE_METHODS = {
    "eval_point",
    "eval_points",
    "eval_derivative",
    "eval_d1",
    "eval_d2",
    "tangent",
    "curvature",
    "curvature_center",
    "bounds",
    "period",
    "is_periodic",
    "is_closed",
    "continuity",
    "min_continuity",
    "extrema",
    "project_point",
    "project_points",
    "contains",
    "parameter_of",
    "bounding_box",
    "intersect",
}

SURFACE_METHODS = {
    "eval_point",
    "eval_points",
    "eval_derivative",
    "eval_d1",
    "eval_d2",
    "u_bounds",
    "v_bounds",
    "u_period",
    "v_period",
    "is_u_closed",
    "is_v_closed",
    "continuity_u",
    "continuity_v",
    "min_continuity_u",
    "min_continuity_v",
    "curvatures",
    "project_point",
    "project_points",
    "contains",
    "parameters_of",
    "bounding_box",
    "intersect",
    "u_iso",
    "v_iso",
    "reversed_u",
    "reversed_v",
    "u_reversed_parameter",
    "v_reversed_parameter",
}

BOUNDED_CURVE_METHODS = {"discretize", "point_count"}
SIDED_METHODS = {"eval_d1_sided", "eval_d2_sided", "eval_derivative_sided"}


def _curves_3d():
    """One small instance of each whole 3D curve class, near the origin."""
    return [
        curves.Line(Point(0.0, 0.0, -2.0), Vector(0.1, 0.2, 1.0)),
        curves.Circle(Point.origin(), Vector.unit_z(), 1.0),
        curves.Ellipse(Point.origin(), Vector.unit_y(), Vector.unit_x(), 1.5, 0.5),
        curves.Parabola(Point.origin(), Vector.unit_x(), Vector.unit_y(), 0.5),
        curves.Hyperbola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 1.0, 0.5),
        curves.BSplineCurve(
            2,
            [Point(-1.0, -1.0, 0.5), Point(0.0, 2.0, -0.5), Point(1.0, -1.0, 0.5)],
            [0.0, 1.0],
            [3, 3],
            False,
        ),
        curves.BezierCurve(
            2, [Point(-1.0, 0.5, -1.0), Point(0.0, 0.5, 2.0), Point(1.0, 0.5, -1.0)]
        ),
    ]


def _curves_2d():
    """One small instance of each whole 2D curve class, near the origin."""
    return [
        curves.Line2D(Point2D(-2.0, 0.1), Vector2D(1.0, 0.05)),
        curves.Circle2D(Point2D(0.0, 0.0), 1.0),
        curves.Ellipse2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0), 1.5, 0.5),
        curves.Parabola2D(Point2D(0.0, -0.5), Vector2D(0.0, 1.0), 0.5),
        curves.Hyperbola2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0), 0.5, 0.5),
        curves.BSplineCurve2D(
            2,
            [Point2D(-1.0, -1.0), Point2D(0.0, 2.0), Point2D(1.0, -1.0)],
            [0.0, 1.0],
            [3, 3],
            False,
        ),
        curves.BezierCurve2D(2, [Point2D(-1.0, 0.5), Point2D(0.0, -1.0), Point2D(1.0, 0.5)]),
    ]


def _surfaces():
    """One small instance of each surface class, near the origin."""
    patch = [
        [Point(-1.0, -1.0, 0.2), Point(-1.0, 1.0, -0.2)],
        [Point(1.0, -1.0, -0.2), Point(1.0, 1.0, 0.2)],
    ]
    return [
        surfaces.Plane(Point(0.0, 0.0, 0.1), Vector(0.0, 0.1, 1.0)),
        surfaces.Cylinder(Point.origin(), Vector.unit_y(), 0.5),
        surfaces.Cone(Point(0.0, 0.0, -1.0), Vector.unit_z(), 0.4, 0.5),
        surfaces.Sphere(Point(0.2, 0.0, 0.0), 0.8),
        surfaces.Torus(Point.origin(), Vector.unit_x(), 1.0, 0.3),
        surfaces.BSplineSurface(
            1, 1, patch, [0.0, 1.0], [2, 2], [0.0, 1.0], [2, 2], False, False
        ),
        surfaces.BezierSurface(1, 1, patch),
    ]


def _piece(curve, trimmed_class):
    """A trimmed piece of `curve` over a finite part of its range."""
    first, last = curve.bounds()
    if first == float("-inf") or last == float("inf"):
        first, last = -1.0, 1.0
    return trimmed_class(curve, first, first + 0.75 * (last - first))


def _operands_3d():
    whole = _curves_3d()
    return whole + [_piece(c, curves.TrimmedCurve) for c in whole[:2]] + _surfaces()


def _operands_2d():
    whole = _curves_2d()
    return whole + [_piece(c, curves.TrimmedCurve2D) for c in whole[:2]]


def _classes(module):
    return {
        name: value
        for name, value in vars(module).items()
        if inspect.isclass(value) and not name.startswith("_")
    }


# ---------------------------------------------------------------- method sets


def test_every_class_has_a_fixture():
    covered = {type(shape).__name__ for shape in _operands_3d() + _operands_2d()}
    assert set(_classes(curves)) | set(_classes(surfaces)) == covered


@pytest.mark.parametrize("name", sorted(_classes(curves)))
def test_every_curve_class_has_the_curve_methods(name):
    assert CURVE_METHODS <= set(dir(_classes(curves)[name]))


@pytest.mark.parametrize("name", sorted(_classes(surfaces)))
def test_every_surface_class_has_the_surface_methods(name):
    assert SURFACE_METHODS <= set(dir(_classes(surfaces)[name]))


def test_bounded_curves_discretize_and_unbounded_ones_do_not():
    for curve in _operands_3d()[:9] + _operands_2d():
        bounded = all(abs(end) != float("inf") for end in curve.bounds())
        has = BOUNDED_CURVE_METHODS <= set(dir(curve))
        assert has == bounded, type(curve).__name__


def test_spline_classes_have_the_sided_evaluations():
    sided = {
        name
        for module in (curves, surfaces)
        for name, cls in _classes(module).items()
        if SIDED_METHODS <= set(dir(cls))
    }
    assert sided == {
        "BSplineCurve",
        "BSplineCurve2D",
        "BezierCurve",
        "BezierCurve2D",
        "BSplineSurface",
        "BezierSurface",
    }


# ---------------------------------------------------------------- intersect


def _family(a, b):
    if isinstance(a, tuple(_classes(surfaces).values())):
        if isinstance(b, tuple(_classes(surfaces).values())):
            return SurfaceSurfaceIntersection
        return CurveSurfaceIntersection
    if isinstance(b, tuple(_classes(surfaces).values())):
        return CurveSurfaceIntersection
    return CurveCurveIntersection


def test_every_3d_operand_intersects_every_other():
    operands = _operands_3d()
    for a in operands:
        for b in operands:
            result = a.intersect(b)
            assert isinstance(result, _family(a, b)), (a, b)
            assert result.first is a and result.second is b
            assert result == a.intersect(b)
            assert repr(result).startswith(type(result).__name__ + "(")


def test_every_2d_operand_intersects_every_other():
    operands = _operands_2d()
    for a in operands:
        for b in operands:
            result = a.intersect(b)
            assert isinstance(result, CurveCurveIntersection2D), (a, b)
            assert result.first is a and result.second is b
            assert result == a.intersect(b)


def test_2d_curves_do_not_intersect_3d_shapes():
    for a in _operands_2d():
        for b in _operands_3d():
            with pytest.raises(TypeError):
                a.intersect(b)
            with pytest.raises(TypeError):
                b.intersect(a)


def test_curve_surface_results_name_the_operands_by_role():
    for surface in _surfaces():
        for curve in _curves_3d():
            assert surface.intersect(curve).surface is surface
            assert curve.intersect(surface).curve is curve


def test_2d_bspline_intersections_run_from_several_threads():
    spline = _curves_2d()[5]
    others = [_curves_2d()[0], _curves_2d()[1], _curves_2d()[6]]
    expected = [spline.intersect(other) for other in others]
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(lambda k: spline.intersect(others[k % 3]), range(24)))
    assert results == [expected[k % 3] for k in range(24)]


# ---------------------------------------------------------------- conversions


def test_every_class_has_a_repr_and_value_equality():
    for shape in _operands_3d() + _operands_2d():
        assert repr(shape).startswith(type(shape).__name__ + "(")
        assert shape == shape
    for shapes in (_operands_3d(), _operands_2d()):
        for a, b in zip(shapes, shapes[1:]):
            assert a != b


def test_every_whole_curve_comes_back_as_its_class():
    for curve in _curves_3d():
        basis = _piece(curve, curves.TrimmedCurve).basis()
        assert type(basis) is type(curve) and basis == curve
    for curve in _curves_2d():
        basis = _piece(curve, curves.TrimmedCurve2D).basis()
        assert type(basis) is type(curve) and basis == curve


def _middle(bounds):
    first, last = (max(end, -1.0) if end < 0 else min(end, 1.0) for end in bounds)
    return (first + last) / 2


def test_every_surface_gives_its_isoparametric_curves():
    for surface in _surfaces():
        u, v = _middle(surface.u_bounds()), _middle(surface.v_bounds())
        on_surface = surface.eval_point(u, v)
        for curve, parameter in ((surface.u_iso(u), v), (surface.v_iso(v), u)):
            assert type(curve).__name__ in _classes(curves), type(surface).__name__
            assert curve.eval_point(parameter).distance(on_surface) < 1e-9


def test_every_surface_is_accepted_where_a_surface_is_expected():
    line = curves.Line(Point(0.0, 0.0, -2.0), Vector.unit_z())
    for surface in _surfaces():
        try:
            line.parametrize_on(surface)
        except geomcore.GeomcoreError:
            pass  # accepted, but the line does not lie on the surface
