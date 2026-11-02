"""Every curve and surface class answers its tier's shared method set and
pairs with every other class.

The geometry is checked elsewhere; this file checks the bindings' shape:
that each class carries the methods its core trait gives it, that the
stubs type and document them as the runtime has them, that `intersect`
accepts every class on either side, and that each class survives the
conversions through the core enums.

The method sets are read from the stubs, where each is a stub-only class
(`_CurveMethods`, `_DiscretizeMethods`, ...) that a class inherits.
"""

import ast
import inspect
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

import geomcore
from geomcore import (
    Axis,
    CurveCurveIntersection,
    CurveCurveIntersection2D,
    CurveSurfaceIntersection,
    Frame,
    Point2D,
    Point,
    SurfaceSurfaceIntersection,
    Transform,
    Vector2D,
    Vector,
)
from geomcore import curves, surfaces

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
        surfaces.ExtrusionSurface(
            curves.Ellipse(Point(0.0, 0.0, -0.5), Vector.unit_z(), Vector.unit_x(), 0.6, 0.4),
            Vector(0.1, 0.0, 1.0),
        ),
        surfaces.RevolutionSurface(
            curves.Circle.from_frame(
                Frame(Point(0.7, 0.0, 0.0), Vector.unit_y(), Vector.unit_x()), 0.25
            ),
            Axis(Point.origin(), Vector.unit_z()),
        ),
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


def _shape_classes():
    return {**_classes(curves), **_classes(surfaces)}


class _StubClass:
    """A class of the curve or surface stubs, as written there."""

    def __init__(self, node):
        functions = [item for item in node.body if isinstance(item, ast.FunctionDef)]
        self.bases = [_base_name(base) for base in node.bases]
        self.methods = [function.name for function in functions]
        # An overloaded method carries its docstring on its first overload.
        self.docstrings = {}
        for function in functions:
            docstring = ast.get_docstring(function)
            if docstring is not None:
                self.docstrings.setdefault(function.name, docstring)


def _base_name(base):
    return base.value.id if isinstance(base, ast.Subscript) else base.id


def _read_stub_classes():
    stub_directory = Path(geomcore.__file__).parent
    classes = {}
    for module in ("curves", "surfaces"):
        tree = ast.parse((stub_directory / f"{module}.pyi").read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, ast.ClassDef):
                classes[node.name] = _StubClass(node)
    return classes


STUB_CLASSES = _read_stub_classes()


def _stub_bases(name):
    """`name` and the stub classes it inherits, nearest first."""
    yield STUB_CLASSES[name]
    for base in STUB_CLASSES[name].bases:
        if base in STUB_CLASSES:
            yield from _stub_bases(base)


def _stub_methods(name):
    """The methods the stub of class `name` gives it, inherited ones included."""
    return {method for stub in _stub_bases(name) for method in stub.methods}


def _stub_docstrings(name):
    docstrings = {}
    for stub in _stub_bases(name):
        for method, docstring in stub.docstrings.items():
            docstrings.setdefault(method, docstring)
    return docstrings


# ---------------------------------------------------------------- method sets


def test_every_class_has_a_fixture():
    covered = {type(shape).__name__ for shape in _operands_3d() + _operands_2d()}
    assert set(_classes(curves)) | set(_classes(surfaces)) == covered


def test_every_class_takes_its_tier_method_set():
    for name in _shape_classes():
        if name in _classes(surfaces):
            tier = "_SurfaceMethods"
        elif name.endswith("2D"):
            tier = "_Curve2DMethods"
        else:
            tier = "_CurveMethods"
        assert tier in STUB_CLASSES[name].bases, name


@pytest.mark.parametrize("name", sorted(_shape_classes()))
def test_every_class_answers_each_method_its_stub_gives_it(name):
    # stubtest checks the signature of every method the stub gives a class,
    # but not that the runtime class has each method set the stub inherits.
    assert _stub_methods(name) <= set(dir(_shape_classes()[name]))


def test_3d_curves_are_transformed_and_2d_curves_are_not():
    for name, cls in _classes(curves).items():
        assert ("transformed" in dir(cls)) == (not name.endswith("2D")), name


def test_bounded_curves_discretize_and_unbounded_ones_do_not():
    discretize = _stub_methods("_DiscretizeMethods")
    for curve in _operands_3d()[:9] + _operands_2d():
        bounded = all(abs(end) != float("inf") for end in curve.bounds())
        has = discretize <= set(dir(curve))
        assert has == bounded, type(curve).__name__


def test_spline_classes_have_the_sided_evaluations():
    sided = {
        name
        for module in (curves, surfaces)
        for name, cls in _classes(module).items()
        if _stub_methods("_SidedEvaluationMethods") <= set(dir(cls))
    }
    assert sided == {
        "BSplineCurve",
        "BSplineCurve2D",
        "BezierCurve",
        "BezierCurve2D",
        "BSplineSurface",
        "BezierSurface",
    }


# ---------------------------------------------------------------- stubs


def test_every_method_set_method_has_a_docstring():
    for name, stub in STUB_CLASSES.items():
        if name.endswith("Methods"):
            missing = sorted(set(stub.methods) - set(stub.docstrings))
            assert not missing, (name, missing)


@pytest.mark.parametrize("name", sorted(_shape_classes()))
def test_stub_docstrings_are_the_runtime_ones(name):
    runtime = _shape_classes()[name]
    for method, docstring in _stub_docstrings(name).items():
        assert docstring == inspect.getdoc(getattr(runtime, method)), (name, method)


def test_type_checkers_see_every_method_of_every_class(tmp_path):
    mypy = pytest.importorskip("mypy.api")
    lines = [
        "from typing_extensions import assert_type",
        "import geomcore",
        "from geomcore import curves, surfaces",
    ]
    for module in (curves, surfaces):
        for name, cls in sorted(_classes(module).items()):
            for method in sorted(m for m in vars(cls) if not m.startswith("_")):
                lines.append(f"{module.__name__}.{name}.{method}")
    lines += [
        "def check(",
        "    line: curves.Line,",
        "    line_2d: curves.Line2D,",
        "    ellipse: curves.Ellipse,",
        "    spline_2d: curves.BSplineCurve2D,",
        "    plane: surfaces.Plane,",
        "    torus: surfaces.Torus,",
        ") -> None:",
        "    assert_type(line.eval_point(0.0), geomcore.Point)",
        "    assert_type(line_2d.eval_point(0.0), geomcore.Point2D)",
        "    assert_type(line.intersect(line), geomcore.CurveCurveIntersection)",
        "    assert_type(line.intersect(plane), geomcore.CurveSurfaceIntersection)",
        "    assert_type(line_2d.intersect(spline_2d), geomcore.CurveCurveIntersection2D)",
        "    assert_type(plane.curve_3d_to_uv(ellipse), curves._WholeCurve2D)",
        "    assert_type(ellipse.discretize(0.1), geomcore.Discretization)",
        "    assert_type(spline_2d.discretize(0.1), geomcore.Discretization2D)",
        "    assert_type(spline_2d.eval_d1_sided(0.0, geomcore.Side.BELOW)[1], geomcore.Vector2D)",
        "    assert_type(plane.u_iso(0.0), curves.Line)",
        "    assert_type(torus.v_iso(0.0), curves.Circle)",
        "    assert_type(torus.intersect(plane), geomcore.SurfaceSurfaceIntersection)",
    ]
    source = tmp_path / "uses_the_stubs.py"
    source.write_text("\n".join(lines) + "\n")
    report, errors, status = mypy.run(
        ["--strict", "--no-incremental", "--cache-dir", str(tmp_path / "cache"), str(source)]
    )
    assert status == 0, report + errors


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
            surface.curve_3d_to_uv(line)
        except geomcore.GeomcoreError:
            pass  # accepted, but the line does not lie on the surface


def _rigid_motion_with_mirror():
    axis = Axis(Point(0.5, -1.0, 0.0), Vector(1.0, 2.0, 2.0))
    mirror = Transform.mirror_plane(
        Frame(Point(0.1, 0.2, 0.3), Vector(1.0, 1.0, 0.0), Vector(0.0, 0.0, 1.0))
    )
    return Transform.rotation(axis, 0.7).then(mirror)


def test_every_shape_is_transformed_into_its_class_and_commutes_with_evaluation():
    t = _rigid_motion_with_mirror()
    for curve in _operands_3d()[:9]:
        moved = curve.transformed(t)
        assert type(moved) is type(curve)
        u = _middle(curve.bounds())
        assert moved.eval_point(u).distance(t.apply_point(curve.eval_point(u))) < 1e-9
    for surface in _surfaces():
        moved = surface.transformed(t)
        assert type(moved) is type(surface)
        u, v = _middle(surface.u_bounds()), _middle(surface.v_bounds())
        image = t.apply_point(surface.eval_point(u, v))
        assert moved.eval_point(u, v).distance(image) < 1e-9, type(surface).__name__


def test_every_shape_rejects_a_collapsing_transform():
    collapse = Transform.scaling(Point.origin(), 0.0)
    for shape in _operands_3d():
        with pytest.raises(geomcore.ConstructionError):
            shape.transformed(collapse)


def test_every_curve_reversed_meets_itself_at_the_reversed_parameter():
    for curve in _operands_3d()[:9] + _operands_2d():
        reversed_curve = curve.reversed()
        assert type(reversed_curve) is type(curve)
        ends = sorted(curve.reversed_parameter(end) for end in curve.bounds())
        assert reversed_curve.bounds() == tuple(ends)
        first, last = (max(end, -1.0) if end < 0 else min(end, 1.0) for end in curve.bounds())
        u = first + 0.3 * (last - first)
        there = reversed_curve.eval_point(curve.reversed_parameter(u))
        assert there.distance(curve.eval_point(u)) < 1e-9, type(curve).__name__
        assert reversed_curve.reversed() == curve
