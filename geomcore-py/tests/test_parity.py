"""Every curve and surface class answers its tier's shared method set and
pairs with every other class.

The geometry is checked elsewhere; this file checks the bindings' shape:
that each class carries the methods its core trait gives it, that the
stubs type and document them as the runtime has them, that `intersect`
accepts every class on either side, and that each class survives the
conversions through the core enums.

The method sets are read from the stubs, where each is a stub-only class
(`_CurveMethods`, `_BoundedCurveMethods`, ...) that a class inherits.
"""

import ast
import importlib.util
import inspect
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import get_args

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

def _curves_3d():
    """One small instance of each whole 3D curve class, near the origin."""
    return [
        geomcore.Line(Point(0.0, 0.0, -2.0), Vector(0.1, 0.2, 1.0)),
        geomcore.Circle(Point.origin(), Vector.unit_z(), 1.0),
        geomcore.Ellipse(Point.origin(), Vector.unit_y(), Vector.unit_x(), 1.5, 0.5),
        geomcore.Parabola(Point.origin(), Vector.unit_x(), Vector.unit_y(), 0.5),
        geomcore.Hyperbola(Point.origin(), Vector.unit_z(), Vector.unit_x(), 1.0, 0.5),
        geomcore.BSplineCurve(
            2,
            [Point(-1.0, -1.0, 0.5), Point(0.0, 2.0, -0.5), Point(1.0, -1.0, 0.5)],
            [0.0, 1.0],
            [3, 3],
            False,
        ),
        geomcore.BezierCurve(
            2, [Point(-1.0, 0.5, -1.0), Point(0.0, 0.5, 2.0), Point(1.0, 0.5, -1.0)]
        ),
        geomcore.QuarticBranch(
            geomcore.Cone(Point.origin(), Vector.unit_z(), 1.0, 0.5),
            geomcore.Cone(Point(0.0, 0.3, 0.0), Vector.unit_x(), 1.0, 0.5),
            0,
        ),
    ]


def _curves_2d():
    """One small instance of each whole 2D curve class, near the origin."""
    return [
        geomcore.Line2D(Point2D(-2.0, 0.1), Vector2D(1.0, 0.05)),
        geomcore.Circle2D(Point2D(0.0, 0.0), 1.0),
        geomcore.Ellipse2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0), 1.5, 0.5),
        geomcore.Parabola2D(Point2D(0.0, -0.5), Vector2D(0.0, 1.0), 0.5),
        geomcore.Hyperbola2D(Point2D(0.0, 0.0), Vector2D(1.0, 0.0), 0.5, 0.5),
        geomcore.BSplineCurve2D(
            2,
            [Point2D(-1.0, -1.0), Point2D(0.0, 2.0), Point2D(1.0, -1.0)],
            [0.0, 1.0],
            [3, 3],
            False,
        ),
        geomcore.BezierCurve2D(2, [Point2D(-1.0, 0.5), Point2D(0.0, -1.0), Point2D(1.0, 0.5)]),
    ]


def _surfaces():
    """One small instance of each surface class, near the origin."""
    patch = [
        [Point(-1.0, -1.0, 0.2), Point(-1.0, 1.0, -0.2)],
        [Point(1.0, -1.0, -0.2), Point(1.0, 1.0, 0.2)],
    ]
    return [
        geomcore.Plane(Point(0.0, 0.0, 0.1), Vector(0.0, 0.1, 1.0)),
        geomcore.Cylinder(Point.origin(), Vector.unit_y(), 0.5),
        geomcore.Cone(Point(0.0, 0.0, -1.0), Vector.unit_z(), 0.4, 0.5),
        geomcore.Sphere(Point(0.2, 0.0, 0.0), 0.8),
        geomcore.Torus(Point.origin(), Vector.unit_x(), 1.0, 0.3),
        geomcore.ExtrusionSurface(
            geomcore.Ellipse(Point(0.0, 0.0, -0.5), Vector.unit_z(), Vector.unit_x(), 0.6, 0.4),
            Vector(0.1, 0.0, 1.0),
        ),
        geomcore.RevolutionSurface(
            geomcore.Circle.from_frame(
                Frame(Point(0.7, 0.0, 0.0), Vector.unit_y(), Vector.unit_x()), 0.25
            ),
            Axis(Point.origin(), Vector.unit_z()),
        ),
        geomcore.BSplineSurface(
            1, 1, patch, [0.0, 1.0], [2, 2], [0.0, 1.0], [2, 2], False, False
        ),
        geomcore.BezierSurface(1, 1, patch),
    ]


def _piece(curve, trimmed_class):
    """A trimmed piece of `curve` over a finite part of its range."""
    first, last = curve.bounds()
    if first == float("-inf") or last == float("inf"):
        first, last = -1.0, 1.0
    return trimmed_class(curve, first, first + 0.75 * (last - first))


def _operands_3d():
    whole = _curves_3d()
    return whole + [_piece(c, geomcore.TrimmedCurve) for c in whole[:2]] + _surfaces()


def _operands_2d():
    whole = _curves_2d()
    return whole + [_piece(c, geomcore.TrimmedCurve2D) for c in whole[:2]]


_CURVE_CLASSES = {
    cls.__name__: cls
    for cls in (
        *get_args(geomcore.Curve),
        geomcore.TrimmedCurve,
        *get_args(geomcore.Curve2D),
        geomcore.TrimmedCurve2D,
    )
}
_SURFACE_CLASSES = {
    cls.__name__: cls for cls in (*get_args(geomcore.Surface), geomcore.TrimmedSurface)
}


def _shape_classes():
    return {**_CURVE_CLASSES, **_SURFACE_CLASSES}


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
    if isinstance(base, ast.Subscript):
        return base.value.id
    return base.attr if isinstance(base, ast.Attribute) else base.id


def _read_stub_classes():
    stub_directory = Path(geomcore.__file__).parent
    classes = {}
    tree = ast.parse((stub_directory / "__init__.pyi").read_text(encoding="utf-8"))
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
    assert set(_shape_classes()) == set(_fixtures())


def test_every_class_takes_its_tier_method_set():
    for name in _shape_classes():
        if name == "TrimmedSurface":
            # A patch has its own methods: the set's isoparametric curves,
            # reversals and transforms are its basis's.
            continue
        if name in _SURFACE_CLASSES:
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
    for name, cls in _CURVE_CLASSES.items():
        assert ("transformed" in dir(cls)) == (not name.endswith("2D")), name


def test_bounded_curves_have_ends_and_discretize_and_unbounded_ones_do_not():
    bounded_methods = _stub_methods("_BoundedCurveMethods")
    for curve in _operands_3d()[:9] + _operands_2d():
        bounded = all(abs(end) != float("inf") for end in curve.bounds())
        has = bounded_methods <= set(dir(curve))
        assert has == bounded, type(curve).__name__


def test_spline_classes_have_the_sided_evaluations():
    sided = {
        name
        for name, cls in _shape_classes().items()
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


def test_every_runtime_docstring_is_in_the_stub():
    # Editors that read only the stub show its docstrings on hover (#334).
    path = Path(__file__).resolve().parents[2] / "scripts" / "sync_stub_docstrings.py"
    spec = importlib.util.spec_from_file_location("sync_stub_docstrings", path)
    sync = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sync)
    edits = sync.stub_edits(sync.STUB.read_text(encoding="utf-8"))
    assert not edits, (
        f"{len(edits)} stub docstrings differ from the runtime ones, the first at "
        f"line {edits[0][0][0]}: run scripts/sync_stub_docstrings.py"
    )


def _runtime_callables():
    for name in dir(geomcore):
        value = getattr(geomcore, name)
        if inspect.isclass(value):
            yield name, value
            for member in vars(value):
                if not member.startswith("__") or member == "__new__":
                    yield f"{name}.{member}", getattr(value, member)
        elif callable(value):
            yield name, value


def test_no_default_prints_as_an_ellipsis():
    # `help()` and `inspect.signature` show the value a parameter takes (#334).
    for name, value in _runtime_callables():
        try:
            signature = inspect.signature(value)
        except (TypeError, ValueError):
            continue
        for parameter in signature.parameters.values():
            assert parameter.default is not Ellipsis, (name, parameter.name)
    stub = (Path(geomcore.__file__).parent / "__init__.pyi").read_text(encoding="utf-8")
    assert "= ..." not in stub


def test_the_printed_tolerance_defaults_are_the_default_tolerance():
    defaults = {
        name: parameter.default
        for name, parameter in inspect.signature(geomcore.Tolerance).parameters.items()
    }
    assert geomcore.Tolerance(**defaults) == geomcore.Tolerance.default()
    assert defaults["confusion"] == inspect.signature(
        geomcore.BSplineCurve.length
    ).parameters["max_length_error"].default


def _stub_trees():
    stub_directory = Path(geomcore.__file__).parent
    return {
        module: ast.parse((stub_directory / f"{module}.pyi").read_text(encoding="utf-8"))
        for module in ("__init__",)
    }


def _union_members(value):
    """The names a union annotation lists, `Union[A, B]` or `A | B`, or
    None when it is not a union of plain names. A method set such as
    `_SurfaceMethods[Any, Any]` counts as its plain name."""
    if (
        isinstance(value, ast.Subscript)
        and isinstance(value.value, ast.Name)
        and value.value.id in STUB_CLASSES
    ):
        return [value.value.id]
    if (
        isinstance(value, ast.Subscript)
        and isinstance(value.value, ast.Name)
        and value.value.id == "Union"
    ):
        elements = value.slice.elts if isinstance(value.slice, ast.Tuple) else [value.slice]
    elif isinstance(value, ast.BinOp) and isinstance(value.op, ast.BitOr):
        left, right = _union_members(value.left), _union_members(value.right)
        if left is None or right is None:
            return None
        return left + right
    elif isinstance(value, ast.Name):
        return [value.id]
    else:
        return None
    if not all(isinstance(element, ast.Name) for element in elements):
        return None
    return [element.id for element in elements]


def _union_aliases(trees):
    """Each union alias of the stubs whose members are all curve or surface
    classes, such as `Curve` or `Surface`, as its class names."""
    members = {}
    for tree in trees.values():
        for node in tree.body:
            value = getattr(node, "value", None)
            if not isinstance(node, ast.Assign) or not isinstance(
                value, (ast.Subscript, ast.BinOp)
            ):
                continue
            names = _union_members(value)
            if names is not None:
                members[node.targets[0].id] = names

    def expand(name):
        return [
            leaf
            for member in members[name]
            for leaf in (expand(member) if member in members else [member])
        ]

    shapes = _shape_classes()
    expanded = {name: expand(name) for name in members}
    return {
        name: leaves
        for name, leaves in expanded.items()
        if all(leaf in shapes for leaf in leaves)
    }


def _union_parameters():
    """Every stub parameter typed with a union alias: (owner, function, the
    function's required parameters as (name, annotation), the union one)."""
    trees = _stub_trees()
    aliases = _union_aliases(trees)
    shapes = _shape_classes()
    # A parameter typed with a method set, such as `_CurveMethods`, takes
    # every class that inherits it.
    for base in STUB_CLASSES:
        inheritors = [
            name
            for name in shapes
            if name != base and any(stub is STUB_CLASSES[base] for stub in _stub_bases(name))
        ]
        if base.startswith("_") and inheritors:
            aliases[base] = inheritors
    found = []

    def annotation_key(annotation):
        # A plain name, or an inline union of shapes such as
        # `Curve | TrimmedCurve`, registered as an alias of its own.
        names = None if annotation is None else _union_members(annotation)
        if names is None:
            # Another annotation, as `Tuple[float, float]`, by its text.
            return None if annotation is None else ast.unparse(annotation)
        if len(names) == 1:
            return names[0]
        leaves = [leaf for name in names for leaf in aliases.get(name, [name])]
        if not all(leaf in shapes for leaf in leaves):
            return None
        key = " | ".join(names)
        aliases[key] = leaves
        return key

    def scan(owner, function):
        arguments = function.args.args
        required = arguments[: len(arguments) - len(function.args.defaults)]
        if owner is not None and not _is_static(function):
            required = required[1:]  # self or cls
        typed = [(argument.arg, annotation_key(argument.annotation)) for argument in required]
        for index, (_, annotation) in enumerate(typed):
            if annotation in aliases:
                found.append((owner, function.name, typed, index))

    for module, tree in trees.items():
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and module == "__init__":
                scan(None, node)
            elif isinstance(node, ast.ClassDef):
                for item in node.body:
                    if isinstance(item, ast.FunctionDef):
                        scan(node.name, item)
    return aliases, found


def _is_static(function):
    names = [getattr(decorator, "id", None) for decorator in function.decorator_list]
    return "staticmethod" in names


UNION_ALIASES, UNION_PARAMETERS = _union_parameters()


def _fixtures():
    """One instance per curve and surface class, by class name."""
    fixtures = {}
    for shape in _operands_3d() + _operands_2d():
        fixtures.setdefault(type(shape).__name__, shape)
    # A patch intersects surfaces only, so it is no operand of the
    # curve-surface parity checks.
    fixtures["TrimmedSurface"] = fixtures["Cylinder"].trimmed((0.0, 3.0), (-1.0, 1.0))
    return fixtures


def _receiver(owner, fixtures):
    """An instance the methods of stub class `owner` run on: its own
    fixture, or that of the first class inheriting a method set."""
    if owner in fixtures:
        return fixtures[owner]
    return next(
        shape
        for name, shape in fixtures.items()
        if any(stub is STUB_CLASSES[owner] for stub in _stub_bases(name))
    )


_OTHER_ARGUMENTS = {
    "Vector": Vector(0.1, 0.0, 1.0),
    "Axis": Axis(Point(0.0, 0.0, 5.0), Vector.unit_x()),
    "float": 0.1,
    "int": 8,
    "Tuple[float, float]": (0.0, 1.0),
}


@pytest.mark.parametrize(
    "owner, function, typed, index",
    UNION_PARAMETERS,
    ids=[
        f"{owner or 'geomcore'}.{function}({typed[index][0]})"
        for owner, function, typed, index in UNION_PARAMETERS
    ],
)
def test_every_member_of_a_stub_union_is_accepted(owner, function, typed, index):
    # stubtest checks that a parameter exists with a compatible annotation,
    # not that the runtime takes each class of a `Union`: #176 found three
    # parameters typed for any curve that raised TypeError on a
    # `TrimmedCurve`.
    fixtures = _fixtures()
    if owner is None:
        call = getattr(geomcore, function)
    elif function == "__new__":
        call = _shape_classes()[owner]
    else:
        call = getattr(_receiver(owner, fixtures), function)

    def argument(annotation):
        if annotation in UNION_ALIASES:
            return fixtures[UNION_ALIASES[annotation][0]]
        return _OTHER_ARGUMENTS[annotation]

    for member in UNION_ALIASES[typed[index][1]]:
        arguments = [argument(annotation) for _, annotation in typed]
        arguments[index] = fixtures[member]
        try:
            call(*arguments)
        except TypeError as error:
            pytest.fail(f"{member}: {error}")
        except geomcore.GeomcoreError:
            pass  # accepted, but these operands have no answer


def test_type_checkers_see_every_method_of_every_class(tmp_path):
    mypy = pytest.importorskip("mypy.api")
    lines = [
        "from typing_extensions import assert_type",
        "import geomcore",
    ]
    for name, cls in sorted(_shape_classes().items()):
        for method in sorted(m for m in vars(cls) if not m.startswith("_")):
            lines.append(f"geomcore.{name}.{method}")
    lines += [
        "def check(",
        "    line: geomcore.Line,",
        "    line_2d: geomcore.Line2D,",
        "    ellipse: geomcore.Ellipse,",
        "    spline_2d: geomcore.BSplineCurve2D,",
        "    plane: geomcore.Plane,",
        "    torus: geomcore.Torus,",
        ") -> None:",
        "    assert_type(line.point_at(0.0), geomcore.Point)",
        "    assert_type(line_2d.point_at(0.0), geomcore.Point2D)",
        "    assert_type(line.intersect(line), geomcore.CurveCurveIntersection)",
        "    assert_type(line.intersect(plane), geomcore.CurveSurfaceIntersection)",
        "    assert_type(line_2d.intersect(spline_2d), geomcore.CurveCurveIntersection2D)",
        "    assert_type(plane.curve_3d_to_uv(ellipse), geomcore.Curve2D)",
        "    assert_type(ellipse.discretize(0.1), geomcore.Discretization)",
        "    assert_type(spline_2d.discretize(0.1), geomcore.Discretization2D)",
        "    assert_type(spline_2d.first_derivatives_sided(0.0, geomcore.Side.BELOW)[1], geomcore.Vector2D)",
        "    assert_type(plane.u_iso(0.0), geomcore.Line)",
        "    assert_type(torus.v_iso(0.0), geomcore.Circle)",
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
    if isinstance(a, tuple(_SURFACE_CLASSES.values())):
        if isinstance(b, tuple(_SURFACE_CLASSES.values())):
            return SurfaceSurfaceIntersection
        return CurveSurfaceIntersection
    if isinstance(b, tuple(_SURFACE_CLASSES.values())):
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
        basis = _piece(curve, geomcore.TrimmedCurve).basis()
        assert type(basis) is type(curve) and basis == curve
    for curve in _curves_2d():
        basis = _piece(curve, geomcore.TrimmedCurve2D).basis()
        assert type(basis) is type(curve) and basis == curve


def _middle(bounds):
    first, last = (max(end, -1.0) if end < 0 else min(end, 1.0) for end in bounds)
    return (first + last) / 2


def test_every_surface_gives_its_isoparametric_curves():
    for surface in _surfaces():
        u, v = _middle(surface.u_bounds()), _middle(surface.v_bounds())
        on_surface = surface.point_at(u, v)
        for curve, parameter in ((surface.u_iso(u), v), (surface.v_iso(v), u)):
            assert type(curve).__name__ in _CURVE_CLASSES, type(surface).__name__
            assert curve.point_at(parameter).distance(on_surface) < 1e-9


def test_every_surface_is_accepted_where_a_surface_is_expected():
    line = geomcore.Line(Point(0.0, 0.0, -2.0), Vector.unit_z())
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
        assert moved.point_at(u).distance(t.apply_point(curve.point_at(u))) < 1e-9
    for surface in _surfaces():
        moved = surface.transformed(t)
        assert type(moved) is type(surface)
        u, v = _middle(surface.u_bounds()), _middle(surface.v_bounds())
        image = t.apply_point(surface.point_at(u, v))
        assert moved.point_at(u, v).distance(image) < 1e-9, type(surface).__name__


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
        there = reversed_curve.point_at(curve.reversed_parameter(u))
        assert there.distance(curve.point_at(u)) < 1e-9, type(curve).__name__
        assert reversed_curve.reversed() == curve


def test_every_surface_says_whether_it_is_periodic_in_each_direction():
    # #197: `is_u_periodic` was a B-spline surface method; `u_period` was shared.
    for surface in _surfaces():
        assert surface.is_u_periodic() == (surface.u_period() is not None), surface
        assert surface.is_v_periodic() == (surface.v_period() is not None), surface


def test_every_spline_converts_to_a_bspline_like_its_siblings():
    # #197: generic code special-cased the B-spline classes, which had no
    # `to_bspline`, when converting a curve or surface.
    shapes = _curves_3d() + _curves_2d() + _surfaces()
    for shape in shapes:
        bspline = (geomcore.BSplineCurve, geomcore.BSplineCurve2D, geomcore.BSplineSurface)
        if isinstance(shape, bspline):
            assert shape.to_bspline() == shape
    splines = ("BSplineCurve", "BSplineCurve2D", "BezierCurve", "BezierCurve2D")
    for name in splines + ("BSplineSurface", "BezierSurface"):
        assert "to_bspline" in dir(_shape_classes()[name]), name


# The classes whose repr lists the values that define the shape, so its
# keywords are the constructor's. The spline, trimmed and swept reprs are
# summaries (pole counts, the basis, the profile) and are not checked.
_ELEMENTARY = {
    "Line", "Circle", "Ellipse", "Parabola", "Hyperbola",
    "Plane", "Cylinder", "Cone", "Sphere", "Torus",
    "Line2D", "Circle2D", "Ellipse2D", "Parabola2D", "Hyperbola2D",
}


def test_elementary_reprs_use_the_constructor_keywords():
    seen = set()
    for shape in _operands_3d() + _operands_2d():
        name = type(shape).__name__
        if name not in _ELEMENTARY:
            continue
        seen.add(name)
        keywords = re.findall(r"(\w+)=", repr(shape))
        parameters = list(inspect.signature(type(shape)).parameters)
        # In the constructor's order; a repr may leave out a direction
        # (the x direction of a conic) that the frame derives.
        assert keywords == [p for p in parameters if p in keywords], (repr(shape), parameters)
    assert seen == _ELEMENTARY
