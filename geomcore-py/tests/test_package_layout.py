"""One path per name, at the package root (#203): the shapes sit beside the
values and results, `Curve`, `Curve2D` and `Surface` name the classes of
the Rust enums, a trimmed curve or surface being a piece of one, and the old
`geomcore.curves` / `geomcore.surfaces` paths warn."""

import pickle
import sys
from typing import get_args

import pytest

import geomcore
from geomcore import (
    Circle,
    Curve,
    Curve2D,
    Line2D,
    Plane,
    Point,
    Surface,
    TrimmedCurve,
    TrimmedCurve2D,
    TrimmedSurface,
    Vector,
)


def _shape_classes_in_all():
    names = (getattr(geomcore, name) for name in geomcore.__all__)
    return {cls for cls in names if isinstance(cls, type) and hasattr(cls, "intersect")}


def test_curve_curve2d_and_surface_hold_every_shape_class_but_the_trimmed_once():
    kinds = [get_args(kind) for kind in (Curve, Curve2D, Surface)]
    kinds.append((TrimmedCurve, TrimmedCurve2D, TrimmedSurface))
    assert sum(len(kind) for kind in kinds) == len(_shape_classes_in_all())
    assert set().union(*kinds) == _shape_classes_in_all()
    assert all(cls.__name__.endswith("2D") for cls in get_args(Curve2D))
    assert not any(cls.__name__.endswith("2D") for cls in get_args(Curve))


def test_isinstance_takes_the_tuple_of_classes_on_every_python():
    circle = Circle(Point.origin(), Vector.unit_z(), 1.0)
    assert isinstance(circle, get_args(Curve))
    assert not isinstance(circle, get_args(Surface))


@pytest.mark.skipif(sys.version_info < (3, 10), reason="isinstance takes a Union from Python 3.10")
def test_isinstance_takes_the_kind_names():
    circle = Circle(Point.origin(), Vector.unit_z(), 1.0)
    assert isinstance(circle, Curve)
    # A trimmed curve is a piece of a curve, as in Rust, not a variant.
    assert not isinstance(TrimmedCurve(circle, 0.0, 1.0), Curve)
    assert not isinstance(circle, Surface)
    assert isinstance(Plane(Point.origin(), Vector.unit_z()), Surface)
    assert isinstance(Line2D.from_two_points(geomcore.Point2D(0, 0), geomcore.Point2D(1, 0)), Curve2D)


def test_every_class_lives_at_the_root():
    for name in geomcore.__all__:
        item = getattr(geomcore, name)
        if isinstance(item, type):
            assert item.__module__ == "geomcore", name
    circle = Circle(Point.origin(), Vector.unit_z(), 2.0)
    assert pickle.loads(pickle.dumps(circle)) == circle


def test_the_native_module_is_private():
    assert "geomcore.geomcore" not in sys.modules
    assert not hasattr(geomcore, "geomcore")
    assert geomcore._native.Circle is Circle
    assert geomcore.__version__ == geomcore._native.__version__


def test_the_old_submodule_paths_warn_and_give_the_root_class():
    with pytest.warns(DeprecationWarning, match=r"geomcore\.curves\.Circle is deprecated: use geomcore\.Circle"):
        from geomcore.curves import Circle as old_circle
    assert old_circle is Circle
    with pytest.warns(DeprecationWarning, match=r"geomcore\.surfaces\.Plane"):
        assert geomcore.surfaces.Plane is Plane
    with pytest.raises(AttributeError):
        geomcore.curves.Point  # noqa: B018
    with pytest.raises(AttributeError):
        geomcore.nothing  # noqa: B018
