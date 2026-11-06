from typing import Union as _Union

from geomcore._native import *
from geomcore import _native

__doc__ = _native.__doc__
__version__ = _native.__version__

# Every 3D curve, 2D curve and surface class, as the union to annotate
# with. `isinstance` takes a union from Python 3.10; `typing.get_args` gives
# the tuple of classes.
Curve = _Union[Line, Circle, Ellipse, Parabola, Hyperbola, BSplineCurve, BezierCurve, TrimmedCurve]
Curve2D = _Union[
    Line2D,
    Circle2D,
    Ellipse2D,
    Parabola2D,
    Hyperbola2D,
    BSplineCurve2D,
    BezierCurve2D,
    TrimmedCurve2D,
]
Surface = _Union[
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    ExtrusionSurface,
    RevolutionSurface,
    BSplineSurface,
    BezierSurface,
]

__all__ = [*_native.__all__, "Curve", "Curve2D", "Surface"]


def __getattr__(name):
    # `geomcore.curves` and `geomcore.surfaces` without importing them
    # first, as before the shapes moved to the root.
    if name in ("curves", "surfaces"):
        import importlib

        return importlib.import_module(f"geomcore.{name}")
    raise AttributeError(f"module 'geomcore' has no attribute {name!r}")
