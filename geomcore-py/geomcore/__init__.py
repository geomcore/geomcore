"""Curves and surfaces for CAD, in Rust: intersect, project, fit, evaluate.

    >>> from geomcore import BSplineCurve, Circle, Line, Plane, Point, Vector
    >>> circle = Circle(Point(0, 0, 0), Vector(0, 0, 1), 2.0)
    >>> hits = circle.intersect(Plane(Point(1, 0, 0), Vector(1, 0, 0)))
    >>> print(len(hits), hits.method)
    2 IntersectionMethod.EXACT
    >>> round(circle.project_point(Point(3, 0, 0)).distance, 9)
    1.0
    >>> circle.extrema(Line(Point(0, 0, 3), Vector(1, 0, 0))).closest().distance
    3.0
    >>> points = [Point(0, 0, 0), Point(1, 2, 0), Point(3, 2, 1), Point(4, 0, 0)]
    >>> BSplineCurve.interpolate(points, 3).degree()
    3

Every result says how it was found (`method`) and whether it is complete
(`completeness`). `Curve`, `Curve2D` and `Surface` are the unions to
annotate with. Documentation: https://geomcore.github.io/geomcore/
"""

from typing import Union as _Union

from geomcore._native import *
from geomcore import _native

__version__ = _native.__version__

# The union to annotate with of each kind of shape. A curve is one of the
# classes of the `Curve` enum's variants, as in Rust; a `TrimmedCurve` is a
# piece of one, so code that takes either says `Curve | TrimmedCurve`.
# `isinstance` takes the union; `typing.get_args` gives the tuple of classes.
Curve = _Union[
    Line, Circle, Ellipse, Parabola, Hyperbola, BSplineCurve, BezierCurve, QuarticBranch
]
Curve2D = _Union[
    Line2D,
    Circle2D,
    Ellipse2D,
    Parabola2D,
    Hyperbola2D,
    BSplineCurve2D,
    BezierCurve2D,
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

# `help(geomcore.Curve)` shows these rather than `typing.Union`'s own
# text, where the union takes a docstring: from Python 3.14 a union is a
# `types.UnionType`, whose docstring is read-only.
for _union, _text in (
    (Curve, "Any 3D curve but a TrimmedCurve: annotate with `Curve | TrimmedCurve` to take both."),
    (Curve2D, "Any 2D curve but a TrimmedCurve2D."),
    (Surface, "Any surface but a TrimmedSurface."),
):
    try:
        _union.__doc__ = _text
    except (AttributeError, TypeError):
        pass
del _union, _text

__all__ = [*_native.__all__, "Curve", "Curve2D", "Surface"]

