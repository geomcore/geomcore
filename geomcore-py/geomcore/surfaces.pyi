# Deprecated: the surfaces are at the package root, `from geomcore import Plane`.
# This module goes in the next release.

from geomcore import (
    BSplineSurface as BSplineSurface,
    BezierSurface as BezierSurface,
    Cone as Cone,
    Cylinder as Cylinder,
    ExtrusionSurface as ExtrusionSurface,
    Plane as Plane,
    RevolutionSurface as RevolutionSurface,
    Sphere as Sphere,
    Torus as Torus,
)

__all__ = [
    "BSplineSurface",
    "BezierSurface",
    "Cone",
    "Cylinder",
    "ExtrusionSurface",
    "Plane",
    "RevolutionSurface",
    "Sphere",
    "Torus",
]
