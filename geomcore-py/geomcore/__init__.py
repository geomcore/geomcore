import warnings as _warnings

from .geomcore import *

__doc__ = geomcore.__doc__
if hasattr(geomcore, "__all__"):
    __all__ = geomcore.__all__

# The 0.36 names, deprecated in 0.37 (ADR 0013): only 2D types carry a
# dimension suffix. Kept for one release so callers migrate with warnings.
_RENAMED = {
    "Point3D": "Point",
    "Vector3D": "Vector",
    "Axis3D": "Axis",
    "Frame3D": "Frame",
    "BoundingBox3D": "BoundingBox",
    "Discretization3D": "Discretization",
}
_RENAMED_CURVES = {
    "Line3D": "Line",
    "Circle3D": "Circle",
    "Ellipse3D": "Ellipse",
    "Parabola3D": "Parabola",
    "Hyperbola3D": "Hyperbola",
    "BSplineCurve3D": "BSplineCurve",
    "TrimmedCurve3D": "TrimmedCurve",
}


def _renamed(module, renamed, module_name):
    def __getattr__(name):
        if name in renamed:
            _warnings.warn(
                f"{module_name}.{name} is renamed {renamed[name]}",
                DeprecationWarning,
                stacklevel=2,
            )
            return getattr(module, renamed[name])
        raise AttributeError(f"module {module_name!r} has no attribute {name!r}")

    return __getattr__


__getattr__ = _renamed(geomcore, _RENAMED, "geomcore")
curves.__getattr__ = _renamed(curves, _RENAMED_CURVES, "geomcore.curves")
