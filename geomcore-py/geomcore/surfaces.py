"""Deprecated: the surfaces are at the package root, `from geomcore import Plane`.

Each name still resolves here, with a DeprecationWarning, until the next
release removes this module.
"""

import warnings
from typing import get_args

import geomcore

__all__ = sorted(cls.__name__ for cls in get_args(geomcore.Surface))


def __getattr__(name):
    if name in __all__:
        warnings.warn(
            f"geomcore.surfaces.{name} is deprecated: use geomcore.{name}",
            DeprecationWarning,
            stacklevel=2,
        )
        return getattr(geomcore, name)
    raise AttributeError(f"module 'geomcore.surfaces' has no attribute {name!r}")
