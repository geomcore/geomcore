"""The 0.36 names still resolve in 0.37, with a DeprecationWarning (ADR 0013)."""

import pytest

import geomcore
import geomcore.curves


def test_old_value_name_warns_and_resolves():
    with pytest.warns(DeprecationWarning, match="renamed Point"):
        assert geomcore.Point3D is geomcore.Point


def test_old_curve_name_warns_and_resolves():
    with pytest.warns(DeprecationWarning, match="renamed Circle"):
        from geomcore.curves import Circle3D
    assert Circle3D is geomcore.curves.Circle


def test_unknown_name_still_raises():
    with pytest.raises(AttributeError):
        geomcore.NoSuchType
