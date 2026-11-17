"""The least-squares fits (#253): the points they take, the tolerance, the
result classes and the errors. The numerics are tested in Rust."""

import math

import numpy
import pytest

import geomcore
from geomcore import (
    Circle,
    Circle2D,
    CircleFit,
    CircleFit2D,
    Cylinder,
    CylinderFit,
    FitError,
    Plane,
    PlaneFit,
    Point,
    Point2D,
    Sphere,
    SphereFit,
    Tolerance,
)


def _ring(count=8, radius=2.0, height=0.0):
    return [
        (radius * math.cos(k * math.tau / count), radius * math.sin(k * math.tau / count), height)
        for k in range(count)
    ]


def _tube():
    return [p for height in (0.0, 1.0, 2.0) for p in _ring(6, 1.5, height)]


FITS = [
    (Plane, PlaneFit, _ring),
    (Circle, CircleFit, _ring),
    (Sphere, SphereFit, lambda: _ring() + [(0.0, 0.0, 2.0)]),
    (Cylinder, CylinderFit, _tube),
]


@pytest.mark.parametrize("shape, result, points", FITS, ids=lambda o: getattr(o, "__name__", ""))
def test_a_fit_takes_points_tuples_or_an_array_and_returns_its_result(shape, result, points):
    as_points = [Point(*p) for p in points()]
    fits = [
        shape.fit(as_points),
        shape.fit(points()),
        shape.fit(numpy.array(points())),
        shape.fit(as_points, Tolerance.default()),
        shape.fit(as_points, tol=Tolerance(1e-9, 1e-12)),
    ]
    for fit in fits:
        assert type(fit) is result
        assert type(fit.shape) is shape
        assert isinstance(fit.root_mean_square_residual, float)
        assert isinstance(fit.maximum_residual, float)
        assert fit == fits[0]


def test_a_circle_in_the_plane_takes_2d_points():
    points = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)]
    fit = Circle2D.fit(numpy.array(points))
    assert type(fit) is CircleFit2D and type(fit.shape) is Circle2D
    assert fit == Circle2D.fit([Point2D(*p) for p in points])
    assert fit.shape.radius() == pytest.approx(1.0)
    assert fit.maximum_residual == pytest.approx(0.0, abs=1e-14)


def test_the_residuals_are_read_from_the_fit():
    # Every point 0.25 off the circle of radius 1, alternately in and out.
    points = [
        Point2D((1.0 + 0.25 * (-1) ** k) * math.cos(k * math.tau / 8),
                (1.0 + 0.25 * (-1) ** k) * math.sin(k * math.tau / 8))
        for k in range(8)
    ]
    fit = Circle2D.fit(points)
    assert fit.root_mean_square_residual == pytest.approx(0.25)
    assert fit.maximum_residual == pytest.approx(0.25)


def test_a_fit_repr_names_its_shape_and_residuals():
    fit = Circle2D.fit([(3.0, 0.0), (0.0, 3.0), (-3.0, 0.0)])
    assert repr(fit).startswith("CircleFit2D(shape=Circle2D(")
    assert repr(fit).endswith(", root_mean_square_residual=0.0, maximum_residual=0.0)")
    assert repr(fit.shape) in repr(fit)


@pytest.mark.parametrize(
    "call, message",
    [
        (lambda: Plane.fit([(0, 0, 0), (1, 0, 0)]), "the fit needs at least 3 points, 2 were given"),
        (lambda: Sphere.fit(_ring()), "the points are coplanar"),
        (lambda: Cylinder.fit(_ring()), "the points are coplanar"),
        (lambda: Circle.fit([(0, 0, 0), (1, 1, 1), (2, 2, 2)]), "the points are collinear"),
        (lambda: Circle2D.fit([(1, 1)] * 3), "the points are confused"),
        (lambda: Plane.fit([(0, 0, 0), (1, 0, 0), (0, math.nan, 0)]), "NaN or infinite"),
    ],
)
def test_a_refused_fit_raises_a_fit_error(call, message):
    with pytest.raises(FitError, match=message):
        call()


def test_a_fit_error_is_a_construction_error():
    assert issubclass(FitError, geomcore.ConstructionError)
    assert issubclass(FitError, ValueError)


def test_a_fit_rejects_points_of_the_wrong_kind():
    with pytest.raises(TypeError):
        Plane.fit([Point2D(0, 0), Point2D(1, 0), Point2D(0, 1)])
    with pytest.raises(ValueError):
        Plane.fit(numpy.zeros((4, 2)))
