"""Every intersection, projection, extrema, curvature, approximation and
discretization result pickles with the operands it holds, and copies
(#318)."""

import concurrent.futures
import copy
import math
import multiprocessing
import pickle

import pytest

import geomcore
from geomcore import (
    BSplineCurve,
    BSplineCurve2D,
    BSplineSurface,
    Circle,
    Circle2D,
    Line,
    Line2D,
    Plane,
    Point,
    Point2D,
    Sphere,
    Vector,
    Vector2D,
    curve_curve_extrema,
)


def _wave() -> BSplineSurface:
    """A bicubic patch over [0, 4]², waved in z."""
    poles = [
        [Point(4 * i / 3, 4 * j / 3, math.sin(4 * i / 3) * math.cos(4 * j / 3)) for j in range(4)]
        for i in range(4)
    ]
    return BSplineSurface(3, 3, poles, [0.0, 1.0], [4, 4], [0.0, 1.0], [4, 4], False, False)


def _results() -> list:
    circle = Circle(Point(0, 0, 0), Vector(0, 0, 1), 1.0)
    line = Line(Point(0, 0, 0), Vector(1, 0, 0))
    sphere = Sphere(Point(0, 0, 0), 2.0)
    wave = _wave()
    curve = wave.v_iso(0.3)
    arch = [Point(i, math.sin(i / 3), 0) for i in range(20)]
    grid = [[Point(i, j, math.sin(i / 3)) for j in range(8)] for i in range(8)]
    tilted = Plane(Point(0, 0, 0.1), Vector(0.1, 0.2, 1))
    results = [
        circle.intersect(line),
        Circle2D(Point2D(0, 0), 1.0).intersect(Line2D(Point2D(0, 0), Vector2D(1, 0))),
        line.intersect(sphere),
        Plane(Point(0, 0, 0), Vector(0, 0, 1)).intersect(sphere),
        # Numeric results: B-spline curves with their curves in parameters,
        # and incomplete searches with regions of each family.
        wave.intersect(tilted),
        curve.extract(0.2, 0.7).intersect(curve),
        wave.u_iso(0.4).intersect(wave),
        wave.intersect(wave.extract((0.2, 0.6), (0.2, 0.6))),
        circle.project_point(Point(3, 0, 0)),
        Circle2D(Point2D(0, 0), 1.0).project_point(Point2D(3, 0)),
        sphere.project_point(Point(0, 0, 5)),
        curve_curve_extrema(circle, Line(Point(0, 0, 3), Vector(1, 0, 0))),
        sphere.curvatures(0.3, 0.4),
        BSplineCurve.approximate_within(arch, 3, 1e-3),
        BSplineCurve2D.approximate_within([Point2D(p.x, p.y) for p in arch], 3, 1e-3),
        BSplineSurface.approximate_within(grid, 3, 3, 1e-3),
        circle.discretize(1e-2),
        Circle2D(Point2D(0, 0), 1.0).discretize(1e-2),
    ]
    parts = []
    for result in results:
        for name in ("hits", "overlaps", "curves", "minima"):
            parts.extend(getattr(result, name, []))
        completeness = getattr(result, "completeness", None)
        if completeness is not None:
            parts.append(completeness)
            parts.extend(completeness.regions)
    return results + parts


_RESULTS = _results()


@pytest.mark.parametrize("original", _RESULTS, ids=lambda o: type(o).__name__)
def test_every_result_pickles_to_an_equal_object(original):
    for protocol in range(2, pickle.HIGHEST_PROTOCOL + 1):
        restored = pickle.loads(pickle.dumps(original, protocol=protocol))
        assert type(restored) is type(original)
        assert restored == original
        assert repr(restored) == repr(original)
    assert copy.deepcopy(original) == original
    assert copy.copy(original) == original


def test_every_result_class_is_covered():
    names = {type(result).__name__ for result in _RESULTS}
    result_classes = {
        name
        for name in geomcore.__all__
        if name.endswith(
            (
                "Intersection",
                "Intersection2D",
                "Projection",
                "Projection2D",
                "Hit",
                "Hit2D",
                "Overlap",
                "Region",
                "Extrema",
                "Extremum",
                "Approximation",
                "Approximation2D",
            )
        )
        or name in ("Completeness", "IntersectionCurve", "SurfaceCurvatures")
        or name.startswith("Discretization") and "Criteria" not in name and "Error" not in name
    }
    assert result_classes <= names, sorted(result_classes - names)
    for name in ("Completeness", "CurveCurveRegion", "CurveSurfaceRegion", "SurfaceSurfaceRegion"):
        assert name in names


def test_a_pickled_result_keeps_its_operands():
    circle = Circle(Point(0, 0, 0), Vector(0, 0, 1), 1.0)
    line = Line(Point(0, 0, 0), Vector(1, 0, 0))
    restored = pickle.loads(pickle.dumps(circle.intersect(line)))
    assert restored.first == circle and restored.second == line
    assert restored.points == circle.intersect(line).points


def test_a_result_pickled_by_another_version_is_refused_with_both_versions_named():
    hit = Circle(Point(0, 0, 0), Vector(0, 0, 1), 1.0).intersect(Line(Point(0, 0, 0), Vector(1, 0, 0))).hits[0]
    constructor, (version, state) = hit.__reduce__()
    assert version == geomcore.__version__
    with pytest.raises(ValueError, match=f"geomcore 0.0.1 and geomcore {geomcore.__version__}"):
        constructor("0.0.1", state)


def _intersect_in_a_worker(radius: float):
    circle = Circle(Point(0, 0, 0), Vector(0, 0, 1), radius)
    return circle.intersect(Line(Point(0, 0, 0), Vector(1, 0, 0)))


def test_a_worker_returns_an_intersection():
    context = multiprocessing.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(1, mp_context=context) as pool:
        result = pool.submit(_intersect_in_a_worker, 2.0).result()
    assert result == _intersect_in_a_worker(2.0)
