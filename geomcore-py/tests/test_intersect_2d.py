"""`a.intersect(b)` on every 2D pair: exact for lines and conics, numeric
for a B-spline, in either order and through `TrimmedCurve2D`."""

import itertools
import math

import pytest

from geomcore import IntersectionMethod, Point2D, Vector2D
from geomcore.curves import (
    BSplineCurve2D,
    Circle2D,
    Ellipse2D,
    Hyperbola2D,
    Line2D,
    Parabola2D,
    TrimmedCurve2D,
)

ORIGIN = Point2D.origin()
X, Y = Vector2D.unit_x(), Vector2D.unit_y()


def line():
    return Line2D(Point2D(0.0, 1.0), X)


def circle():
    return Circle2D(ORIGIN, 2.5)


def ellipse():
    return Ellipse2D(ORIGIN, X, 3.0, 2.0)


def parabola():
    return Parabola2D(ORIGIN, X, 1.0)


def hyperbola():
    return Hyperbola2D(ORIGIN, X, 1.0, 0.5)


def spline():
    return Ellipse2D(ORIGIN, Y, 3.0, 2.0).to_bspline()


EXACT = {
    "line": line,
    "circle": circle,
    "ellipse": ellipse,
    "parabola": parabola,
    "hyperbola": hyperbola,
}


def points(result):
    return sorted((round(h.point.x, 7), round(h.point.y, 7)) for h in result.hits)


@pytest.mark.parametrize(
    ("first", "second"), list(itertools.product(EXACT, EXACT)), ids=lambda name: name
)
def test_every_pair_of_lines_and_conics_is_exact_in_both_orders(first, second):
    a, b = EXACT[first](), EXACT[second]()
    forward, backward = a.intersect(b), b.intersect(a)
    assert forward.method == IntersectionMethod.EXACT
    assert backward.method == IntersectionMethod.EXACT
    assert forward.is_complete() and backward.is_complete()
    assert points(forward) == points(backward)
    assert sorted(h.multiplicity for h in forward.hits) == sorted(
        h.multiplicity for h in backward.hits
    )


def test_line_meets_each_conic():
    r3 = 3.0 * math.sqrt(3.0) / 2.0
    result = line().intersect(ellipse())
    assert points(result) == [(round(-r3, 7), 1.0), (round(r3, 7), 1.0)]
    vertical = Line2D(Point2D(4.0, 0.0), Y)
    assert points(vertical.intersect(parabola())) == [(4.0, -4.0), (4.0, 4.0)]
    two = Line2D(Point2D(2.0, 0.0), Y).intersect(Hyperbola2D(ORIGIN, X, 1.0, 1.0))
    root = math.sqrt(3.0)
    assert points(two) == [(2.0, round(-root, 7)), (2.0, round(root, 7))]


def test_tangencies_have_multiplicity_two():
    tangent = Line2D(Point2D(0.0, 2.0), X).intersect(ellipse())
    assert [h.multiplicity for h in tangent.hits] == [2]
    touching = Circle2D(ORIGIN, 2.0).intersect(ellipse())
    assert [h.multiplicity for h in touching.hits] == [2, 2]
    apex = Line2D(ORIGIN, Y).intersect(parabola())
    assert [h.multiplicity for h in apex.hits] == [2]


def test_two_conics_meet_up_to_four_times():
    a = Ellipse2D(ORIGIN, X, 3.0, 1.0)
    b = Ellipse2D(ORIGIN, Y, 3.0, 1.0)
    result = a.intersect(b)
    t = 3.0 / math.sqrt(10.0)
    expected = sorted((round(sx * t, 7), round(sy * t, 7)) for sx in (-1, 1) for sy in (-1, 1))
    assert points(result) == expected
    assert all(h.multiplicity == 1 for h in result.hits)


def test_coincident_conics_overlap():
    result = ellipse().intersect(ellipse())
    assert result.method == IntersectionMethod.EXACT
    assert len(result.overlaps) == 1
    overlap = result.overlaps[0]
    assert overlap.first_range == pytest.approx((0.0, math.tau))
    assert overlap.same_direction

    a = TrimmedCurve2D(parabola(), -2.0, 2.0)
    b = TrimmedCurve2D(parabola(), 0.0, 5.0)
    shared = a.intersect(b).overlaps
    assert len(shared) == 1
    assert shared[0].first_range == pytest.approx((0.0, 2.0))
    assert shared[0].second_range == pytest.approx((0.0, 2.0))


def test_trimmed_conics_clip_in_both_orders():
    vertical = Line2D(Point2D(1.5, 0.0), Y)
    upper = TrimmedCurve2D(ellipse(), 0.0, math.pi)
    root = math.sqrt(3.0)
    for result in (vertical.intersect(upper), upper.intersect(vertical)):
        assert result.method == IntersectionMethod.EXACT
        assert points(result) == [(1.5, round(root, 7))]
    arc = TrimmedCurve2D(hyperbola(), 0.0, 2.0)
    assert arc.intersect(upper).is_complete()
    assert upper.intersect(arc).is_complete()


@pytest.mark.parametrize("name", list(EXACT))
def test_bspline_against_every_exact_curve_is_numeric_in_both_orders(name):
    curve = EXACT[name]()
    if name in ("parabola", "hyperbola"):
        curve = TrimmedCurve2D(curve, -6.0, 6.0)
    forward, backward = spline().intersect(curve), curve.intersect(spline())
    assert forward.method == IntersectionMethod.NUMERIC
    assert backward.method == IntersectionMethod.NUMERIC
    assert forward.is_complete() and backward.is_complete()
    assert points(forward) == points(backward)


def test_bspline_against_bspline_and_a_trimmed_one():
    a = Circle2D(ORIGIN, 1.0).to_bspline()
    b = Circle2D(Point2D(1.0, 0.0), 1.0).to_bspline()
    result = a.intersect(b)
    h = math.sqrt(0.75)
    assert result.method == IntersectionMethod.NUMERIC
    assert result.is_complete()
    assert points(result) == [(0.5, round(-h, 7)), (0.5, round(h, 7))]
    half = TrimmedCurve2D(Ellipse2D(ORIGIN, X, 3.0, 2.0).to_bspline(), 0.0, math.pi)
    cut = half.intersect(Line2D(Point2D(1.5, 0.0), Y))
    assert cut.method == IntersectionMethod.NUMERIC
    assert points(cut) == [(1.5, round(math.sqrt(3.0), 7))]


def test_miss_is_empty_and_complete():
    far = Ellipse2D(Point2D(10.0, 0.0), X, 3.0, 2.0)
    result = ellipse().intersect(far)
    assert result.is_empty()
    assert not result
    assert result.is_complete()
