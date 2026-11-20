# A 2D curve offsets exactly or within a deviation, reporting its cusps

The offset of a 2D curve `c` at the signed distance `d` is the curve `c(t) + d n(t)`, `n` the unit left normal, the tangent turned a quarter counterclockwise. A positive distance goes to the side a positive signed curvature turns toward, so a counterclockwise circle shrinks. Tool-path compensation needs it (#250). It is an operation, as a surface's offset is (ADR 0016):

```rust
pub trait Offset2D: ParametricCurve2D {
    type Offset;
    fn offset(&self, distance: f64, max_deviation: f64, tol: Tolerance)
        -> Result<CurveOffset2D<Self::Offset>, CurveOffset2DError>;
}
```

It is implemented for every 2D curve type, `Curve2D` and `TrimmedCurve2D<C>`, in `operations/offset/`, and Python has it in the 2D curve method set.

## Exact where the offset stays in its class

- A **line** offsets to the line through `origin + d n`, with the same direction, so the same parameter. A whole line is unbounded and its offset is a whole line.
- A **circle** offsets to the circle of the same frame, so the same parameter, of radius `r − d` when it runs counterclockwise (its centre on its left) and `r + d` when it runs clockwise. A radius at or below `tol.confusion` is refused with `NonPositiveRadius { radius }`. Below zero the formula still gives a circle, run backwards, but a tool path never wants it.
- A **trimmed curve** offsets its basis over its range and trims the result to the same range: a segment to a segment, an arc to an arc, across the seam too. A hidden routing hook, `offset_over(bounds, …)`, gives the basis's offset over the piece's range.

An exact offset reports a deviation of `0`, proven, and `is_exact()`.

## Within a deviation otherwise

The offset of an ellipse, a parabola, a hyperbola, a Bézier curve or a B-spline is no curve of its class, nor rational. It comes back as one non-rational B-spline in the source curve's parameter, within `max_deviation` of `c(t) + d n(t)` at the same `t`. The same-parameter deviation bounds the Hausdorff distance both ways, and keeps the cusps' parameters meaningful on the result.

It reuses the Hermite spans of ADR 0035, moved from `uv_curve/approximation.rs` to `geometry/approximation/hermite_spans.rs` beside the result types, so that a curve operation does not reach into the surfaces module. The offset is one more `SpanSource`, its 2D jets carried in the `x` and `y` of the 3D ones as `curve_3d_to_uv_within` does. The jets come from the curve's derivatives by Taylor series: with `V(s)` the series of `c'`, the unit tangent is `V (V · V)^(−1/2)`, the power rule giving the series of the inverse square root, so an order-`n` jet takes the curve's derivatives to order `n + 1`. Every curve evaluates every order. The degree follows ADR 0035: degree 9 along a curve with no knot, unless degree 5 takes no more poles, and degree 5 along a spline with knots. The normal takes one derivative, so at a knot where the spline is `Cᵏ` the offset is `Cᵏ⁻¹`, with one-sided jets. At a `C⁰` knot the offset stays continuous if the tangent keeps its direction within `tol.angular`, and is refused with `UndefinedNormal` where it turns a corner. So is a point where the speed vanishes, by the rule of `tangent`.

The deviation is sampled, as `curve_3d_to_uv_within`'s is: each span is checked at nine points against half the budget, and the result says it is not proven. No closed-form bound on the offset's high derivatives exists for a conic or a spline, so a proven bound would need interval arithmetic on the series. The offset of a whole parabola or hyperbola is unbounded and has no finite B-spline, so it is refused with `Unbounded`: offset a trimmed piece.

## Cusps are reported, loops kept

Where `|d|` crosses the radius of curvature on its side, `1 − d κ` changes sign: the offset's speed, `|c'| (1 − d κ)`, passes through zero, and the offset stops and turns back. The result lists those parameters on the source curve, `cusp_parameters()`, in increasing order. Each is bracketed where the sign of `|c'|³ − d (c' × c'')` changes between the breakpoints and the span samples, then bisected to the precision of a double. The parametric offset is smooth through a cusp, so the spans need no special treatment there. A root where `1 − d κ` touches zero without changing sign is no cusp and is not reported. Between two cusps the offset runs against the curve and forms a loop; trimming loops needs the offset's self-intersections and is left to the caller, who knows which side of the tool path to keep.

## The result

`CurveOffset2D<C = Curve2D>` holds the curve, its deviation, whether it is proven and the cusp parameters. Its curve type is the trait's `Offset`: `Line2D`, `Circle2D`, `BSplineCurve2D`, `Curve2D` for a `Curve2D`, and `TrimmedCurve2D<C::Offset>` for a trimmed curve.

## Considered options

- **`CurveApproximation2D` as the result** (#199): rejected. Its curve is a B-spline by type, and the offset of a whole line is no B-spline. It also has no place for the cusps. `CurveOffset2D` keeps its getters (`curve`, `deviation`, `is_proven`), so a caller reads both alike.
- **One result curve type, `Curve2D`, for every input**: rejected. A segment's offset would lose its range, as `Curve2D` has no trimmed variant, and a caller offsetting a `Line2D` would match on a variant it knows.
- **Exact offsets in a separate method**: rejected. One name returning one shape is simpler for a caller offsetting mixed pieces of a profile, and `is_exact()` says which it got.
- **Trimming the loops**: not now (YAGNI). The first version reports the cusps; trimming needs the offset's self-intersections, and which side to keep depends on the use.
- **A proven deviation**: not now. It needs bounds on derivatives of order 6 or 10 of `c + d n`, which grow with the curvature near a cusp; sampling is what ADR 0035 does for every pair without a closed form.
- **Refusing every unbounded curve, a whole line too**: rejected. A line's offset is exact and unbounded like the line, so there is nothing to approximate.

## Consequences

`hermite_spans.rs` holds the spans, their join and the choice between degree 9 and degree 5 that `uv_curve_to_3d_within` made inline, `sampled_spans_of_fewest_poles`, and the 2D join `curve_3d_to_uv_within` wrote out, `joined_2d`. Both mappings call them unchanged.

`tests/similarity_invariance.rs` compares every 2D sample's offset, whole and trimmed, under the judge's similarities (ADR 0051): the outcome, the degree, the pole count, the cusp count and positions, and points of the offset curve.

A 3D offset of a planar curve in a given plane, if a caller needs it, would read the curve in that plane and reuse this one.
