# A UV curve on a spline surface is composed exactly

`uv_curve_to_3d` (ADR 0030) answered a spline surface only along `u` or `v`. Any other UV curve whose parameter is polynomial or rational now gives its exact 3D curve on a B-spline or Bézier surface, or on the extrusion of a B-spline or Bézier profile: the composition `S∘c`, as a B-spline in the UV curve's own parameter. The code is `surfaces/uv_curve/composition.rs`. The pair table reaches it through the `composed` hook of `IsoCurves` in `to_curve.rs`, once the existing rows (the affine map, then lines and straight splines along `u` or `v`) have refused with `NotAnalytic`. Those rows answer first, so an iso line keeps its path and its result.

## Scope

| UV curve | Surface | 3D curve |
|---|---|---|
| line, parabola | B-spline or Bézier surface, extrusion of a B-spline or Bézier profile | B-spline of degree `k (p + q)`, `k` = 1 or 2, over the part inside the domain |
| B-spline or Bézier curve, polynomial or rational | the same | B-spline of degree `k (p + q)`; a Bézier curve when a Bézier UV curve on a Bézier surface gives one piece |
| circle, ellipse, hyperbola | the same | refused, `NotAnalytic` |

`(p, q)` are the surface's degrees and `k` the UV curve's. A direction in which every UV pole has the same coordinate, exactly, does not count: the net is collapsed at that parameter by de Casteljau first, so a UV curve at constant `v` on a bicubic surface composes to degree `3k`. The result is rational when the UV curve or the surface is. Its parameter is the UV curve's `t`, over the UV curve's bounds, or for an unbounded line or parabola over the interval where it lies inside the domain.

**Circles, ellipses and hyperbolas are refused.** Their parameter is an angle, or a hyperbolic angle, and no rational function of `t` equals `cos t` over an interval. A rational form exists only in another parameter (the half-angle substitution), so it cannot keep the UV curve's parameter, which the method promises.

**The extrusion of a B-spline or Bézier profile is included.** `ExtrusionSurface::to_bspline_patch` clones the profile's spline as the `u` direction, unchanged when the range is the profile's bounds, and adds a degree-1 `v` direction whose knots are the range's ends. Both parameters are kept, which `the_bspline_form_of_an_extruded_spline_keeps_both_parameters` checks against the extrusion at many `(u, v)`, within `1e-14`. The composition reads that patch over the `v` the UV curve's poles cover. A profile that is any other curve has no spline form, and stays `NotAnalytic` off the iso rows.

## Algorithm

This is the B-spline composition of DeRose, Goldman, Hagen and Mann, *Functional composition algorithms via blossoming* (ACM Transactions on Graphics 12(2), 1993), section 4.4: a tensor-product B-spline over a curve is a polynomial between the points where the curve crosses the surface's knot lines.

1. **Crossings.** For each Bézier piece of the UV curve, each distinct knot line `u = s` (or `v = s`) that meets the piece's control hull gives the polynomial `a(t) - s ω(t)` in Bernstein form, `a` and `ω` the homogeneous numerator and weight. Its roots are bracketed with the crate's Bernstein sign-change search (`math::bernstein::rising_roots`) and bisected to the last bit; a line is solved directly. An unbounded line or parabola is first cut to its interval inside the domain, by the roots of a linear or quadratic, solved without cancellation.
2. **Split.** The crossings, merged with the UV curve's own breakpoints (a crossing within a few ulps of a UV knot is the same breakpoint), are inserted as knots, and the UV curve is cut into Bézier pieces.
3. **Patch.** Each piece's midpoint names the surface's knot span in each direction, and the surface's Bernstein net of that patch (`BSplineSurface::patch_net`, the crate's Bézier extraction) is read in homogeneous coordinates `(P w, w)`. The piece's `u(t)` and `v(t)` are mapped affinely into the patch's `[0, 1]²`.
4. **Product.** In the Bernstein basis of the piece,
   `H(t) = Σ_ij C(p, i) C(q, j) U^i (1 - U)^(p - i) V^j (1 - V)^(q - j) (P_ij w_ij, w_ij)`.
   Every product is a Bernstein-basis convolution (`math::bernstein::add_product`), and inside the patch every factor is non-negative, so the sum is a convex combination: no cancellation.
5. **Rational UV curves.** With `u = a / ω`, `v = b / ω`, the factor `B_i^p(a / ω) ω^p` is `C(p, i) a^i (ω - a)^(p - i)`, which keeps the numerators polynomial. The `ω^(p + q)` that multiplies the whole is absorbed by the homogeneous weight, and the degree stays `k (p + q)`.
6. **Join.** The pieces are joined into one B-spline whose interior knots, at the breakpoints, have full multiplicity `n = k (p + q)`. The knots the curve is smoother across are then removed with the crate's knot removal, bounded by `16 ε` times the largest pole coordinate: rounding only. Removal stops at the true continuity. A slanted line across a simple knot line of a bicubic surface, C², keeps 4 of its 6 copies, and across a double one, C¹, keeps 5.
7. **Kind.** The result is a `BezierCurve` when it is one piece and the input was a Bézier curve on a Bézier surface, and a `BSplineCurve` otherwise.

**Periodic directions.** A periodic surface direction has no bounds; its knots are unwrapped across the periods the UV curve covers, and a piece's span is found modulo the period. A slanted line or a parabola that is not held by a bounded direction never ends. It is refused with `NotAnalytic`: a line on a surface periodic in both directions, or on the extrusion of a periodic profile.

**Non-positive weights.** On a rational surface whose weights vary, a rational piece whose control polygon leaves its patch can come out with a non-positive weight coefficient, though the curve itself stays inside. That piece is halved, up to 24 times, until every weight is positive, and refused with `NotAnalytic` after that. With equal surface weights the weight coefficients are those of `ω^(p + q)`, always positive.

## Refusals and the degree rule

- `OutsideSurfaceDomain`, new: a bounded UV curve not inside the surface's domain within `tol.parametric`, or an unbounded line or parabola whose part inside the domain is empty, a single point, or not one piece. A parabola can leave and re-enter the domain; two pieces are two curves, which one `Curve` cannot hold, so the pair is refused rather than one of them chosen.
- `DegreeAboveLimit { degree }`, new: `k (p + q)` above 25, the crate's largest B-spline degree. A cubic UV curve on a biquintic surface gives `DegreeAboveLimit { degree: 30 }`. No other cap is applied: every degree the crate can hold is answered.
- `NotAnalytic`: circles, ellipses and hyperbolas; an unbounded composed interval; a piece whose weights stay non-positive; a non-finite composed pole.
- `Collapsed`: a UV curve whose coordinates are both constant, a single point.

## Considered options

- **The blossom sums of the paper**: rejected. Section 4.4 computes the composed control points as sums of blossom values, which needs one blossom evaluation per multi-index. The Bernstein product form reaches the same coefficients with convolutions of short coefficient lists, and the crate already has them.
- **Sampling the composition and interpolating at degree `k (p + q)`**: rejected. Interpolation in the Bernstein basis solves a Bernstein–Vandermonde system, whose condition number grows exponentially with the degree, so at degrees 12 to 25 the poles lose many digits. Sampling also needs the crossings anyway, since the composition is only piecewise polynomial. The product form has no system to solve, and its points match the surface's to rounding.
- **Approximation or degree reduction**: left out. The result is exact, often at a high degree; a caller that wants a lower one asks for a tolerance. That belongs to a future "approximate within tolerance" operation, where "compose exactly, then reduce degree with a computed bound" is the planned path, so the exact composition is its first step, not a competitor.
- **A cap below the crate's largest degree**: rejected. Any cap would be arbitrary; the crate's `MAX_DEGREE` is the only limit a B-spline has.

## Consequences

`UvCurveTo3DError` gains `OutsideSurfaceDomain` and `DegreeAboveLimit { degree }`. It is unreleased, so this is not a break. Python raises them as `UvCurveTo3DError` with their messages, and adds no API. `BSplineSurface::patch_breaks` exposes the knot spans that `patch_net` numbers its patches by, inside the crate. The composition costs a root search per knot line crossed and a product per piece. `benches/uv_curve_composition.rs` tracks a slanted line on a bicubic surface with interior knots.
