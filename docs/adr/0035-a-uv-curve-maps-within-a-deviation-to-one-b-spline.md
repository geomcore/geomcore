# A UV curve maps within a deviation to one B-spline

`uv_curve_to_3d` gives the exact 3D curve of a UV curve or refuses the pair (ADR 0030). The pairs it refuses include the commonest one in exchanged models: a slanted line on a cylinder or cone, whose 3D curve is a helix or a conical spiral. Those curves are transcendental, so no rational curve matches any arc of them. A circle or ellipse in the `(u, v)` plane of a cylinder, sphere or torus composes to `cos(c + r cos t)`, which has no exact B-spline either. A B-rep still needs a 3D curve for each such edge, so a second operation gives one within a deviation the caller chooses:

```rust
pub fn uv_curve_to_3d_within<S, C>(surface: &S, uv_curve: &C, max_deviation: f64, tol: Tolerance)
    -> Result<CurveApproximation, UvCurveTo3DApproximationError>
where S: ParametricSurface + ?Sized, C: ParametricCurve2D + ?Sized
```

It is a generic free function, as its mirror `curve_3d_to_uv_within` is. It works over the UV curve's own bounds (ADR 0005), so a line is trimmed first. The two operations stay separate. A caller that wants the exact curve where there is one calls `uv_curve_to_3d` first and this one on `NotAnalytic`, which keeps the choice visible (ADR 0030).

## The result

The result is **one curve for one UV curve**: a non-rational, non-periodic B-spline in the UV curve's parameter, so its point at `t` is near `S(c(t))`. The edge it belongs to keeps one parameter for its 3D curve and its UV curve. The deviation is the same-parameter one, `max_t |Q(t) − S(c(t))|`, which is what a B-rep's edge tolerance means. A Hausdorff deviation would be weaker, because it lets the two parametrizations drift apart.

Its degree is one of two:

- **Degree 9, C⁴**, where `f = S∘c` is analytic: on an elementary surface along a UV curve with no knot, which is a line, a conic, a Bézier curve or a spline with no interior knot. A straight spline counts as a line.
- **Degree 5, C²** wherever the UV curve is, everywhere else. At a knot where a spline UV curve is only C¹ or C⁰, the 3D curve keeps the same continuity, through a knot of multiplicity 4 or 5.

Degree 5 is also kept where it takes no more poles than degree 9, as on a short arc that one span covers at either degree. The degree is read from the curve.

`CurveApproximation` (ADR 0012) gives the curve, its deviation, and whether that deviation is proven. `max_deviation` is a contract. When it cannot be met the call is refused with `DeviationNotReached { reached }`, never answered with a looser curve. The reported deviation is the bound or the measured value; it is never replaced by the asked one.

## The algorithm

1. **Spans of Hermite interpolants.** On each span `[a, b]` the curve is the polynomial of degree `p = 2n + 1` that matches `f`, and its derivatives up to order `n`, at both ends: `n = 2` at degree 5, `n = 4` at degree 9. The derivatives come from composing Taylor series. With `δu(s)` and `δv(s)` the UV curve's series from `t`, `f(t + s)` is the sum of `S_ij δuⁱ δvʲ / (i! j!)`, so the order-`n` jet of `f` needs the surface's partials up to total order `n` and the UV curve's derivatives up to order `n`. At order 2 this is the chain rule, `f'' = S_uu u'² + 2 S_uv u'v' + S_vv v'² + S_u u'' + S_v v''`. The span's Bézier poles follow in closed form: pole `k ≤ n` is `Σⱼ C(k, j) (p − j)! / p! hʲ f⁽ʲ⁾(a)`, and pole `p − k` is the same from `b` with `−h`. Neighbouring spans share their end jets, so the joints are `Cⁿ` by construction.

   Every curve evaluates derivatives of any order. The plane, cylinder, cone, sphere and torus do too, from `cos(u + iπ/2)` and `sin(u + iπ/2)`. The other surfaces evaluate up to order 2, so they stay at degree 5. A spline UV curve stays at degree 5 as well: its knots are breakpoints of `f`, and at degree 9 a C² knot would need multiplicity 7, more poles than degree 5 gives.

2. **Proven span count for a straight UV curve on an elementary surface.** For a two-point Hermite interpolant of degree `2n + 1`, the error at each coordinate is `f⁽ᵐ⁾(ξ) ((t − a)(t − b))ⁿ⁺¹ / m!` with `m = 2n + 2`. On a span of length `h` it is at most `max |f⁽ᵐ⁾| (h / 2)ᵐ / m!`. Along `u = u₀ + a t`, `v = v₀ + c t` every coordinate of the surface in its frame is a product of cosines, so `f⁽ᵐ⁾` has a closed-form bound:
   - cylinder: `r aᵐ`
   - cone: `max|ρ| aᵐ + m |c sin α| aᵐ⁻¹`, by Leibniz, the radius `ρ` being linear in `v`
   - sphere: `R ((a + c)ᵐ + (a − c)ᵐ) / 2` across, and `R cᵐ` along the axis
   - torus: `R aᵐ + r ((a + c)ᵐ + (a − c)ᵐ) / 2` across, and `r cᵐ` along the axis
   - plane: `0`

   The bounds of the frame's coordinates combine into the norm, which does not depend on the frame. The span count at each degree is then the fewest equal spans within the budget: computed, not searched, and the deviation is proven up to rounding. The degree with fewer poles is kept, `(n + 1)(N + 1)` for `N` spans, and degree 5 on a tie. For a helix it is the bound of Lu (2012), *On polynomial approximation of circular arcs and helices* (Computers and Mathematics with Applications 63), applied to the curve's own parametrization. On a radius-10 cylinder at `1e-7`, a turn takes 5 spans at degree 9, against 24 at degree 5.

3. **Sampled spans for every other pair.** The search starts from the UV curve's knots, where the jets are one-sided, and from at least four spans, so that a closed UV curve cannot fit at its samples by chance. Each span is checked at nine evenly spaced interior points, its midpoint included, where the error factor `((t − a)(t − b))ⁿ⁺¹` is widest. A span is accepted when that sampled deviation is within half the budget. Otherwise it asks for the number of pieces its error predicts, `1.1 (deviation / accepted)^(1/m)`, a real number up to 64. Each run of consecutive failing spans that no knot of the UV curve breaks gets the sum it asks for, rounded up. The new breakpoints are spread so that each piece holds an equal share of what was asked. Four spans that each ask for 1.2 pieces become five, not eight: at degree 9, halving a span divides its error by 1024, so splitting each span on its own would overshoot. This is a measurement, so the result says it is not proven. Kinks along a surface's knot lines are not breakpoints; the splitting finds them, since the error near one still falls with the span's length.

   An analytic pair is sampled at degree 9 first. Degree 5 is then tried with a span limit, so that it is kept only when it takes no more poles: `3(N₅ + 1) ≤ 5(N₉ + 1)`. A limit it passes stops the trial within its first pass.

4. **Join.** The spans become one B-spline directly, with no knot removal. Each pole of a degree-`p` B-spline is the blossom of the polynomial on a span it bears on, at the `p` knots after it. Knots over two neighbouring breakpoints give a span's Bézier pole, and `p` copies of one breakpoint give its point. Knots over three breakpoints, `α` copies of the one below, a whole joint `t_i` and `β` copies of the one above, only occur around a joint of multiplicity `n + 1`, whose spans share one jet. Expanded from `t_i`, the blossom of `(t − t_i)ᵏ` is `e_k / C(p, k)`. Here `e_k` is the elementary symmetric polynomial of the knots' offsets from `t_i`: `−h_l` `α` times and `h_r` `β` times, `h_l` and `h_r` being the lengths of the spans on either side. It is zero past `k = α + β ≤ n`, so the shared jet is enough: the pole is `Σₖ f⁽ᵏ⁾(t_i) e_k / (k! C(p, k))`. At degree 5 around a knot of multiplicity 3, it is `f + (h_r − h_l) f' / 5 − h_l h_r f'' / 20`. A knot where a spline UV curve is only C¹ or C⁰ has a jet on each side. It gets multiplicity 4 or 5, read from the UV curve's own multiplicity, and its poles are all Bézier poles. The join is exact up to rounding, so the deviation is the spans' own. A first version built the curve at full multiplicity and removed knots within a small share of the budget. Its cost grew quadratically with the span count: 238 µs instead of 36 µs for ten turns of a helix at degree 5.

**Refusals.**

- `InvalidDeviation`: a `max_deviation` that is not positive and finite.
- `Unbounded`: a UV curve with an infinite bound.
- `Collapsed`: a 3D curve that is one point, such as a line along a sphere's pole or a straight spline with coincident poles. It is refused, not approximated, as in ADR 0030.
- `OutsideSurfaceDomain`: a UV curve that leaves a non-periodic direction of the domain by more than `tol.parametric` at any point evaluated. For a straight UV curve the breakpoints include both ends, so the check is exact.
- `DeviationNotReached`: more than 100 000 spans, or a span shorter than `tol.parametric` that still fails.

## Considered options

- **A helix curve type.** Rejected for now. A kernel can carry an exact helix, but it would cover one row of the table, the cylinder or cone with a straight UV curve. Every consumer (intersection, projection, exchange) would need a new variant, and exchange formats have no such entity, so an edge would still need a B-spline to be written out. It can come later beside this operation, which every pair needs anyway.
- **Degree 5 everywhere.** The first version of this ADR. Rejected: over ten turns of a helix at `1e-7` it takes 723 poles, against 235 at degree 9 and 162 from a production kernel fitting degree 11 to 14. Its spans grow only as the sixth root of the budget, against the tenth root at degree 9. Degree 9 also builds the curve faster, 18 µs against 34 µs, with five times fewer spans.
- **A higher degree, 11 or 13.** Rejected: each order of derivatives saves fewer poles than the one before. From degree 9 to 13, the poles per turn fall from 25 to 14 at `1e-7`, while the cost of evaluating the curve doubles, the surfaces must give sixth-order partials, and the poles are worse conditioned. Degree 9 stays within 1.5 times a production kernel's pole count, 150 against 98 over ten turns at `1e-5`, and its points evaluate with their first derivative in about 100 ns, against 200 to 260 ns at that kernel's degrees.
- **A degree chosen by the caller.** Rejected: the degree follows from the pair, and the fewest poles is the only reason to choose. A caller that needs a lower degree can reduce the curve.
- **Projection onto orthogonal polynomials at Gauss points, then conversion to a B-spline.** This is the usual approach of production kernels. Rejected as the core: it measures its error at the Gauss points only, and it has no closed-form span count to prove. The Hermite spans get the same order with fewer evaluations, and the joints come out `Cⁿ` without a constrained fit.
- **Chebyshev interpolation per span with coefficient chopping** (Aurentz and Trefethen 2017). The error estimate comes from the series tail, which is again a measurement. The conversion to Bernstein form adds a step for no gain.
- **Composing an approximated UV curve exactly on the surface** (Yang et al. 2008, *Approximate computation of curves on B-spline surfaces*). That gives a curve exactly on the surface, but in a different parameter, and it needs an exact rational composition that the elementary surfaces only have through their rational forms, at degree `p + q`. Same-parameter is the deviation an edge needs.
- **Several curves for one UV curve.** Rejected: an edge has one 3D curve. The spans are pieces of one B-spline, not separate curves.
- **A lazy curve-on-surface type** that evaluates `S(c(t))` and caches an approximation. Deferred, until a consumer needs evaluation without a stored curve. It would hold this operation's result as its cache.
- **Approximating inside `uv_curve_to_3d`.** Rejected (ADR 0030): one name, one promise.
- **Copying one span of a cylinder's helix**, since its spans are congruent. Deferred, pending profiling (YAGNI).

## Consequences

`derivative` on the plane, cylinder, cone, sphere and torus evaluates every order; it panicked past second order. `ElementarySurface` gives all partials up to an order at once, from one sine and cosine per angle, for the order-4 jets.

A new hidden routing hook, `ParametricCurve2D::as_line`, tells a straight UV curve from any other through a generic `C`. `Line2D` returns itself; `Curve2D`, `&C` and `TrimmedCurve2D` forward it (ADR 0022). `uv_curve_to_3d` reads straight UV curves through the same detection, `StraightUvCurve::of`, now generic.

A crate surface gets proven spans and degree 9 as soon as `as_elementary` returns it. Any other surface, including a user's, gets sampled spans at degree 5 with no further code. An extrusion or revolution with an elementary profile is analytic too, but stays at degree 5 until its derivatives go past second order and a consumer asks for it.

The Python bindings expose it as the surface method `surface.uv_curve_to_3d_within(uv_curve, max_deviation, tol)`, and its mirror as `surface.curve_3d_to_uv_within(curve, max_deviation, tol)`.

Amended (#185): in Python both operations were first module functions, `geomcore.uv_curve_to_3d_within(surface, …)`, mirroring Rust. Users looked for them on the surface, beside `uv_curve_to_3d` and `curve_3d_to_uv`, in IDE completion. The reason for a free function in Rust is genericity over any `ParametricSurface`, including a caller's own. That reason does not hold in Python, where every surface is a concrete class. The methods replace the functions, which warn with a `DeprecationWarning` for one release. The Rust free functions stay.

## Amendment: the mirror, 3D curve to UV curve (#186)

`curve_3d_to_uv_within(surface, curve, max_deviation, tol) -> Result<CurveApproximation2D, Curve3DToUvApproximationError>` replaces the sampled `curve_3d_to_uv_numeric`, which returned raw `(u, v)` samples that wrapped at the seam, for a sample count instead of a deviation. It shares this operation's spans through a `SpanSource` trait: the jets come from inverting the 3D curve's points, and the UV curve's derivatives from the surface's Jacobian, `c' = (JᵀJ)⁻¹ Jᵀ C'` and the same for `C''` less the surface's curvature terms. The UV coordinates ride in the `x` and `y` of the 3D jets, so the Hermite poles and the join are unchanged. A span's deviation is measured in 3D, `|S(c(t)) − C(t)|`, the same-parameter deviation.

- **Degree 5 only.** Order-4 jets would need third and fourth derivatives of the inverse, and the deviation is always sampled: no closed-form bound exists for an inversion.
- **The seam.** Each inversion is moved by whole periods to the value nearest a lift: 129 samples of the UV curve along the curve, each predicted from the one before and its derivative. The UV curve is continuous, and a closed curve around a cylinder ends one period from its start.
- **Refusals.** As here, plus `CurveNotOnSurface { distance }` when an inversion is farther than half the deviation from the curve, and `SingularParametrization` when a breakpoint falls where `JᵀJ` is singular, as a curve that starts or ends at a sphere's pole (amended below, #220). A curve through a pole between its breakpoints is still within the deviation: every `u` is the same point there, so the UV curve turns steeply in `u`. Splitting at the pole, which a B-rep does for its vertex, is the caller's.

## Amendment: a curve to a pole takes the limit from inside (#220)

A B-rep edge that ends at a sphere's pole or a cone's apex is a curve whose breakpoint falls where `JᵀJ` is singular, and its pcurve is what a kernel builder needs, so refusing it with `SingularParametrization` left the case to the caller. At such a breakpoint the UV curve's jet is now its limit from the side of the span it bounds: the point and its first two derivatives are read a short step inside, `1e-4` of the curve's range and at most half the way to the curve's next knot, and carried back along their Taylor series. The error is of the order of the step's cube; the inversion's own error does not grow as the step shrinks, as the step's powers balance its conditioning. A meridian keeps its own `u` up to the pole, and a small circle through a pole leaves it at the angle of its tangent there.

- **Singular means at the point.** The parametrization is taken as singular where `det(JᵀJ)` is within `1e-24` of the square of its largest diagonal term, a point within `1e-12 r` of a sphere's pole, instead of `1e-12`, within `1e-6 r`. A breakpoint just beside the pole, which the span search lands on as it shrinks the spans around a curve that crosses it, was read as a limit from its lower side, across the pole, and got the other side's `u`. Nearer than `1e-6 r` the derivatives lose precision as the inverse of the distance and its square, which the span's length and the distance multiply away in the surface's point.
- **The lift.** A lift sample at a singular point has any `u` and predicts nothing, so a point after it is moved by periods toward the next sample instead.
- **What is still refused.** `SingularParametrization` remains for a curve that lingers at the point, nearing it as the fourth power of its parameter or slower, where the step inside is still singular. A curve through a pole at an interior breakpoint takes the limit from below on both sides, and the spans above turn steeply in `u`, as between breakpoints.

## Amendment: one result for every approximation (#199)

The two operations had a result type each, `UvCurveApproximation` and `Curve3DToUvApproximation`, and the B-spline fits within a deviation returned the bare spline, so a caller could not tell how close a fit came. They now share one result per tier, `GenericCurveApproximation<P>` with the aliases `CurveApproximation` and `CurveApproximation2D` (ADR 0014), and `SurfaceApproximation` for `BSplineSurface::approximate_within`. Each holds the spline, its deviation and whether that deviation is proven. What the curve is and how its deviation is measured stays in each operation's documentation: the result type is the same, the contract is the operation's.

- **Errors stay per operation.** Their variants differ (a fit fails to construct, a mapping leaves the domain or is not on the surface), so one error type would carry variants most calls cannot return. `UvCurveApproximationError` becomes `UvCurveTo3DApproximationError`, which names the operation as its mirror's error does.
- **A fit is proven.** A fit is measured at every point it was given, which is the whole of what it promises, so its deviation is exact and `is_proven()` holds. `curve_3d_to_uv_within` is always sampled; `uv_curve_to_3d_within` is proven where the closed-form bound above applies.
