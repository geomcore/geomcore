# Evaluation past the bounds clamps a spline

Every `point_at` and `derivative` evaluates any finite parameter, inside the shape's bounds or not. What it gives past them depends on the shape (issue #179):

- A line, a conic, and the plane, cylinder, cone, sphere and torus continue along their own equation, which is defined everywhere. A periodic direction wraps.
- A trimmed curve evaluates its basis curve, so a trimmed line goes on along the line and a trimmed spline gives what the spline gives.
- An extrusion or revolution continues along its sweep, with its profile evaluated as that curve is past its own bounds.
- A non-periodic direction of a B-spline or Bézier curve or surface **clamps** the parameter to the nearer bound. The point and every derivative are those at the bound, as if the parameter were the bound.

Operations that build a new shape from a parameter refuse one outside the bounds: `split_at`, `extract`, `insert_knot` and the iso-curves `u_iso` and `v_iso` of a spline surface return `OutOfRange` (or `KnotOutOfRange`). They already return a `Result`, and a curve built at a clamped parameter would be degenerate (a split at the end leaves an empty piece) or would silently stand for another place.

So the inconsistency reported in #179, `point_at(99, 99)` answering where `u_iso(99)` raises, is two contracts on purpose. Evaluation is total, and construction is checked.

## Why clamping

Evaluation is infallible and sits on every numeric hot path: projection, discretization, the curve–surface and surface–surface searches, marching. Those routines step and sample near the bounds and rely on a sample past the bound staying on the shape. Extrapolating the last polynomial span was tried: the curve–surface search then found hits on the continuation outside the patch where there are none (`curve_surface_hull::tests::test_misses_have_no_hits`), and the marching ran off the patch on a loop that touches its edge (`marching::tests::test_marching_small_loop`). Clamping keeps every sample on the shape the caller holds.

## Considered options

- **Extrapolate along the last span's polynomial**: rejected. It is what a trimmed elementary curve does, but a spline's continuation has no meaning beyond its knots, and the searches above would need their own clamping to stay on the patch.
- **Return `Result` or `Option` from evaluation**: rejected. It costs a branch and an unwrap in every inner loop, and every caller that evaluates inside the bounds, which is nearly all of them, would carry an error it never sees.
- **Return NaN past the bounds**: rejected. It poisons sums silently, and the searches that sample a hair past a bound by rounding would lose their samples.
- **Clamp in construction too, so `u_iso(99)` gives the boundary iso-curve**: rejected. A caller asking for an iso-curve or a split at a parameter outside the shape has a wrong parameter, and building at the bound would hide it. Construction already reports the error.
- **Zero derivatives past a clamped bound**: rejected. The clamped point does not move, so its derivative is zero in the strict sense, but every caller that reads a tangent at a parameter a hair past the end by rounding wants the end tangent. Point and derivatives are evaluated at the clamped parameter together, so they stay consistent with each other.
