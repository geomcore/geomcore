# A cone and a cylinder are read along the cylinder

A cone and a cylinder whose axes are neither the same line nor meet where both touch one sphere meet in a space quartic. Until now it was marched (issue #237). It has a closed form along the generators of either surface, as the quartics of ADR 0036 do. Read along the cylinder's, it is the same walk, components and fit.

## Along the cylinder's generators

The cylinder's generator at angle `t`, `B(t) + s a₁`, meets the double cone of apex `O`, axis `a₂` and half-angle `φ` where `A s² + 2β s + γ = 0`, with `A = (a₁·a₂)² − cos² φ` the same on every generator. A cone's own generators change direction with `t`, so read along them `A` would vary and vanish where a generator of the cone is parallel to the cylinder's axis. Along the cylinder's, it does not.

`Δ` does not change when the base of the generator slides along it. Take the base level with the apex: `W(t)`, from the apex to it, runs across the cylinder's axis, around a circle of the cylinder's radius `r`. Then `β = (a₁·a₂)(W·a₂)` and

`Δ = cos² φ ((W·a₂)² + A |W|²)`.

- `A > 0`: the cylinder's axis lies within the cone's opening, and so does every generator. `Δ` is a sum of positive terms, so every generator crosses the cone twice: two loops around the cylinder.
- `A < 0`: `Δ` is a quadratic form in `W`, positive in one direction and negative across it. On the plane across the cylinder's axis it vanishes on two lines through the apex, the edges of the cone's shadow, and `Δ = cos² φ (n₀·W)(n₁·W)`, with `nₖ` normal to each edge and `|nₖ| = √(1 − (a₁·a₂)²)`. As `W` runs around its circle each factor is `r |nₖ| (cos(t − θₖ) − cₖ)`. The generators that meet the cone are those whose base lies in the shadow, and the arcs of them end where the circle crosses an edge: `t = θₖ ± acos cₖ`.

The issue expected `Δ` not to factor, and the arc ends to come from a quartic root through `solve_trig`. It does factor for `A < 0`, over the reals. So the arc ends are arccosines, and the walk keeps ADR 0036's exact half-angle factors at every end: a factor vanishing at an end is evaluated from the parameter's exact distance to it, not from `cos(t − θₖ) − cₖ`. The walk changes in two places. Each factor has its own reference angle `θₖ`: the arcs of a pipe around the apex run from one edge of the shadow to the other. And with `A < 0` the distance along a generator, `(√Δ − β) / A`, falls as `√Δ` grows, so an arc leaves the edge against the generator's direction; the fit's tangent takes its sign from that rate, not from the direction. For `A > 0` the sum of squares is evaluated as such, and it never cancels.

## What the closed form leaves to the search

`intersect_exact` reports `NotAnalytic`, and `intersect` marches, in three configurations:

- **A generator of the cylinder parallel to one of the cone's**, `A` within `tol.angular` of zero. The curve runs off to infinity, as the march's box shows (`IncompleteReason::Unbounded`). A curve without bounds cannot be fitted. Trimming the surfaces to the patches the caller needs (ADR 0039) bounds the march, and the curve with it. This settles the issue's first question for this pair: no closed form for unbounded branches, because a trimmed surface already gives them bounds.
- **The apex on the cylinder**, within confusion. The curve is a rational quartic with a node at the apex.
- **The surfaces touching**, where an edge of the shadow touches the circle within confusion. The curve crosses itself, or is a point.

These are the issue's singular cases. Each needs its own parametrization, as the cylinders' node did (#234). They stay with the march, which handled them before, until a model needs them.

Near them the closed form holds. A cone around a sphere `2 × confusion` wider than the cylinder's is fitted within `0.3 × confusion` of its curve, which nearly splits into the reducer's two ellipses. ADR 0036's distance to the curve, not to the surfaces, keeps it there.

## Considered options

- **Along the cone's generators**: rejected. `A` varies with the generator, the loops of the arcs need the zeros of a varying leading coefficient, and the walk would need a second kind of root. Two cones will need it; a cone and a cylinder do not.
- **`Δ` as a trigonometric polynomial with `solve_trig` roots**: rejected. The arc ends would carry the solver's error, and the fit near them would read `√Δ` from a value that cancels.
- **Close the unbounded branches at a box**: rejected. The march already does so and reports `Incomplete`. A trimmed surface bounds them for a caller who needs them.
