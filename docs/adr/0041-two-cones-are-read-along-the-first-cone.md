# Two cones are read along the first cone

Two cones whose axes are neither parallel with equal half-angles, nor the same line, nor meet where both touch one sphere, and whose apexes differ, meet in a space quartic. Until now it was marched (issue #237). Read along the first cone's generators it has the closed form of ADR 0040: the same walk, discriminants and fit, the ruled quadric a cone.

## Along the first cone's generators

Every generator of the first cone runs from its apex `O₁` along a unit `d`, and meets the second cone, of apex `O₂`, axis `a₂` and half-angle `φ₂`, where `A s² + 2β s + γ = 0`. With `w = O₁ − O₂`:

- `A = (d·a₂)² − cos² φ₂` turns with `d`. It vanishes where a generator of the first cone is parallel to one of the second;
- `β = d·m`, with `m = (w·a₂) a₂ − cos² φ₂ w`;
- `γ = (w·a₂)² − cos² φ₂ |w|²`, the second cone's equation at `O₁`, the same on every generator.

So `Δ = β² − Aγ` is a quadratic form in `d`, `Δ = dᵀ M d` with `M = m mᵀ − γ (a₂ a₂ᵀ − cos² φ₂ I)`. It vanishes along `w`: the line through both apexes meets the second cone at its apex only, twice. `M w = 0`, so `M` has rank two, and across `w` it is one of two kinds:

- `γ > 0`, `O₁` inside the second cone: `M` is positive. Every line through `O₁` crosses the second cone twice, and `Δ` is the sum of two squares, `λ₁ (e₁·d)² + λ₂ (e₂·d)²`: two loops around the first cone.
- `γ < 0`, `O₁` outside: `M` is indefinite, `Δ = (n₀·d)(n₁·d)`, the planes `nₖ·x = 0` through both apexes that touch the second cone. Along the first cone each factor is `aₖ (cos(t − θₖ) − cₖ)`.

The issue expected `Δ` not to factor, and the arc ends to come from a quartic root through `solve_trig`, whose error the fit would carry (its third question). It factors in both pairs, the cylinder being a cone whose apex is at infinity: the arc ends are arccosines, and the walk keeps ADR 0036's exact half-angle factors at every end, with ADR 0040's reference per factor. The arcs are built by the same function for both pairs.

The positive quadratic form of a cone and a cylinder, `cos² φ ((W·a₂)² + A |W|²)` with `A > 0`, is also a sum of two squares, along and across the part of `a₂` across the cylinder's axis. One discriminant, `TwoSquares`, serves both pairs, in six numbers instead of eight.

`A` varies, so the distance along a generator, `(√Δ − β) / A`, and its rate read `A`, `β` and their derivatives from the generator's direction. The cone as the other quadric does so for both pairs; against a cylinder, `A` stays the same.

## What the closed form leaves to the search

`intersect_exact` reports `NotAnalytic`, and `intersect` marches, in four configurations:

- **A generator of one cone parallel to one of the other's**: `d·a₂ = cos φ₁ (a₁·a₂) + sin φ₁ ρ(t)·a₂` ranges over an interval, and `A` keeps one sign only if `|d·a₂|` stays on one side of `cos φ₂`, by more than `tol.angular` in `A`. Otherwise the curve runs off to infinity.
- **Either apex on the other cone**, within confusion. The curve passes through it and crosses itself there.
- **The cones touching**, where a factor vanishes twice: its plane touches the first cone along a generator, within confusion at the point where that generator meets the second cone.
- **A factor that does not turn**, a tangent plane across the first cone's axis. Rare, and marched rather than given its own case.

Near them the closed form holds. With the second apex `2 × confusion` off the first cone, the two loops pass `2 × confusion` apart where they would cross, and are fitted within `0.3 × confusion` of their closed form. The march, within confusion of both, traces them as one.

## Unbounded branches

Two cones in general position often run off to infinity: their direction cones overlap, as for any two wide cones turned apart. The exact curves that run to infinity already in the library, the lines of a plane through a cone's apex, the hyperbolas of a plane or of two parallel cones, are returned whole, with an infinite `parameter_range`, exact and complete. A trimmed surface narrows that range to the patch (ADR 0039); it is not what makes the curve exact. Unbounded quartic branches will follow the same rule, in their own change: the curve, not a piece of it cut at a box. What a quartic lacks is a representation: `intersect` returns its fit (ADR 0036), and a fit has bounds. ADR 0040's reason for leaving them to the march, that a trimmed surface bounds them, is withdrawn.

## Considered options

- **Along the second cone's generators**: equivalent. The first operand is ruled, so each order has its own closed form; both meet in the same loops.
- **`Δ` as a trigonometric polynomial with `solve_trig` roots**: rejected. It factors exactly, and the roots would carry the solver's error into the fit near the arc ends.
- **Unbounded branches in this change**: deferred, as the issue suggested. They need a component that runs to infinity in the walk and a representation in `intersect`, for a cone and a cylinder too.
