# An unbounded quartic branch is an exact curve

Two cones in general position often meet in a curve that runs off to infinity: whenever a generator of the first cone is parallel to a generator of the second. So does a cone and a cylinder whose axis runs along a generator of the cone. ADR 0041 left those to the march and decided that, like the lines and hyperbolas already returned whole, they would come back as the curve, not as a piece cut at a box. A fit has bounds (ADR 0036), so they needed a representation. This is it: `Curve::QuarticBranch`, an exact curve with an infinite range, for two cones first (issue #237).

## Where the walk runs to infinity

Along a generator `d` of the first cone, `A = (d·a₂)² − cos² φ₂` vanishes where `d·a₂ = ± cos φ₂`. With `d·a₂ = P + R cos(t − θ)` along the first cone, those generators are arccosines, at most four. There `Δ = β² − Aγ = β²`, so of the two roots `s = (±√Δ − β) / A` one stays finite, `−γ / 2β`, and the other runs off on the side `−sign β`. Each such generator is a point at infinity on one component of the walk: on a loop around the first cone, the one of that branch; on an arc, the parameter where the walk passes that generator on that root. A component closes through its points at infinity, as a projective curve. In space it falls apart into the pieces between them, each running to infinity at both ends.

The pieces are `QuarticCurve::branches`. A component with no point at infinity is still fitted, exact and complete within confusion, with its curves in parameters. Each piece of one that has them is an `IntersectionCurve` holding a `QuarticBranch`, with the parameter range `(−∞, ∞)` and no curves in parameters.

The finite root is read as `−γ / (√Δ + β)` where `√Δ` and `β` share a sign: near a point at infinity `A → 0`, and `(√Δ − β) / A` would cancel to rounding.

## A cylinder along a generator of the cone

Read along the cylinder's generators, `A = (a₁·a₂)² − cos² φ` is the same on every one, and vanishes when the axis `a₁` runs along a generator of the cone. Within `tol.angular` of that, the quadratic is linear: each generator of the cylinder meets the cone once, at `s = −γ / 2β`, and the curve is one loop around the cylinder (`Discriminant::MetOnce`). The point is read from that quotient directly, not through a root of `Δ = β²`, so that no rounding of the signs near `β = 0` can pick the root that runs off.

With `W(t)` from the apex to the generator's base, `β = n·W(t)` with `n = (a₁·a₂) a₂ − cos² φ a₁`, across the axis: `k₀ + k₁ cos t + k₂ sin t`, vanishing at two arccosines or none. Where it vanishes the point runs off to infinity along that generator, and the loop falls apart into two branches. Where it never does, the cylinder is out of reach of the cone's generator along its axis and the loop is bounded, fitted like any other. Where `β` touches zero, within `tol.angular`, the point runs off and comes back on the same side: one branch, growing as `u²` rather than `|u|`. `γ` and `β` vanish together only on a generator of the cylinder that lies on the cone, through the apex: the apex on the cylinder, still left to the search. *Amended (ADR 0045):* read as that generator and a cubic, in one branch.

The true `A` within `tol.angular` of zero puts the second root about `2β / A` away, and the point read as met once is off the cone by about `A s`: nothing at the distances a model reaches.

`QuarticBranch::from_cone_and_cylinder(cone, cylinder, index, tol)` rebuilds such a branch. `QuarticBranch::surfaces` gives the two surfaces back, two cones or a cone and a cylinder, and `QuarticBranch::tolerance` the tolerance it was read with (see below).

## The parameter

The walk parameter `w` reaches a point at infinity at a finite value, the piece lying between two of them, `w₀ < w₁`. `A` has a simple zero there, so the distance along the generator grows as `1 / (w − w₀)`. The branch's parameter is `u = tan((w − m) π / 2h)`, `m` the middle of `(w₀, w₁)` and `h` its half width: `u` runs over the whole real line, and the distance grows about as `|u|`, as on a line or a parabola. A trimmed branch over a range of `u` is a piece of about that size.

## What is exact and what is measured

Every point is exact, a root on a generator of the first cone, mapped by the branch's placement. The rest has no closed form:

- **Derivatives**: central differences of exact points in `u`, five points for the first two orders, with steps relative to `1 + |u|` near the balance of truncation and rounding: within about `1e-9` relative.
- **Projection**: the least distances among 64 points spread over the whole branch (`u = tan(π (q − ½))`, `q` in `(0, 1)`), refined by Brent's search between neighbours and polished by Newton steps. The distance grows without bound at both ends, so the closest point is among them.
- **Box, discretization, length**: measured on exact points, where every other curve has a proof (ADR 0003). The box of a range goes through the coordinates' extremes found among 64 samples, refined by Brent's search; the discretization takes the fewest equal steps whose measured deviations meet the criteria; the length integrates the speed. As the fit of ADR 0036 is measured against the closed form, these are measured against the exact points.

Far out the points lose precision: `A` cancels to about `1e-16` in its terms, and the distance carries that error relatively, about `1e-16 u²` times the cones' size. That is within the default confusion up to `|u|` of about `1e4`.

## A curve like the others

`QuarticBranch` is a `Curve` variant: it evaluates, projects, reverses (`u` negated), transforms (its placement composed with the transform), and intersects through the numeric search, as a B-spline does. It has no bounding box and no B-spline (`Unbounded`), and no curve in any surface's parameters (`NotAnalytic`). A swept surface takes it as a profile: an extrusion projects onto it by searching along it, a revolution by searching out to where it passes the reach.

`QuarticBranch::from_cones(first, second, index, tol)` builds one from the cones it lies on, its place among the quartic's branches and the tolerance it is read with; `from_cone_and_cylinder` does the same for a cone and a cylinder. With `surfaces`, `tolerance`, `is_reversed` and `placement`, that is what the Python class pickles.

*Amended (issue #283):* the curve keeps the tolerance it was read with, and a branch is rebuilt with it. Within a tolerance the walk reads a factor as touching zero, an apex as on the other cone, `A` as vanishing twice and `β` as vanishing with it, and those readings are not monotone in the tolerance: a rebuild with zero tolerances, or a looser one, can find another curve, and a branch that came from an intersection would not come back. The tolerance adds 24 bytes to `QuarticCurve`, which every exact result carries by value: about 20 to 40 instructions more for an exact quartic, nothing measurable where it is fitted. `QuarticBranch::LOOSEST_ANGULAR`, which the rebuild from a cone and a cylinder read within, is gone.

## What is still left to the search

- **`A` vanishing twice**, a generator of the first cone touching the second's directions: `P ± R` within `tol.angular` of `± cos φ₂`. *Amended:* read as vanishing once, on that generator. `A` keeps its sign on either side of it, so the root that runs off is the same one from both sides, on the side `−β`: the curve runs to infinity there once, and the loop through it is one branch, both its ends out along that generator. Where `A` only comes within `tol.angular` of zero, the curve's points near that generator lie out beyond about `|β| / tol.angular` along it, and are read as at infinity. Where `β` vanishes on that generator too, the next case.
- **`β` vanishing at a point at infinity** too, within confusion: both roots run off together, along a branch like a parabola's. *Amended:* read on the walk. With `A` and `β` vanishing on a generator, so does `Δ = β² − A γ`, and it changes sign there with `A`, as `−A γ`: the generator ends an arc, and both roots run off as it nears, `±√Δ / A`, one each way along it. That end is the curve's one point at infinity there, at the arc's walk parameter `0` or `π`, the branches meeting at it. The generator is matched to the arc's end within `1e-6` in angle: they are `β² / |A′ γ|` apart with `β` within confusion, and the end's arccosine is good to about `1e-8` near `± 1`. The other way round, along the second cone, the same point is one root running off as usual. Where `A` also vanishes twice there, the plane through the generators along the same direction touches both cones, and the pair stays with the search.
- Each apex on the other cone, the line through both on both. (The cones touching: ADR 0043; one apex on the other cone: ADR 0044.) *Amended (ADR 0045):* read as the line and a cubic, whose pieces between its points at infinity are branches like these.

## Considered options

- **A piece cut at a box**: rejected, as ADR 0041 decided. The curve is exact and unbounded, as a cone's lines and hyperbolas are, and a trimmed surface narrows it to a patch (ADR 0039).
- **A fit over a long finite range**: rejected. It has bounds, and the range would be a guess at what the caller needs.
- **The walk parameter as the curve's parameter**: rejected. Its range is a finite open interval whose ends are points at infinity, so the bounds would be where the curve is not defined, and equal steps in it crowd the points near the ends.
- **Exact derivatives**: deferred. They follow through the walk, the root and `A` by the chain rule; central differences are within `1e-9` relative, and nothing yet needs more. Profile first.
- **Transforming the cones instead of keeping a placement**: rejected. A mirror turns a cone's frame, and with it the angles of the walk; the placement maps the exact points and keeps the walk as it was found.
- **A variant for the whole projective loop**: rejected. A loop through infinity is not a curve in space; its pieces are.
