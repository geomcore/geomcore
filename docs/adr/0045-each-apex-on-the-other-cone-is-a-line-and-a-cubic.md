# Each apex on the other cone is a line and a cubic

Where each cone's apex lies on the other, the line through both apexes is a generator of each. The quartic splits into that line and a cubic. ADR 0044 left the pair to the march (issue #237). It now has a closed form, read along the first cone's generators.

## Reading the apexes on the cones

Each apex within confusion of the other cone is not enough. Apexes a few confusions apart are each within confusion of the other cone, yet the line through them points anywhere. The line must also lie along a generator of each within `tol.angular`: each apex's offset from the other cone is at most `tol.angular` times the distance between the apexes. The rest of the walk then holds to about that angle, so a point at distance `R` is off the cones by about `R tol.angular`. Otherwise the pair goes to the search, as before. An apex on a cone built through it is off by rounding, about `1e-16` of the distance, and passes.

## The cubic on the walk

The first apex lies on the second cone, so `γ` vanishes on every generator of the first, and the quadratic `A s² + 2β s` has the root `s = 0`, the apex, and `s = −2β / A`. On the shared generator, at `t₀`, `A` and `β` both vanish, and the ratio is nought over nought.

The common factor is cancelled exactly. Write `Q = a₂ a₂ᵀ − cos² φ₂ I`, so that `A = Dᵀ Q D` and `β = Dᵀ Q w`, with `w` from the second apex to the first. A quadratic form gives `q(v) − q(u) = (v − u)ᵀ Q (v + u)`, and on a cone `D(t) − D(t₀) = 2 sin δ D'(t₀ + δ)`, with `δ` half the angle between the generators. Both coefficients vanish on the shared generator, so each is its difference from there. That divides both by `2 sin δ`, leaving

`s = −2 (D'·Q w) / (D'·Q (D + D₀))`, with `D'` taken midway.

This is exact through `t₀`, where the cubic passes through the second apex. It crosses the line again at the first apex, where `D'·Q w` vanishes.

The cubic runs to infinity where `A = 2 sin δ (D'·a₂)(a₂·(D + D₀))` vanishes away from `t₀`:

- **`D'·a₂ = 0`**: once, `t = 2h − t₀` with `h` the angle of `a₂` across the first axis. This is the other zero of the factor `a₂·D − a₂·D₀`.
- **`a₂·D = −a₂·D₀`**: twice, or never. Where it only comes within `tol.angular` of that value, it is once, and the branch goes out and comes back along one generator on the same side, as for `A` vanishing twice (ADR 0042).

The two factors differ by the constant `2 a₂·D₀`, which is never zero, so they never vanish together. A real cubic has at least one point at infinity, so the cubic always comes as one to three `QuarticBranch` pieces, never fitted. Each branch is its pieces between consecutive points at infinity, as for the quartic.

The cubic is one component of a `QuarticCurve`, `OffSharedGenerator { shared: t₀ }`, under the discriminant `MetOnceOffTheApex`. Its parameter is the generator's angle, as for a loop around. Through the second apex the second cone's normal vanishes, so tangents are central differences of exact points, and the distance to the curve is measured by sliding along it (ADR 0044).

## The result

`ConeConeIntersection::LineAndCubic(Line, QuarticCurve)`: the line from the first apex toward the second, and the cubic. `intersect` returns the line and the cubic's branches, exact and complete. The two apexes are its contact points, where the curves cross, as for the crossing conics of two cones about one sphere. `QuarticBranch::from_cones` rebuilds a branch of the cubic, as of a quartic.

Cones that touch along the shared line split the cubic into the line again and a conic. Their axes then lie in the plane through the line normal to both cones, and meet. Both cones touch the sphere centred there, so the common-sphere conics read them first. A guard for a touch within `tol.angular` hands the pair to the search if it reaches this reading.

## A cylinder along a generator through the apex

A cylinder whose axis runs along a generator of the cone, with the apex on its wall, shares that generator with the cone (ADR 0042 left it to the search). It is read the same way, along the cylinder's generators, as another component of the same kind. Along the cylinder `A` vanishes on every generator, and the apex on the cylinder puts `γ` and `β` both to zero on the shared one, at `t₀`. Now the generator's base turns and its direction `D` does not, so `γ = −2β s` gives the one root. With `W(t) − W₀ = 2 sin δ B'` and `B'` the base's derivative taken midway, the difference identity cancels `2 sin δ` from both:

`s = −(B'·Q (W + W₀)) / (2 D·Q B')`.

`β = n·W(t)`, `n = (a₁·a₂) a₂ − cos² φ a₁`, vanishes at `t₀` and once more, at `2h − t₀` with `h` the angle of `n` across the axis. That is the cubic's one point at infinity, so the cubic is always one branch, never fitted.

The apex is read on the cylinder within confusion of its wall. Unlike two cones' apexes, no angular test is needed: the line runs along the cylinder's axis, so its offset from both surfaces does not grow along it, and `A` already holds within `tol.angular`. A cylinder touching the cone along that generator, its axis in the plane of the generator and the cone's axis, has axes that meet, and the common-sphere conics read it first, as for two cones; a guard on the normals there hands it to the search if it reaches this reading.

`ConeCylinderIntersection::LineAndCubic(Line, QuarticCurve)`: the line through the apex along the cylinder's axis, and the cubic. `intersect` returns both, exact and complete, with the apex as the contact point. `QuarticBranch::from_cone_and_cylinder` rebuilds the branch.

## Considered options

- **A cubic type of its own**, rational in the tangent of the half angle: rejected for now. The quartic's walk already evaluates, splits at infinity, rebuilds and pickles. A rational cubic would duplicate all of that for one pair. It is worth revisiting if a consumer needs the cubic's rational form, for example to convert it to a NURBS.
- **The cubic as a `Quartic` alone, the line dropped**: rejected. The line is part of the intersection.
- **The line as an `IntersectionCurve` inside the `QuarticCurve`**: rejected. `QuarticCurve` is a walk along generators, and the line is the one generator it cannot walk.
- **Reading through `along` with the root zero, as for a cone's apex on a sphere**: rejected. `−β / A` loses all precision near `t₀` and is not finite on it.
