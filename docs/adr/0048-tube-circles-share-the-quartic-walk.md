# Tube circles share the quartic walk

A torus met a plane or a sphere in closed form only in special positions, and every other position was marched (#270). Each tube circle of the torus, `C + (R + r cos v) e(t) + r sin v a`, meets a plane or a sphere where `α cos v + β sin v = γ`. It meets it in two points where `Δ = α² + β² − γ² > 0`, and they merge where `Δ` vanishes. For both surfaces `α` and `γ` are linear in `cos τ`, `τ` the angle from the tube circle towards the plane's normal or the sphere's centre, and `β` is constant. So `Δ` is a quadratic in `cos τ`. Its leading coefficient, `(r² − R²)` times the square of that direction's part across the axis, is negative on a ring torus, and `Δ` factors as `|L| (cos τ − c₁)(c₂ − cos τ)`. That is the same form as the discriminant of the quadric quartics (ADR 0040): two factors, each `cos(t − θₖ) − cₖ`, the second read from `θ + π`.

So the components follow from `factored_loops` as they do for a cone and a cylinder: two loops around, arcs between grazing tube circles, a loop crossing itself where a factor is `1 ± cos`, a point where the surfaces touch and part. The surfaces touch only on the tube circle towards the other surface or across from it, where the torus's normal is in the plane of the tube circle. The fit through exact points is the quadrics' too.

The quadrics' walk was written for straight generators: `Rulings` held the ruled quadric and the other one, and computed both the signed root of `Δ` and the point it gives on a generator. We split it into three parts, each shared where it can be:

- `Roots`, the signed root of `Δ` at a parameter of a component: the arcs' exact ends, the sign through a crossing. It does not know what a member of the family is.
- `CurveFamily`, a trait whose provided methods are the fit: segment errors, refinement, distances, the curves in parameters. Its required methods give the point, the tangent and the sample at a parameter. `Rulings` implements it for generators and `TubeCircles` for tube circles. A tube circle's point is `(cos v, sin v)` along `(αγ − β root, βγ + α root)`, on the torus at every root. A tangent through a crossing is a central difference of exact points, as on an arc of the quadrics.
- `SecondSurface`, the other surface's normals, distances and parameters, split from `Other`, which keeps what only a line meets: the distance along a generator. A plane is a `SecondSurface`, and so is the sphere the quadrics already read.

The fit also learned that a plane's parameters are not angles. It unwraps each parameter of each surface across its seam only where it is an angle, and both of the torus's are.

A tube circle within confusion of grazing the other surface on the tube circle towards it is read as touching. The plane or the sphere is moved by that gap so it touches the torus exactly, as a sphere's radius is moved against a cylinder. The point then lies on the curve's crossing. Within confusion on both tube circles, `Δ` vanishes at both ends and is positive between. The curve passes through both points, crossing itself at each, and splits into two circles: the Villarceau circles of a bitangent plane, or two circles of a sphere touching the torus twice. Each runs through both points and one of the two points on the tube circle a quarter turn away, and is checked against the torus before it is returned. A sphere centred in the torus's plane `√(ρ² − r²)` from the two tube circles at `cos τ = −R / q` holds both: there `α` and `γ` vanish together and `Δ` touches zero inside `(−1, 1)`, so that case is read before.

Reading the sign of `Δ` on an arc broke for symmetric curves. `factored_loops` took it at the arc's middle, which for a plane parallel to the axis is the crossing itself, where a touching factor vanishes. A touching factor never changes sign, so it is read as its sign.

A horn or spindle torus (`r ≥ R`) keeps `NotAnalytic` away from its special positions. Its tube circles cross on the axis, where its normal fails, and the leading coefficient turns positive.

The `marching_intersection` bench and the `marching` fingerprints timed these torus pairs. They now march the same plane and sphere against the torus's B-spline form, which no closed form will reach.

## Considered options

- **Generalise `Ruled` into a family both implement**, the point on a member read from a distance or an angle: rejected. The quadrics read the other surface's quadratic along a line (`Other::along`, its rates, the cone's cases off a shared generator and met once). None of that applies to a circle, and the trait would carry it as `unreachable!` arms on one side.
- **A separate walk for the torus**: rejected. The fit is most of the code, and its refinements (exact arc ends, the crossing, the first-order distance) apply unchanged.
- **Solve the quartic in `tan(v/2)` on each tube circle**: rejected. The sine wave gives both points in closed form, with no polynomial to solve.
