# Touching quadrics cross on the walk

A cone and a cylinder, or two cones, that touch at a point meet in a curve that crosses itself there, or at that point alone. ADRs 0040 and 0041 left them to the march (issue #237). Two cylinders that touch from inside already had a closed form: one loop twice around through the crossing (#234). This extends it to the pairs whose discriminant factors with a shift between the factors, and to the crossing that falls inside an arc.

## Where a factor vanishes twice

`Δ = scale · |cos(t − θ₀) − c₀| · |cos(t − θ₁) − c₁|` along the ruled quadric. A factor with `|cₖ|` within its touch tolerance of 1 nearly vanishes twice. Its plane, or its edge of the cone's shadow, touches the ruled quadric along a generator, and on that generator the quadrics touch within confusion. That factor is read as touching exactly: `cₖ` is set to `± 1`, and the crossing is the generator at `θₖ`, or `θₖ + π`. The point where they touch is the double root there, `−β / A`. Setting `cₖ` moves `Δ` by no more than the touch tolerance, chosen so that it moves the points within confusion. Two cylinders move their radii instead, as they can. Both factors touching, the quadrics touching twice, stays with the march. *Amended (ADR 0044):* for two cones the touch is now measured as a distance to the second cone (issue #277).

Where `Δ` is negative on both sides of the crossing, the quadrics touch at an isolated point:

- with no other arc, the result is `TangentPoint`, a new variant of `ConeCylinderIntersection` and `ConeConeIntersection`, as two cylinders and a sphere already have;
- with arcs elsewhere, it is the curve's contact point.

Where `Δ` is positive on both sides:

- **The other factor never vanishes**: every generator meets the other quadric twice. That is the cylinders' node, one loop twice around through the crossing. Its other factor is now read from its own reference, as the factors of a cone and a cylinder, or of two cones, have their own.
- **The crossing inside an arc of the other factor**: the arc runs between that factor's zeros and passes the crossing on its way. Its square root is read signed there (`Ends::Crossing`): `|cos(t − θₖ) ∓ 1| = 2 sin²((t − t_c) / 2)`, so `√Δ` is `√(2 scale g) sin((t − t_c) / 2)`, smooth through zero. The arc then runs out from one end through the crossing on one sign and back through it on the other, as a single loop.

The quadrics' normals are parallel at the crossing, so their cross product is no tangent there: it vanishes and turns over. An arc through a crossing takes its tangent from a central difference of its exact points in the walk parameter, as the node takes it from its velocity. Its points are exact, and the fit is the same as for any other arc.

## Considered options

- **Two arcs meeting at the crossing**: rejected. Each would end at a point where its tangent jumps to the other branch's, and the fit of each would see a corner. Read signed, the loop is smooth everywhere and passes the crossing twice, once on each branch.
- **The exact velocity, as the node has**: deferred. Through the arc's own parameter it is the chain rule through `sin²(p / 2)` and the signed factor. The difference is within the fit's need, and nothing needs more. Profile first.
- **Moving the surfaces to make the touch exact**, as the cylinders do: rejected for cones. No single parameter of a cone moves only the touching factor.
