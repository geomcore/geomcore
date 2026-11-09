# An apex on the other surface is a crossing on the walk

A cone whose apex lies on a cylinder, or on another cone, meets it in a quartic with a double point at the apex: two branches crossing there, or the apex alone beside a loop elsewhere. ADRs 0040 and 0041 left these to the march (issue #237). Read on the walk they are one more crossing, as where the surfaces touch (ADR 0043).

## The second cone's apex on the walked quadric

The generator of the walked quadric through the other cone's apex meets that cone at the apex alone, twice: `Δ` vanishes there, and every factor of it does. Each edge of the cone's shadow on the cylinder's circle, or each plane through both apexes that touches the second cone, passes through the apex.

That generator is found from the apex: `ρ(t) = −p / |p|` on the cylinder, and along the line from the first apex to the second on a cone, turned over where that line runs into the other nappe. Each factor is then read as vanishing exactly there. With two factors, `cₖ = cos(t_a − θₖ)`, and each factor's other zero is the reflection of `t_a` across its reference. With two squares, the constant of each linear form is set so that it vanishes at `t_a`. The apex is within confusion of the walked quadric, so this moves `Δ` by no more than the touch does in ADR 0043.

- **Two factors**: three zeros, `t_a` twice and one other zero of each factor. `Δ` keeps its sign at `t_a` and changes it at the other two, so one arc runs between those:
  - Positive across `t_a`: the arc runs through it. Each factor vanishes once there, so `√Δ` is read with the sign of `sin((t − t_a) / 2)`. The arc crosses itself at the apex, smoothly on each branch.
  - Negative: the arc runs away from it, and the apex is a point of its own. The tangent plane of the walked quadric there misses the cone but at its apex.
- **Two squares**: `Δ` is positive everywhere but at `t_a`, and both roots meet at the apex. Read with the sign of `sin((t − t_a) / 2)`, the two loops around are one loop twice around (`ThroughApex`), through the apex each time. This is the node of ADR 0043 without a factor to read the sign from.

The apex is the curve's contact point either way. The cone's angle around its axis is undefined at the apex, and its rate along the curve is nought over nought there. So the curve in that cone's parameters is not fitted, as for a sphere through whose pole the curve passes, and tangents are central differences of exact points (ADR 0043). The distance to the curve is measured by sliding along it, as near a sphere's cone apex.

## The first cone's apex on the second

Every generator of the first cone then meets the second at the shared apex point. Read along the second cone instead, this is the case above, with the operands swapped in the result. Each apex on the other cone puts the line through both on both cones: the curve is that line and a cubic, and it stays with the march, as does a cylinder along a generator of a cone through its apex (ADR 0042). *Amended (ADR 0045):* the two cones are read in closed form, as the line and the cubic.

## Touching cones, measured

ADR 0043 read a factor as touching where `|cₖ|` came within a tolerance of one. For two cones that tolerance divided by the distance to `−β / A` on the generator at the factor's reference angle. That distance vanishes where `β` does, which happens far from any touch, and the pair went to the march (issue #277). It was also read on the wrong generator where `cₖ < 0`. Now the factor touches where the generator on which it is least, `θₖ` or across from it, comes within confusion of the second cone at the middle of its roots: `|Q| / |∇Q|`, with `Q = −Δ / A` and `Δ` read as that factor times the other. The bound stays finite where `β` vanishes. It belongs to the factor, so another factor vanishing on that generator, an ordinary arc end, does not make this one touch.

## Considered options

- **Moving the surfaces onto the apex**, as two cylinders move their radii: rejected. The result would report surfaces the caller did not pass, and the factors are as easily read as vanishing there.
- **The apex as two arcs meeting there**: rejected, for the reason of ADR 0043.
- **Fitting the cone's parameters through the apex**: deferred. The angle is continuous along each branch, `v` changing sign, but its derivative is a limit there. No consumer needs it yet.
