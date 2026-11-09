# A trimmed surface clips the intersection of its basis

A face of a model lies on a patch of its surface, not on the whole surface: two faces on crossing cylinders meet along an arc of the quartic, not the whole loop, and two patches of a cone and a cylinder along one of its generators meet along a bounded piece of a curve that, untrimmed, runs to infinity and makes the march cut it and report `Incomplete` (issue #266, after the box parameter of #264 and #265 was closed). `TrimmedSurface<S>`, built by `surface.trimmed(u_range, v_range)`, is the surface restricted to a parameter rectangle, and `intersect` takes it with any surface, in either order.

The intersection of two surfaces, one or both trimmed, is computed in three steps:

1. **Solve on the bases.** The pair is routed as its two bases are, through the closed-form pair table, so an exact pair stays exact. A march is held `within` the box of the trimmed sides (the intersection of the two boxes when both are trimmed), grown by 1 % of its diagonal so a branch is traced a little past the ranges: a branch ends on the box and is not reported as cut. Boxes that do not meet return an empty exact result without solving.
2. **Clip each curve through its curve in parameters.** A curve's curve in a trimmed side's parameters (the exact one, the march's fit, or one mapped within a tenth of the confusion distance when the result has none) is intersected with the sides of the parameter rectangle, and with its copies shifted by whole periods along a periodic direction. The curve is cut at the crossings and at the ends of a stretch along a side, and each piece is kept when its middle lies in the rectangle, within confusion in space past an end. Pieces adjacent after the cut are merged, a piece shorter than confusion is dropped, and a periodic curve cut on both sides of its seam keeps the piece across the seam whole, its range running past the period. An unbounded curve, a line, is first bounded to where it leaves the trimmed sides' box. A curve that cannot be mapped to parameters, such as one lingering at a cone's apex, is cut where 64 samples of it enter or leave the rectangle, refined by bisection.
3. **Filter the rest.** A contact point is kept when its parameters lie in the ranges of each trimmed side. Two patches of one surface are coincident when a 9 × 9 grid of either rectangle has a point in the other; otherwise the result is empty.

The result keeps the curve of the bases and says which piece of it lies on both patches: `IntersectionCurve` gains `parameter_range`, the curve's own bounds for an untrimmed pair (infinite for a line) and the piece for a trimmed one. One curve of the bases may give several pieces, each its own `IntersectionCurve` sharing the curve. The method and completeness of the bases' result are kept.

## Considered options

- **A box parameter on `intersect`** (#264, #265): rejected and closed. A box is not what a face carries, it does not clip exact curves, and the caller had to compute it.
- **Trim the curve itself, a `TrimmedCurve` in `IntersectionCurve`**: rejected. `Curve` has no trimmed variant, and the field's type would change for every untrimmed result. A range beside the curve keeps the untrimmed result as it was, and `TrimmedCurve::new(curve, range)` builds the piece when one is wanted.
- **Clip by sampling the 3D curve**: kept only as the fallback. Samples miss a curve that leaves the rectangle between two of them; the crossings of the curve in parameters with the rectangle's sides are found by the 2D intersection, which does not.
- **Make `TrimmedSurface` a `SurfaceOperand` with its own basis kind**: rejected. Every exact solver would have to clip, and a trimmed surface of a trimmed surface would be an operand too. The clipping is one function over trait objects, compiled once, after any solver.
- **Shift each curve in parameters by periods into the ranges**: deferred. The curves in parameters stay those of the bases, which the trimmed surface evaluates the same; `Curve2D` has no translation yet.

## Consequences

- Breaking in Rust: a struct literal of `IntersectionCurve` needs `parameter_range`.
- A curve and a trimmed surface do not intersect yet: the Python `intersect` raises `TypeError` for the pair, and Rust has no implementation. A curve against the basis, filtered by the patch's `contains`, does it meanwhile.
- A trimmed surface's projection is the basis's closest point when it lies in the ranges, otherwise the closest point on the boundary; another local minimum inside the ranges is not searched.
