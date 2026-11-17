# Why ten special cases

The [elementary geometry reference](../geometry/elementary.md) describes
ten shapes that a kernel could, in principle, stop distinguishing. Every
one of them — lines and conics included — can be written as a NURBS
curve or surface: a conic as a rational quadratic. Why keep ten types?

Keeping the ten buys four things:

- **Exactness that survives editing.** A radius is stored, not derived.
  Moving a circle's centre or changing its radius edits a number; the
  shape stays a circle instead of becoming a fitted approximation.
- **Closed-form inversion.** Point to parameter is `atan2` and a dot
  product, not a Newton iteration. The same inverse map answers
  containment, projection and extrema.
- **Closed-form intersections.** Plane ∩ cylinder is a conic, not a
  march. The exact solver knows the pair and returns the curve's own
  type, with tangency and coincidence classified by distance.
- **Recognisability.** STEP names these entities, and a hole stays a
  cylinder through a round trip. Downstream tools and users can name
  what they have rather than carrying a generic spline.

The cost is the pairing table in the intersector: every pair of ten
shapes needs an owner, even when many pairs share a solver. The
[intersection guide](../guides/intersection.md) describes how the router
uses that table, and
[ADR 0011](../adr/0011-intersections-route-to-the-exact-solution-else-numeric.md)
records the decision to route each pair to its exact solution when one
exists and to name the numeric fallback otherwise.

B-splines are the other half of the model, not a replacement for this
one. Exact conversions bridge the two: `to_bspline()` preserves every
point of a circle or an ellipse, and the
[NURBS explanation](../nurbs/index.md) covers why rational splines can
represent the conics in the first place.
