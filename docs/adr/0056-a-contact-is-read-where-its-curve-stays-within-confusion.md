# A contact is read where its curve stays within confusion

The closed forms answer two surfaces near a contact with the topology of the contact: a tangent point, a tangent line or circle, a loop that crosses itself, two conics or two circles that cross. They used to decide it on the offset alone: surfaces within `tol.confusion` of touching were read as touching (#438, #439). Near a contact two surfaces part quadratically, so offset by `ε` along their common normal, the curve where they meet lies `√(2ε / κ)` from the contact point, `κ` their relative curvature across it. For `ε` near confusion that is far beyond it: a plane `5e-8` off a saddle level of a torus of radii 10 and 2 was answered with a crossing `1.3e-3` from the curve it stood for, and a thin cylinder `1e-7` from touching a thick one inside with a crossing `1e-3` from its two loops. Every point of those answers was within `ε` of both surfaces; the answer still was not within confusion of their intersection, nor of its topology.

A contact is now read where the curve of the surfaces as given stays within confusion of it: `√(2 |ε| / κ) ≤ confusion`, that is `|ε| ≤ κ confusion² / 2`. Beyond that, the closed form gives the curve of the surfaces as given: the parted or joined loops, the small loop around the contact, the two lines or circles a hair apart. The rule has one home, `reads_as_touching` in `operations/intersect/exact/contact.rs`, and every closed form asks it:

- **Where the surfaces stand apart**, no curve near the contact, a gap within `tol.confusion` still reads as touching. The contact point is then within confusion of both surfaces and no curve of theirs is left out (ADR 0034). Where they pass through each other, or where the curve lies on both sides of the contact, a crossing that parts one way and joins the other, the offset counts by its size.
- **`κ`** is the smaller principal curvature, in size, of one surface relative to the other at the contact: the eigenvalue nearest zero of the difference of their second fundamental forms read along a common normal. It bounds the curve's reach every way from a point contact. A contact along a line or a circle, of an elementary pair, takes the curvature across it in closed form: `1/r` for a plane on a cylinder or a torus's top circle, `1/r1 ± 1/r2` for parallel cylinders or spheres.
- **An offset within its rounding** reads as no offset: within 16 roundings of the lengths it is computed from and of the contact point's coordinates. An exact contact placed anywhere in the model reads as one, and an offset the rounding of the operands cannot settle is answered with the contact.

The rule is invariant under similarities: `ε` scales with the model, `κ` inversely and `tol.confusion` with it (ADR 0051).

It amends ADR 0043 and ADR 0044 on one point: a factor of the discriminant within its touch tolerance of vanishing twice is read as touching only where it also passes this rule, the gap measured as before (`r (1 − |cₖ|)` for a cone and a cylinder, `|Q| / |∇Q|` for two cones, the tube circle's reach for a torus). Cylinders whose radii move to make a touch exact move them by no more than the rule allows.

## Considered options

- **Keep the offset rule**: rejected. It answers configurations a hair from a contact with a topology and a point thousands of confusions from their intersection, and a caller has no way to tell.
- **Report an `Unseparated` region where the closed form cannot tell the strands apart**, as the curve–surface search does (#429): not needed. Every closed form that snapped already reads the parted or joined loops of the surfaces as given, each point exact; the fit (ADR 0036) follows their pinch like any other turn.
- **The larger relative curvature, or the one on the side of the offset**: rejected for a single number. At a saddle the curve reaches `√(2ε / κ)` with `κ` of the offset's sign, and around an elliptic contact the loop's longer axis is set by the smaller curvature. The smaller one in size bounds both; the band it loses is a factor of the ratio of the curvatures, a few confusions squared over the radius.
- **Curve–curve and curve–surface hits**: unchanged. A hit is a point, and ADR 0034 reads two curves within confusion over a stretch as one contact on purpose.

## Consequences

Within confusion of a contact, a result now depends on which side the surfaces are and by how much, and the loops near a contact may pass within `2√(2 confusion / κ)` of each other. The fits of those loops take more segments near the pinch.

The octic curves of a torus and a cylinder, a cone or a second torus off its axis decide their touching points with a tolerance of their own (`touching_exactly` in `quartic/octic.rs`), and a sphere holding a tube circle of a torus is decided on its offset; both are left to a follow-up.
