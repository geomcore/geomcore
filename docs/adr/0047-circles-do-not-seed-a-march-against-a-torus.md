# Circles do not seed a march against a torus

A march seeds from the points where the iso-curves of an elementary side meet the other surface. Circles were seeded only against a surface the closed-form pair table paired them with, read through `ExactFamily::has_closed_form`. A torus had no circle entry, so a sphere's or a torus's circles never seeded against a torus, nor did a cylinder's or a cone's latitude circles.

Issue #269 gives the circle–torus pair a closed form, a quartic in `tan(t/2)`, and the rule would then seed those circles too. We tried it. On the `marching_intersection/sphere_torus` bench it cost 3.7% more instructions: the quartic solves, plus deduplicating the extra seeds and checking them against the traced curves. The curves came out the same, and so did every fingerprint. Probes on main and on the branch gave the same curves in every case:

- 18 small spheres biting the tube between two of the torus's meridians, as small as a radius of 0.03 sunk 0.005 into the tube;
- small tori biting a tube;
- five pairs of thin tori whose tubes cross;
- a cylinder across a thin tube.

The other seeds of those pairs already find the curves: the torus's meridians against the sphere or cylinder, the cylinder's rulings, and the patch-pair seeds of two bounded surfaces.

So circles still do not seed against a torus. The rule is now stated by surface (`SeedTarget::takes_circle_seeds`) instead of read from the pair table, and `has_closed_form` lost its last caller outside the tests and was removed. Routing tests check the table by the method each pair routes to.

## Considered options

- **Seed the circles, as the table now allows**: rejected. It costs 3.7% on a sphere and a torus and buys no curve we could find. A pair where they would find a loop the other seeds miss reopens the question with that pair as its test.
- **Keep `has_closed_form` for the tests**, behind `#[cfg(test)]`: rejected. The routing tests already assert which pairs route exact, through the method of the result.
