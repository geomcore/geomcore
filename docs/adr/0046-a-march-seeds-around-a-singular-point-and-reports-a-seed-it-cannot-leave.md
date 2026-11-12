# A march seeds around a singular point and reports a seed it cannot leave

The march steps along `N1 × N2` and corrects onto both surfaces through their partial derivatives. Where a surface's partials span no plane, neither works. Issue #260 found three ways this lost curves silently. Two cones with one apex, placed at their apex (reference radius 0), had every seed land on the apex, every trace lost, and came back `Sampled` with no curves. The same cones placed at a circle of radius 1 gave their four, three or two lines. A touch at a sphere's pole was lost the same way. A lost trace was dropped without a trace, so the result claimed what the march had not seen.

## A seed no trace leaves is a region

A seed whose trace is `Lost` is kept. Once every seed is traced, each lost seed that no traced curve passes through becomes a region of the result, its parameters and point, with reason `Unseparated`. The reason already covers "a hit at a singular point" for curves. A curve may leave that point, and the march could not tell, so the result says where it could not look instead of claiming nothing is there. A lost seed that a traced curve passes through is explained, and is dropped. A seed met again, within confusion of one lost already, is not traced or seeded around again.

The issue also proposed the same honesty for a pair of unbounded surfaces that gets no patch-pair search. That was not done. Every numeric result is already `Sampled`, which says the search may miss a loop between its seeds. A patch-pair search that did not run is not a known miss at a known place, and an `Incomplete` result needs regions to hold where the miss is.

## Seeds around a singular point

At a lost seed, the march checks whether a surface's partials span no plane there, with one of them vanishing. Only a lost seed is checked: a seed that traces a curve, or a stall, costs the march nothing more. Checking those too was tried, cost up to 0.5% of the budgeted marches, and no test needed it. Where the partials span no plane, the surface near the point is swept by a ring. The parameter of the vanishing partial runs over its range, and the other moves by a first step's length either way, which gives one ring on each nappe of a cone. The distance from the ring to the other surface is sampled 48 times, and each local minimum is refined by Brent's method. A minimum within confusion is a branch crossing the ring, converged onto both surfaces and traced like any seed.

This reads the tangent cone at the point without naming the surface. At a cone's apex the branches leaving it are the generators that lie on the other surface's tangent cone there: the shared generators of another cone through the apex (the closed form of #235), or the generators in a smooth surface's tangent plane. A generator along which the two surfaces only touch is a double direction, where the distance has a minimum of zero without changing sign. It is decided within tolerance: kept where the minimum comes within confusion. A spindle torus's points on its axis are read the same way.

So that a branch through a singular point does not depend on a seed landing there by chance, the known singular points of an elementary surface are seeded too, when they lie within confusion of the other surface. These are a cone's apex, and a torus's points on its axis when its minor radius is at least its major radius. They go after the other seeds, so a point a traced curve already passes through is skipped.

A cone placed at its apex has no size of its own, its reference radius being zero. The seeded trace box was grown by that radius, so it stopped the lines a first step from the apex. Such a cone now counts as having no size, as a plane does, and the box is grown from the seeding windows.

## A pole has a normal

A sphere's pole is smooth, and only its parameters degenerate. `ParametricSurface::normal` already recovers the limiting normal there (ADR 0010), and the crossing direction uses it. So a curve crossing a pole was already traced through it, from the pole or across it. What failed was a touch at the pole: `grazes` read the angle between the surfaces from the second forms, which need the partials, and said "crossing". The crossing direction then vanished, and the trace was lost. `grazes` now falls back, where the second forms are missing, to the normals of the geometry: the touch is read when they are parallel within `tol.angular`. That limiting normal is taken as the surface's normal only when it is the same a third and two thirds of the way round the vanishing parameter's range. At a pole it is. At an apex it turns with the generator, and the apex keeps no normal, so a crossing there stays lost and is seeded around.

A touch at a pole is then a contact point. A touch along a curve through a pole has no contact direction there, from the same missing second forms, so its trace is one contact point. That case is left until an input needs it.

## Considered options

- **Projecting the step's guess onto a surface whose partials degenerate**, so that a step from a pole leaves it: tried, and not kept. Newton's method already left the pole from the unchanged parameters, and no test failed without it.
- **Closed forms for the tangent cone at each kind of singular point** (the apex's generators against a plane, a sphere, another cone): rejected. The ring reads every kind of singular point with one rule, and the closed forms of those pairs already answer them before the march.
- **A region for every unbounded pair without a patch-pair search**: rejected, as above.
