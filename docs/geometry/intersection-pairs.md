# Closed-form intersection pairs

See the [intersection guide](../guides/intersection.md) for calls,
result fields and numeric fallback.

## Planes and quadrics

Closed-form pairs include plane/quadric/torus surfaces against each other
in their symmetric configurations — among them the tee of two equal
cylinders whose axes meet, two ellipses that cross at the two points
where the cylinders touch, both listed in `contact_points` — and any two
cylinders whose axes cross, a sphere and a cylinder or a cone whose
axis misses its centre, a cone and a cylinder off its axis, or two
cones in general position. Those meet in a quartic, in one loop or two:
the bicylindric curve of two cylinders, or, for a sphere, Viviani's
curve among others. Where a generator of one cone runs along one of the
other's, the curve runs off to infinity there: each piece between two
such generators comes back whole as an exact `QuarticBranch`, unbounded
both ways, the loops that stay bounded fitted as usual. A cylinder whose
axis runs along a generator of a cone meets it once on each of its own
generators, on one loop that is bounded, or runs off to infinity in
branches like those. A cone and a cylinder, or two cones, that touch
at a point meet in a curve crossing itself there, one loop through the
point twice, or at that point alone; so does a cone whose apex lies on
the other surface, at its apex. Where each apex of two cones lies on
the other, the line through both lies on both, and the rest is a cubic
through both apexes: `intersect` returns the line and the cubic's
pieces between its points at infinity as `QuarticBranch` curves, with
both apexes as contact points. So does a cylinder along a generator of
a cone through its apex: that generator and a cubic through the apex,
in one branch, the apex the contact point. Surfaces that touch twice
and meet nowhere else, such as a pipe beside a cone's apex touching
both nappes, give those two points as `contact_points`.

`intersect_exact` returns its closed form
(`QuarticCurve`, exact points on both surfaces), and `intersect` fits
each loop within confusion and reports the result exact and complete,
with the point where the surfaces touch, if any, in `contact_points`:
two cylinders with axes `|r1 − r2|` apart meet in one loop that crosses
itself there.

Two cones with one apex and axes apart split into lines
instead: up to four through the apex, the generators they share, with the apex in
`contact_points`. A cone and a cylinder, or two cones, whose axes meet
at a point as far from both, so that both touch one sphere centred
there, split into two plane conics (the ellipses of a reducer), with the
points where they touch in `contact_points`. Two cones with parallel
axes and equal half-angles, each the other moved sideways, meet in one
plane conic: the hyperbola of two neighbouring chamfers.

## Tori

A torus meets a plane or a sphere in any position in closed form too,
read along its tube circles instead of a quadric's generators. A plane
cuts a spiric section: two loops around the axis, one loop or two
(Cassini's ovals, parallel to the axis at the minor radius from it), or
one crossing itself where the plane touches the inside of the tube, that
point in `contact_points`. A sphere off the axis meets it in a quartic of the
same shapes. A sphere touching it twice meets it in two circles through
both points, and one centred in the torus's plane can hold two tube
circles. A plane or a sphere that touches the torus and meets it nowhere
else gives that point. A horn or spindle torus, whose tube circles cross
on its axis, runs the numeric search away from the special positions.
A cylinder or a cone off the torus's axis meets it in a curve of degree
eight, an `OcticCurve`, read along the same tube circles: each holds up
to four of its points, the roots of a quartic. A pipe through the tube
cuts two loops, one along a diameter four, and one touching the torus
crosses itself there, that point in `contact_points`, as does a cone
whose apex lies on the torus. One touching it and meeting it nowhere
else gives that point. A horn or spindle torus, or a cylinder or a cone
holding a tube circle, runs the numeric search.
Two tori off one axis meet in an `OcticCurve` too, read along the first
one's tube circles: two loops for rings linked like a chain, loops
crossing where tori side by side touch, a point where they touch and
part. Pipe bends joined end to end share a tube circle and touch all
along it: it comes back as a tangent curve (`OcticCurve::shared_circle`),
with the loops of the rest of the curve. A circle of one torus lying on
the other otherwise, a latitude or Villarceau circle, comes back as one
of the loops, the rest of the curve crossing it, in either order: where
a tube circle of the first torus lies on the second, the pair is read
along the second's tube circles instead.

## Curves

Lines and conics intersect each other and elementary surfaces in closed
form, with three exceptions. An ellipse, parabola or hyperbola against
a torus takes the numeric route because substitution gives a degree-eight
equation. A
circle against a torus gives one of degree 4, the torus being
bicircular: up to four hits, or the whole circle where it lies on the
torus (a tube, latitude or Villarceau circle). Any two conics —
coplanar (up to four hits) or in crossing planes (met on the shared
line). The 2D pairs are exact too: any two of `Line2D`, `Circle2D`,
`Ellipse2D`, `Parabola2D` and `Hyperbola2D`, whole or trimmed, solved
by the coplanar solvers with the same tangency multiplicities and
overlaps; only a pair holding a `BSplineCurve2D` or a `BezierCurve2D`
runs the numeric search. A Bézier curve or surface is searched as the
B-spline it is.

## Plane cuts of swept surfaces

A plane meets an extrusion or a revolution in closed form in four
positions: generator lines where it holds the extrusion's direction, the
profile carried onto it otherwise, circles where it is perpendicular to
the revolution's axis, and the profile turned into it where it holds the
axis. The curves are exact types, a `Line`, a `Circle` or a curve of the
profile's kind. They sit where the profile crosses the plane, though,
and that crossing is exact only for an elementary profile. On a
B-spline or Bézier profile it comes from the numeric curve–plane search,
so a horizontal cut of a spline revolution is a list of `Circle`s with
`method` `Numeric` and completeness `Complete`.
