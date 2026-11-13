# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html). The Rust crate and
the Python distribution are versioned in lockstep.

## [Unreleased]

### Breaking

- `Tolerance::new` returns `Result<Tolerance, ToleranceConstructionError>`
  (#306): append `?` or `.unwrap()`. Its error, new, names the scale that
  is not a positive finite number. A tolerance built field by field is
  still not checked.
- `SphereConeIntersection` gains `CircleAndApex(Circle, Point)` (#312).
  A `match` on it needs an arm for the new variant, as for `TwoCircles`
  with the point added to the contacts.
- `BSplineConstructionError` and `BezierConstructionError` gain
  `NonFinite`, and `PoleEditError` gains `NonFinitePole` (#307). A
  `match` on them needs an arm for the new variant.

### Fixed

- `numpy.array(discretization)` is the `(N, 3)` (or `(N, 2)`) array of its
  points, not a 0-d array of one object (#319), and `reversed()` works on
  a discretization and on every intersection and extrema result, which
  iterate and index but raised. Intersection results stay unconvertible,
  as `docs/python.md` now says: they hold hits, not points.
- Every value and shape repr pastes back as an equal object (#320):
  `eval(repr(x)) == x` held for 6 of 32 classes. Shapes printed points
  as tuples the constructors refuse; 3D ellipses, parabolas and
  hyperbolas left out their x direction; B-splines and Bézier curves
  printed no poles; `TrimmedCurve` hid its basis; `Axis` and `Frame` had
  no repr; `BoundingBox` printed `min=`/`max=`; a moved `QuarticBranch`
  printed like the original; `Transform` showed rounding noise. A repr
  is now the constructor called with keywords, or the static constructor
  pickling uses when the keywords cannot state the placement, with lists
  past ten items cut by an explicit `...`.
- Every enumeration member prints as `Class.NAME` (#321).
  Since Python 3.11 an `IntEnum` prints its number, so
  `print(result.completeness.status)` showed `2`, and an f-string showed
  it on every version, while `PointParametrization` printed its name.
  `str()` and f-strings now give `CompletenessStatus.SAMPLED` on every
  enumeration and every Python version. Members still compare equal to
  their integer or string value.
- numpy on the left of a point or vector no longer takes over its
  operators (#325). numpy read the value as a sequence and ran first:
  `numpy.float64(2) * Vector(1, 2, 3)` gave an array, and
  `Point(1, 2, 3) + numpy.array([1.0, 0, 0])` an array, not an error.
  `Point`, `Vector`, `Point2D` and `Vector2D` set `__array_ufunc__ =
  None`, so the number scales the vector, an array on either side of an
  operator raises `TypeError`, and `==` against an array is `False`.
  `numpy.array(p)` converts as before.
- A NaN or infinite pole, weight, knot or data point is refused by every
  B-spline and Bézier constructor and fit, curve and surface, 2D and 3D
  (#307). It built a curve of NaN without complaint: one bad sample in
  six gave `BSplineCurve.interpolate` the knots `[0, nan, nan, 1]`, and
  `approximate_within` called the NaN fit proven within any deviation.
  Each now raises `ConstructionError` (Rust: `NonFinite`), as the conic
  constructors already did; `with_pole` refuses a non-finite pole too.
- A trimmed curve against a revolution or extrusion surface, or any
  surface without a closed form, no longer returns hits beyond its ends
  (#308). The numeric search's Newton iteration is not held to the
  curve's range, and from a seed near an end it reached hits of the
  basis curve beyond it: a segment from x = 0 to 1 met an extruded
  circle of radius 4 at parameter 3.99. A hit now counts only within the
  range, or within `tol.confusion` beyond an end, where it is snapped to
  the end.
- Far from the origin, a march through a crossing of two branches could
  run from one branch onto the other: rounding opens the crossing into
  a bend a few millionths wide, and the turns through it were too small
  to end the trace. A plane through a generator of a hyperboloid 1e7
  from the origin gave three lines where it gives four near it. Where
  the surfaces graze, a turn larger than the sine between them now ends
  the trace whatever its size. The last three sample comparisons with
  an unbounded surface now hold.
- A cone whose apex lies on a sphere centred on its axis meets it in a
  circle and the apex (#312). The apex came back as a second circle, of
  radius 0, and with no contact point, unlike every other apex on the
  other surface. `intersect` now reports one curve and the apex as a
  contact point. A `Circle` of radius 0 stays valid, documented as the
  point at its centre, as the iso curve through a sphere's pole is.
- A negative, zero, NaN or infinite tolerance is refused instead of
  answering as if nothing were close enough (#306): `tol=-1.0` gave a
  circle and a line through its centre an exact, complete, empty
  intersection. `Tolerance::new`, Python's `Tolerance(...)` and a number
  passed as `tol=` raise `ConstructionError` (a `ValueError`) naming the
  scale.
- A curve's point extrema are the local minima of the distance, and only
  those (#309), as the projection guide promised. An ellipse, a
  hyperbola and a parabola listed every stationary point, the farthest
  included: from (4, 1, 0) an ellipse of radii 3 and 2 gave its nearest
  and its farthest point, so `extrema(p)[1]` read as the second closest
  approach was the farthest. A B-spline or a trimmed curve listed both
  ends always. A stationary point is now kept where the distance grows
  on both sides of it, and an end where it grows into the curve.
- The Python docstrings say which way surface normals point and that
  `to_bspline` changes the parametrization (#333). `normal()` is
  `du × dv`: along a plane's `normal`, away from a cylinder's or cone's
  axis, out of a sphere and out of a torus's tube, flipped by
  `reversed_u()` and `reversed_v()`. Each `to_bspline` states whether it
  keeps the parameters and names `parameters_of`/`parameter_of` as the
  bridge when it does not. `x`, `y` and `z` on the points and vectors
  have docstrings; `Ellipse` explains `x_direction` and the major ≥
  minor rule, `Cone.half_angle` its unit and range, and `Tolerance` the
  units and defaults of its three scales. The elementary-geometry page
  gains both sections.
- A `TrimmedCurve` given as the profile of a `RevolutionSurface` or an
  `ExtrusionSurface` names the way out (#328): pass
  `profile.to_bspline()`, the piece as a whole curve, or trim the
  surface. It raised a list of the accepted classes.

## [0.46.0] - 2026-11-12

Tori meet everything in closed form. A torus against a circle, a plane,
a sphere, a cylinder, a cone or a second torus, in any position, is now
intersected exactly instead of marched: the curves of degree eight where
a cylinder, a cone or a torus meets it are read along its tube circles,
up to four points each, and come back as `OcticCurve`, a new exact
curve; pipe bends joined end to end share their tube circle. A result is
now meant not to depend on the model's size or placement (ADR 0051): a
judge checks it on every sample pair, and the curve–curve search, its
extrema, the march with an unbounded surface and the closed forms'
tangency decisions no longer read lengths in model units. The march
passes through a cone's apex and a spindle torus's point. The breaking
changes are new variants of the torus result enums; each entry under
Breaking says what to write.

### Added

- Two tori meet in closed form in any position (#272): pipe bends
  meeting at a corner, a fillet against a bend, rings linked like a
  chain. It used to be only with collinear axes, the rest documented as
  "up to degree 16" and marched.
  - The torus being bicircular, a second torus's equation along the
    first's tube circles is of degree two in each angle, as a quadric's:
    up to four exact points per tube circle, and a curve of degree eight
    read by the same sweep (ADR 0050).
  - `TorusTorusIntersection` gains `Octic(Box<OcticCurve>)` and
    `TangentPoint(Point)`: loops, crossing where the tori touch, or a
    point where they touch and part. A torus touching within confusion
    is moved by one translation or one change of its minor radius.
  - Pipe bends joined end to end share a tube circle and touch all along
    it. It is divided out of the equation and comes back as
    `OcticCurve::shared_circle`, routed as a tangent curve, with the
    loops of the rest; alone, it is a `TangentCircle`.
  - A circle of one torus lying on the other as a latitude or a
    Villarceau circle stays `NotAnalytic`, as does a horn or spindle
    torus off the other's axis.
  - `intersect` fits each loop within confusion, with curves in both
    tori's parameters, and reports the pair `Exact` and complete, in
    Rust and Python. Timed under `intersect_routed/torus_octic/tori_*`
    and `pipe_bends_*`.
  - Fitting a loop through a crossing measures the distance to it to
    first order away from the crossing, sliding along it only where the
    normals are near parallel: 2–3 times cheaper on loops that cross.

- A torus meets a cylinder or a cone in closed form in any position
  (#271): a pipe through a fillet, a drill axis off a bend, a cone
  through a toroidal edge round. It used to be only coaxial; every other
  position ran the numeric search.
  - On each tube circle the other surface's equation is
    `A cos 2v + B sin 2v + C cos v + D sin v + E`, a quartic in
    `tan(v/2)`: up to four exact points. The curve is of degree eight,
    read between the tube circles where a pair of points appears or
    vanishes and the points where the surfaces touch. Both are common
    zeros of two trigonometric polynomials, found by subdivision and
    certified (ADR 0049).
  - `TorusCylinderIntersection` and `TorusConeIntersection` gain
    `Octic(Box<OcticCurve>)` and `TangentPoint(Point)`:
    - two loops where a pipe or a cone passes through the tube, four
      where a pipe runs along a diameter;
    - a loop crossing itself where the surfaces touch, or where a cone's
      apex lies on the torus, the point in `contact_points`;
    - a point where either touches the torus and meets it nowhere else.
  - A cone's apex, or a point where the surfaces touch, within confusion
    is moved to touch exactly, by one translation or a cylinder's radius.
  - `intersect` fits each loop within confusion, with curves in both
    surfaces' parameters, and reports the pair `Exact` and complete, in
    either order, in Rust and Python.
  - `NotAnalytic` stays for a horn or spindle torus, a cylinder or a cone
    holding a tube circle, contact of higher order than a crossing, and
    near-contacts no single move makes exact.
  - Timed under `intersect_routed/torus_octic`.

- A torus meets a plane or a sphere in closed form in any position
  (#270). It used to be only through or across its axis, in a Villarceau
  plane, or with the sphere centred on its axis; every other position
  ran the numeric search.
  - Each tube circle meets the plane or the sphere where
    `α cos v + β sin v = γ`. Its two points, `(cos v, sin v)` along
    `(αγ ∓ β√Δ, βγ ± α√Δ)`, merge where `Δ = α² + β² − γ²` vanishes.
    `Δ` is a quadratic in the cosine of the angle around the axis, which
    factors like the discriminant of the quadric quartics. So the loops
    are read and fitted by the same code, along the torus's tube circles
    instead of a quadric's generators.
  - `TorusPlaneIntersection` and `TorusSphereIntersection` gain
    `Quartic(QuarticCurve)` and `TangentPoint(Point)`:
    - a spiric section: two loops around the axis, one or two loops
      (Cassini's ovals among them, parallel to the axis at the minor
      radius from it), or one crossing itself where the plane touches
      the inside of the tube;
    - a sphere off the axis: a quartic of the same shapes;
    - a point where either touches the torus and meets it nowhere else.
  - A sphere touching the torus twice meets it in two circles through
    both points, and one centred in the torus's plane can hold two tube
    circles. Both come back as `TwoCircles`.
  - `intersect` fits each loop within confusion, with curves in both
    surfaces' parameters, and reports the pair `Exact` and complete, in
    either order, in Rust and Python.
  - A horn or spindle torus keeps `NotAnalytic` away from the special
    positions: its tube circles cross on its axis, where its normal
    fails.
  - The `marching_intersection` bench and the `marching` fingerprints
    now march the torus pairs against the torus's B-spline form
    (`bspline_torus_plane_oblique`, `sphere_bspline_torus`). The closed
    form is timed under `intersect_routed/torus_quartic`.
- A circle meets a torus in closed form (#269):
  `IntersectExact<Torus> for Circle` returns a `CircleSurfaceIntersection`.
  - Along a circle, the torus equation is
    `A cos 2t + B sin 2t + C cos t + D sin t + E = 0`, a quartic in
    `tan(t/2)`. This gives up to four hits, with tangency decided by
    distance as for the other circle pairs.
  - Where the circle lies on the torus (a tube, latitude or Villarceau
    circle), the result is the whole circle.
  - `intersect` reports the pair `Exact` and complete, in Rust and
    Python, where it used to run the numeric search.
  - A march still does not seed from circles against a torus. They cost
    3.7% more on a sphere and a torus and found no curve the other seeds
    missed (ADR 0047).
  - An ellipse, a parabola or a hyperbola against a torus still runs the
    numeric search: substituted into the torus, they give an equation of
    degree 8.

### Breaking

- `TorusTorusIntersection` gains `Octic(Box<OcticCurve>)` and
  `TangentPoint(Point)` (#272), and is no longer `Copy`.
- `TorusCylinderIntersection` and `TorusConeIntersection` gain
  `Octic(Box<OcticCurve>)` and `TangentPoint(Point)` (#271): a `match` on
  them needs the arms. `TorusCylinderIntersection` is no longer `Copy`.
- `TorusPlaneIntersection` and `TorusSphereIntersection` gain
  `Quartic(QuarticCurve)` and `TangentPoint(Point)` (#270): a `match` on
  them needs the arms.

### Fixed

- A curve against a surface without a closed form no longer misses a
  crossing that stands alone (#310). A hyperbola against a revolved vase
  gave 3 of its 4 hits, `Sampled`, since 0.44.0: it ran close to the
  vase around three crossings, the 24 closest samples all crowded there,
  and the fourth showed as one dip of the distance only. The search now
  also starts Newton from every dip of the sampled distance that the 24
  closest samples leave out.
- Marching no longer drops a curve through a point where a surface has
  no tangent plane (#260, ADR 0046). Two cones sharing an apex, each
  placed at its apex (reference radius 0), came back `Sampled` with no
  curves, every seed landing on the apex; they now give the same lines
  as the same cones placed elsewhere: four, three (one where they touch)
  or two. The march seeds the branches leaving a cone's apex, or a
  spindle torus's point on its axis, a first step away from it, both
  where a seed lands there and where the point lies on the other
  surface. A seed that no trace can leave and that no curve passes
  through now makes the result `Incomplete` with an `Unseparated` region
  at it, instead of being dropped. A surface touching a sphere at its
  pole gives the contact point, where nothing was reported.
- `project_point` and `extrema` on a B-spline or Bézier curve far from
  the origin: at 1e7 model units the interior minima were dropped and
  the closest point came back up to four times too far. The rounding
  floor of the search now follows the distances from the query to the
  poles, not the coordinates.
- Closed forms that took two surfaces for tangent, or for apart, by
  their size: a sphere on a cone's axis, two spheres, a torus and a
  cone on one axis, and two tori on one axis. A sphere and a cone
  crossing in two circles came back as one tangent circle a million
  times larger, and two spheres crossing 0.01 short of touching as a
  tangent point a million times smaller. Each is now decided by a
  distance against `tol.confusion`: the surfaces are tangent when they
  are within `tol.confusion` of touching. For the torus this widens the
  tangent case, which was decided by how close two roots were.
- `curve_curve_extrema` lost minima when the curves were scaled or
  moved far from the origin, for 33 sample pairs. Three causes. The
  minimiser's shift of an indefinite Hessian was measured on a
  Hessian mixing an angle and a length, so a parabola scaled by 2^20
  barely moved; it is now measured on the Hessian scaled to a unit
  diagonal. Far from the origin the squared distance rounds long
  before the parameters converge; a full Newton step that stays convex
  and at least halves the gradient is now taken as progress, as in one
  parameter. And a `QuarticBranch` differenced its world points
  for its derivatives, so 1e7 away its second derivative was wrong in
  the third digit; it now differences its points before its placement.
- The intersection of two B-spline or Bézier curves ran out of its
  budget of halvings once the curves sat 1e3 from the origin, where it
  had settled them near the origin: the result went from complete to
  incomplete. Its halvings rounded about the origin, so the edge of a
  tangent contact, where the gap between the curves reaches
  `tol.confusion`, was chased into rounding. The pieces' control points
  are now taken about the first curve's first pole.
- The surface–surface march lost curves, or found extra ones, when the
  surfaces were scaled with `tol.confusion`, for 51 sample comparisons. A
  plane and a ring scaled by 2^20 met in no curve and still said
  `Sampled`. Five causes. Its Newton steps solved on derivatives
  mixing an angle and a length, so at 2^20 the singularity test
  rejected every step; each derivative is now divided by its length, and
  the step measured along the surfaces (an angle's step is taken modulo
  its period, which keeps its precision next to a pole). Each trace dropped the points
  closer than 1e-7 to the one before, so at 2^-20 the polyline cut
  corners and a curve already traced was traced again; it now drops
  those within `tol.confusion`. And where branches cross, a trace could
  turn onto the other branch through steps of a few `tol.confusion`,
  whose turns are rounding noise; it now ends at the crossing, which
  seeds every branch through it. A trace dropped its last step onto a
  parameter bound when shorter than 1e-7, so at 2^-20 curves ended short
  of their bounds and were traced twice; it now keeps it beyond
  `tol.confusion`. And where the surfaces cross at a small angle, a seed
  could sit at the edge of the band where they lie within confusion of
  each other, from where no step reached the crossing; seeds are now
  polished onto it, and two spheres 1e-5 apart no longer come back with
  an `Unseparated` region.
- The march between two unbounded surfaces (a cylinder, a cone, a
  plane, an extrusion, the revolution of a line) sampled its seeds in
  the parameters `[-10, 10]` of each unbounded direction, whatever the
  size of the model: a cylinder cut by the extrusion of a line scaled
  by 2^20 came back as pieces ending on that window, `Incomplete`, where
  at scale 1 it is one closed ellipse. The seeds now span a box around
  a point of the other surface, five sizes of the pair wide (a radius,
  a profile's size, or the distance between the surfaces), or the
  model's cube where the pair has no length of its own. 14 sample
  comparisons now hold.

### Changed

- A result is meant not to depend on the model's size or placement
  (ADR 0051). `tests/similarity_invariance.rs` checks it on every
  sample shape and pair, and
  `tests/fixtures/similarity_invariance_known.txt` lists the
  comparisons that do not hold yet.
- A line and a circle of the plane are decided by the distance from
  the centre to the line, without trigonometry: the typed intersection
  takes about 60% fewer instructions, the routed one about 40%.

## [0.45.0] - 2026-11-10

Cylinders, cones and spheres meet in closed form. Every pair of them,
in any position but one, is now intersected exactly instead of
marched: the space quartic where two of them meet is fitted loop by
loop within confusion, and the pieces that run to infinity come back
as `QuarticBranch`, a new exact curve. The special positions have
readings of their own: surfaces touching at one or two points, an apex
on the other surface, a shared apex, a common inscribed sphere, two
parallel cones of one half-angle. One cone pair stays with the search:
a plane through parallel generators touching both cones (ADR 0042).
`TrimmedSurface` restricts any surface to the patch a face lies on,
and an intersection cuts its curves to the patch. The breaking changes
are new variants of the closed-form result enums and a new field of
`IntersectionCurve`; each entry under Breaking says what to write.

### Added

- `TrimmedSurface`, a surface restricted to a parameter rectangle, the
  patch a face lies on, built by `surface.trimmed(u_range, v_range)` on
  every surface, in Rust and Python (#266). A periodic direction's range
  may cross the seam. It evaluates as its basis, its bounds are the
  ranges, and its projection and bounding box keep to the patch.
  `intersect` takes it with any surface, or another trimmed surface, in
  either order: the bases are intersected, exactly where a closed form
  exists, and each curve is cut to the pieces inside the ranges of the
  trimmed sides, which each `IntersectionCurve` gives as its new
  `parameter_range`; contact points outside the ranges are dropped, and
  two patches of one surface are coincident only where they overlap. A
  march is held to the patches' box, so two patches of unbounded
  surfaces whose intersection runs to infinity are no longer cut and
  incomplete (ADR 0039).

- Two cones in general position meet in an exact quartic instead of a
  marched one (#237, ADR 0041). `intersect`, in Rust and Python, fits
  each loop within confusion, with both curves in parameters, and
  reports the result exact and complete: two loops around the first cone
  when its apex is inside the second, otherwise arcs that end on the
  planes through both apexes touching the second cone, or nothing. In
  Rust, `intersect_exact` returns the closed form as
  `ConeConeIntersection::Quartic(QuarticCurve)`, read along the first
  cone's generators, and cones that miss as `Empty`. Each apex on the
  other cone gives the line through both and a cubic through both
  apexes, `ConeConeIntersection::LineAndCubic(Line, QuarticCurve)`
  (ADR 0045): `intersect` returns the line and the cubic's pieces
  between its points at infinity as exact `QuarticBranch` curves, with
  both apexes, where they cross, as contact points.
- `QuarticBranch`, a new `Curve` variant and Python class: a piece of
  the quartic where two cones meet that runs to infinity at both ends,
  where a generator of one cone is parallel to one of the other's (#237,
  ADR 0042). Two such cones no longer go to the march: `intersect`
  returns the loops that stay bounded fitted as before, and each piece of
  a loop that runs to infinity as an exact `QuarticBranch` with the
  parameter range `(−∞, ∞)`, the result exact and complete;
  `QuarticCurve::branches` gives the pieces. Every point of a branch is
  exact; its parameter runs over the real line, the distance growing
  about as fast. Its derivatives are central differences of exact points,
  its projection a sampled search refined by Brent's method and Newton
  steps, and a trimmed piece's box, discretization and length are
  measured on exact points. It reverses, transforms, has no bounding box,
  B-spline or curve in parameters, and intersects other shapes through
  the numeric search. `QuarticBranch::from_cones(first, second, index,
  tol)` rebuilds one, `QuarticBranch::surfaces` and
  `QuarticBranch::tolerance` give back its surfaces and the tolerance it
  was read with, which rebuilds it as it was (#283), and Python pickles
  it; it rebuilds a branch of the cubic of `LineAndCubic` too. A
  generator of one cone touching the other's directions is one where the branch runs off at both ends, and two
  branches meet at infinity where both roots run off together along one
  generator. Still marched: those two at once, the plane through the
  generators touching both cones.
- A cone and a cylinder whose axis runs along a generator of the cone
  meet in an exact curve instead of a marched one (#237, ADR 0042): each
  generator of the cylinder meets the cone once, on one loop around it.
  Out of reach of the cone's generator along the axis the loop is bounded
  and fitted; otherwise it runs off to infinity where a generator of the
  cylinder is parallel to the cone, and `intersect` returns its pieces as
  exact `QuarticBranch` curves, the result exact and complete.
  `QuarticBranch::from_cone_and_cylinder(cone, cylinder, index, tol)`
  rebuilds one, and the Python `QuarticBranch(first, second, index,
  tol=None)` takes a cone and a cylinder too. Trimmed patches of the
  pair cut those branches at their ranges. The `marching_intersection` bench no longer times this
  pair; `intersect_dispatch` does, as `cone_cylinder/along_*`. With the
  apex on the cylinder, that generator of the cone lies on the cylinder,
  and the rest is a cubic through the apex,
  `ConeCylinderIntersection::LineAndCubic(Line, QuarticCurve)` (ADR
  0045): `intersect` returns the line along the axis through the apex
  and the cubic's one exact `QuarticBranch`, with the apex as the
  contact point.
- A cone and a cylinder off its axis meet in an exact quartic instead of
  a marched one, in about a quarter of the time (#237, ADR 0040).
  `intersect`, in Rust and Python, fits each loop within confusion, with
  both curves in parameters, and reports the result exact and complete:
  two loops around a cylinder steeper than the cone, or the arcs of a
  pipe across it, which run between the edges of the cone's shadow, or
  nothing. In Rust, `intersect_exact` returns the closed form as
  `ConeCylinderIntersection::Quartic(QuarticCurve)`, read along the
  cylinder's generators, and a pipe beside the cone as `Empty`.
- A cone and a cylinder, or two cones, that touch at a point meet in an
  exact curve instead of a marched one (#237, ADR 0043). Where the
  surfaces part on both sides of the point, `intersect` reports it as a
  contact point, alone or beside the loops elsewhere; where they cross
  there, the curve crosses itself, as one loop through the point twice,
  fitted within confusion, the point among its `contact_points`. In
  Rust, `intersect_exact` returns a lone point as the new
  `TangentPoint(Point)`. Surfaces that touch twice and meet nowhere else
  give both points, `TwoTangentPoints(Point, Point)` in Rust: a pipe
  beside a cone's apex touching both nappes, or two cones tangent to the
  same two planes through both apexes, in opposite sectors between them.
- A cone whose apex lies on a cylinder, or on another cone, meets it in
  an exact curve instead of a marched one (#237, ADR 0044): the curve
  crosses itself at the apex, one loop through it twice, or the apex is
  a point of its own beside a loop elsewhere. `intersect` reports it
  exact and complete, the apex among its `contact_points`; the curve in
  the cone's own parameters is left out, its angle undefined at the
  apex.
- Two cones with parallel axes and equal half-angles, such as the
  chamfers of two neighbouring parallel holes, meet in one plane conic
  instead of being marched (#262): each cone is the other moved sideways,
  so their equations differ by a plane. `intersect`, in Rust and Python,
  reports it exact and complete: a hyperbola when each apex is outside
  the other cone, an ellipse when inside, and the line through both
  apexes, marked tangent, when each apex is on the other cone. In Rust,
  `ConeConeIntersection` gains `ParallelSection(PlaneConeIntersection)`.
- A cone and a cylinder, or two cones, whose axes meet at a point as far
  from both surfaces, so that both touch one sphere centred there, meet
  in two plane conics instead of being marched (#236): the ellipses of a
  reducer, or ellipses, parabolas and hyperbolas between two cones.
  `intersect`, in Rust and Python, reports them exact and complete,
  each with its curves in both surfaces' parameters, and the points
  where the surfaces touch along the sphere in `contact_points`. In Rust,
  `ConeCylinderIntersection` and `ConeConeIntersection` gain
  `TwoConics(CommonSphereConics)`, the plane section of the cylinder or
  of the first cone in each of the two planes. The cone–cone conics are
  boxed so the result of every other cone pair stays small, so
  `ConeConeIntersection` is no longer `Copy`.
- Two cones with one apex but different axes meet in the lines through
  it that lie on both, up to four, instead of being marched, which
  found none of them (#235). `intersect`, in Rust and Python, reports
  them exact and complete, each with its curves in both cones'
  parameters and marked tangent where the cones touch along it, and the
  apex in `contact_points`; with no line in common, the apex alone. In
  Rust, `ConeConeIntersection` gains `Generators(SharedGenerators)`, and
  `ApexPoint` covers cones that share only the apex. A line within
  confusion of a cone's apex now maps onto the cone's generator in
  `curve_3d_to_uv`.
- A sphere meets a cylinder or a cone whose axis misses its centre in an
  exact quartic (#232, #233), in Rust and Python, instead of a marched
  one. `intersect` fits each loop within confusion and reports the
  result exact and complete, with the point where the surfaces touch in
  `contact_points`. Against a cylinder the sphere's radius decides
  among two loops around the cylinder, one loop over an arc of its
  generators, Viviani's curve (the loop that crosses itself where the
  sphere touches the cylinder from inside), a tangent point, or nothing.
  Against a cone: a loop on each nappe around the apex inside the
  sphere, a loop through the apex on it, or, on each nappe facing the
  sphere, two loops, a crossing loop, an arc, or a tangent point. In
  Rust, `SphereCylinderIntersection` and `SphereConeIntersection` gain
  `Quartic(QuarticCurve)` and `TangentPoint`, and the cone pair
  `TwoTangentPoints`. A curve through a pole of the sphere has no curve
  in the sphere's parameters, whose longitude jumps there.
- Two cylinders whose axes cross meet in an exact quartic, the
  bicylindric curve, instead of a marched one, in 2–3× less time.
  `intersect`, in Rust and Python, fits each of its one or two loops
  within confusion, with both curves in parameters and the `deviation`
  they reach, and reports the result exact and complete (ADR 0036). In Rust, `intersect_exact`
  returns its closed form as
  `CylinderCylinderIntersection::Quartic(QuarticCurve)`, whose
  `point_at(component, parameter)` lies on both cylinders to rounding.
  Cylinders touching from outside give their contact point
  (`TangentPoint`); axes farther apart than the sum of the radii
  give nothing. Axes `|r1 − r2|` apart, where the thin cylinder touches
  the thick one from inside, give one loop that crosses itself there,
  with the touching point in `contact_points` (#234).

### Breaking

- `IntersectionCurve` gains the public field `parameter_range`, the
  range of its curve along which the surfaces meet: a struct literal
  needs it, or `..` from another curve (#266).
- `ConeConeIntersection` gains `Quartic(QuarticCurve)`,
  `LineAndCubic(Line, QuarticCurve)` and `Empty` (#237): a `match` on it
  needs the arms.
- `ConeCylinderIntersection` gains `Quartic(QuarticCurve)`,
  `LineAndCubic(Line, QuarticCurve)` and `Empty` (#237): a `match` on it
  needs the arms.
- `ConeCylinderIntersection` and `ConeConeIntersection` gain
  `TangentPoint(Point)` and `TwoTangentPoints(Point, Point)` (#237): a
  `match` on them needs the arms.
- `SphereCylinderIntersection` and `SphereConeIntersection` lose the
  variant `NotAnalytic`, every configuration having a closed form, and
  gain `Quartic` and `TangentPoint`, the cone pair also
  `TwoTangentPoints` (#232, #233): a `match` on them needs the arms.

## [0.44.0] - 2026-11-06

The API reads the way a newcomer guesses it. Every Python name is at the
package root, as in the crate, and `Curve` is the Python counterpart of
the Rust enum. Methods and arguments say what they compute in plain
words (`point_at`, `first_derivatives`, `max_length_error`,
`half_angle`, `start_parameter`, `min_corner`), and a plane takes the
`normal(u, v, tol)` every surface has. Results are sequences of their
items. The value classes gain the predicates callers wrote by hand
(`is_parallel`, `angle`, `signed_distance`, `between_frames`). Python
also gets a numpy bridge for bulk queries, and every `tol=` takes a
number. Many changes are breaking; each entry under Breaking says what
to write instead. The modules `geomcore.curves` and `geomcore.surfaces`
and the module functions `curve_curve_extrema` and
`curve_3d_to_uv_within` still work, with a `DeprecationWarning`.

### Added

- Every Python shape is importable from the package root (#203):
  `from geomcore import Circle, Plane, Point`, as every Rust item is at
  the crate root. `geomcore.Curve`, `geomcore.Curve2D` and
  `geomcore.Surface` name the union of the classes of the Rust enum of
  the same name, to annotate with and, from Python 3.10, to pass to
  `isinstance` (`typing.get_args` gives the tuple on 3.9). As in Rust, a
  `TrimmedCurve` is a piece of a curve, not one of them: a method that
  takes either is typed `Curve | TrimmedCurve`. The stubs are one file.
  ADR 0038.
- Python 3D curves take a curve in `extrema(other, tol)` (#201), as
  well as a point. A curve or `TrimmedCurve` gives the
  `CurveCurveExtrema` of the pair, holding this curve as `first`. A
  point gives its projections as before. Rust keeps the free function
  `curve_curve_extrema` (ADR 0033).
- Python surfaces have `curve_3d_to_uv_within(curve, max_deviation, tol)`
  and `uv_curve_to_3d_within(uv_curve, max_deviation, tol)` methods
  (#185), beside the exact `curve_3d_to_uv` and `uv_curve_to_3d`. The
  Rust free functions are unchanged (ADR 0035, amended).
- B-spline curves built and read the textbook way (#184), in Rust and
  Python, on both tiers. `from_flat_knots` takes the flat knot vector and
  `flat_knots()` gives it back. `clamped_uniform(degree, poles)` spaces
  clamped knots evenly over `[0, 1]`. `with_pole(index, point)` and
  `with_weight(index, weight)` return the curve with one pole moved or
  reweighted (Python indexes like a list and raises `IndexError`; in Rust
  the new `PoleEditError`). `geometric_continuity(u, tol)` compares the
  unit tangents and curvature vectors on each side of a knot, so a knot
  of multiplicity `degree` between poles in line reports G1 or G2 where
  `continuity` reports C0. The Python repr shows the knots and
  multiplicities of a curve with at most 8 knots.
- Every Python `tol=` takes a number as well as a `Tolerance` (#183):
  `tol=1e-3` is `Tolerance(confusion=1e-3)`, with the default angular and
  parametric tolerances. Any other type raises a `TypeError` naming what
  `tol` accepts. The intersection guide says how the confusion sets the
  size of a fitted intersection curve.
- Value protocols in Python (#182). Points and vectors unpack, index and
  convert to numpy (`x, y, z = p`, `p[-1]`, `numpy.array(points)` of shape
  `(N, 3)`); vectors divide by a number and give their length to `abs`.
  Every value, curve and surface pickles and deep-copies to an equal
  object, so shapes reach `multiprocessing` workers. `Transform` has
  `from_matrix` and `to_matrix` for 4×4 homogeneous matrices, `inverse()`,
  a `repr`, and `a @ b` for `b.then(a)`, the matrix product's order.
  `Circle2D.from_frame` builds a clockwise circle. In Rust:
  `Transform::from_matrix`, `to_matrix` and `inverse`, the new
  `TransformMatrixError`, and `Vector / f64` and `Vector2D / f64`.
- `tight_bounding_box(tol)` on every curve and surface (#178), in Rust
  (`HasBoundingBox`, `HasBoundingBox2D`) and Python: the box that holds
  the shape with each side within `tol.confusion` of it, for collision
  pre-checks. `bounding_box()` stays the box of a B-spline's control
  points, which may be several times wider; the tight box halves the
  Bézier pieces until their control points come that close, in a few
  microseconds for a curve. Elementary shapes' boxes already touch them;
  a revolution surface's tight box comes from its profile's.

- `curve_3d_to_uv_within` maps a curve that starts or ends at a sphere's
  pole or a cone's apex (#220), in Rust and Python: the UV curve takes the
  limit of its point and derivatives from inside the curve, so a meridian
  keeps its own `u` up to the pole and a B-rep edge to a pole vertex gets
  its pcurve. `SingularParametrization` now only refuses a curve that
  lingers at the point.

- `curve_3d_to_uv` maps the iso-curves of swept surfaces in closed form
  (#221), in Rust and Python. On an extrusion, a generator line maps to a
  line along `v` and the profile translated along the direction to a
  line along `u`; on a revolution, a parallel circle maps to a line along
  `u` and the profile turned about the axis to a line along `v`. The
  plane cuts of swept surfaces now carry both `first_uv_curve` and
  `second_uv_curve`, the generator lines through a spline profile's
  numeric crossings included. The carried profile on an oblique plane
  still has none on the swept side: it climbs the surface as it runs.
- Conic and cone conveniences for teaching and drafting (#189), in Rust
  and Python:
  - `eccentricity()` on every ellipse, hyperbola and parabola; `foci()` and
    `directrices()` on ellipses (no directrices for a circle) and
    hyperbolas, ordered so each directrix goes with its focus; `focus()`
    and `directrix()` on parabolas;
  - `Hyperbola::other_branch()` and `Hyperbola2D::other_branch()`, the
    branch through the center that no parameter of the first one reaches;
  - `Parabola::from_focus_and_directrix(focus, &directrix)` and its 2D
    twin, the parameter running along the directrix;
  - `Cone::from_apex(apex, axis_direction, half_angle)`;
  - `geomcore.__version__` in Python;
  - the elementary geometry chapter shows how `Plane.uv_curve_to_3d`
    places a 2D curve in space.
- Every surface intersection curve reports how far it strays (#188), in
  Rust and Python: `IntersectionCurve.deviation`, the largest distance
  between the curve's point and each surface's point at its UV curve,
  in the same parameter. That is the tolerance an edge along the curve
  needs. A marched curve samples it at the middle of each knot span,
  where its fits stray furthest (about 1e-8 at the default tolerance);
  an exact curve is `0.0`.
- Edges as trimmed curves (#187), in Rust and Python:
  - arcs from points: `TrimmedCurve::arc_through(start, middle, end)`
    and `TrimmedCurve::arc(center, start, end, normal, tol)`, with
    `TrimmedCurve2D` twins that take a `counterclockwise` flag in place of
    the normal. They fail with the new `ArcConstructionError`
    (`ConstructionError` in Python), and an `end` confused with `start`
    makes the full turn;
  - an arc's `sweep_angle()` and turning sense: the circle's `normal()`,
    and the new `Circle2D::is_counterclockwise()`;
  - `length(max_length_error)` on every trimmed curve, exact on segments and
    circular arcs, and on `Ellipse` and `Ellipse2D`; `Circle2D.length()` in
    Python;
  - `ParametricSurface::trimmed_curve_3d_to_uv`: the UV curve of an edge,
    its basis's image over the same range. Python's `curve_3d_to_uv` takes
    a `TrimmedCurve` and returns a `TrimmedCurve2D`, and
    `TrimmedCurve.basis()` is typed as a whole curve rather than `Any`.
- Two equal cylinders whose axes meet, the tee of two equal pipes, meet
  in two exact ellipses (#175), in Rust and Python, instead of four
  marched B-splines of some 560 poles each: one in each plane through the
  meeting point that bisects the axes, at any angle between them. Both
  points where the cylinders touch, where the ellipses cross, are in
  `contact_points`; the march found only one.
- Every operation within a deviation reports the deviation it reached
  (#199), in Rust and Python. One result, `CurveApproximation` (3D) or
  `CurveApproximation2D`, comes back from `uv_curve_to_3d_within`,
  `curve_3d_to_uv_within` and `approximate_within` on B-spline curves;
  `SurfaceApproximation` from `BSplineSurface::approximate_within`. Each
  gives its `curve` (or `surface`), its `deviation`, and `is_proven()`. A
  fit measures its deviation at every point, so it is proven.
- Surface intersection curves carry their curves in parameters (#193), in
  Rust and Python: `IntersectionCurve.first_uv_curve` and
  `second_uv_curve`, the curve in each surface's `(u, v)` in the same
  parameter as `curve`. A marched curve always has both, fitted on its
  knots from the parameters the march walked, continuous across seams;
  an exact curve has the closed form where one exists, and `None`
  otherwise (map it with `curve_3d_to_uv_within`).
- `curve_3d_to_uv_within(surface, curve, max_deviation, tol)` (#186), in
  Rust and Python: the UV curve of a 3D curve with no exact one, as one
  degree-5 `BSplineCurve2D` in the curve's parameter, continuous across the
  surface's seams, whose image on the surface is within the deviation of
  the curve at the same parameter. The mirror of `uv_curve_to_3d_within`,
  it returns a `CurveApproximation2D` (its `curve` and sampled
  `deviation`) or a `Curve3DToUvApproximationError`.
- A numpy bridge for bulk queries in Python (#181): every curve and
  surface gains `points_at_array` and `project_points_array`, which take an
  array-like and return float64 arrays (`(N, 3)` or `(N, 2)` points, and
  the parameters and distances of the projections) with the interpreter
  lock released during the computation. The B-spline fitting methods take
  an `(N, 3)`, `(N, 2)` or `(n, m, 3)` array in place of point objects.
  numpy stays optional: it is imported when one of these is called, and
  `pip install geomcore[numpy]` installs it.
- `PointParametrization` in Python (#202): `CHORDAL` or `CENTRIPETAL`,
  taken by `interpolate`, `approximate` and `fit_within` on B-spline
  curves and surfaces. A member equals its lowercase name, and the strings
  `"chordal"` and `"centripetal"` are still accepted; the stubs type them
  as literals, so a typo is caught before the call.

- `axis()` on `Circle`, `Cone` and `Torus`, as `Cylinder` and
  `RevolutionSurface` have it, so the axis of any rotational shape is one
  `Axis` (#196). A circle's axis runs along its `normal`.

- `CurveSurfaceOverlap`, in Rust and Python: a curve–surface overlap is a
  typed item with a `curve_range`, as a curve–curve overlap is a
  `CurveCurveOverlap`, so it can gain a field later without breaking
  callers (#200).

- Python: the results are sequences of their primary items. `len()`,
  iteration and indexing (slices included) run over the hits of
  `CurveCurveIntersection`, `CurveCurveIntersection2D` and
  `CurveSurfaceIntersection`, the curves of `SurfaceSurfaceIntersection`,
  the minima of `CurveCurveExtrema` (which had only `len()`) and the points
  of `Discretization` and `Discretization2D`. The hit-based intersections
  gain `points`, the hit points in the order of `hits`. Truth is
  unchanged: an intersection is true when it holds anything, so one with
  only an overlap or a contact point is true with length 0 (#191).

- `CurveCurveExtremum` carries its two points, `first_point` and
  `second_point`, next to its parameters and distance (Python too), so
  drawing or reporting the closest approach no longer evaluates both
  curves again. `CurveCurveExtrema::closest()` (Python `closest()`) gives
  the closest minimum, or `None`. Breaking in Rust only for code that
  builds a `CurveCurveExtremum` literal, which now needs the two points
  (#192).

- `is_u_periodic()` and `is_v_periodic()` are methods of every surface,
  on the `ParametricSurface` trait (Python: the shared surface methods),
  where `u_period()` and `v_period()` already were. They were
  `BSplineSurface` methods only; in Rust, calling them now needs the trait
  in scope (`use geomcore::prelude::*`).
- `BSplineCurve`, `BSplineCurve2D` and `BSplineSurface` gain
  `to_bspline()`, returning the shape itself, as the other curves and
  surfaces with an exact B-spline form give theirs, so converting code
  needs no special case for a B-spline (#197).

- The value classes answer the geometric questions callers wrote out by
  hand (#198, first part). The caller gives the tolerance: an `angular`
  bound (Python: the default tolerance's when omitted) or a `Tolerance`.
  - `Vector` and `Vector2D`: `angle(other)`, the unsigned angle in
    `[0, π]`; `is_parallel(other, angular)` (either orientation) and
    `is_perpendicular(other, angular)`, which bound the sine and the cosine
    of the angle for vectors of any length and are false for a zero vector;
    `component_along(direction)` and `component_perpendicular_to(direction)`.
    `Vector2D` also gets `signed_angle(other)` and `rotated(angle)`.
  - `Axis` and `Axis2D`: `point_at(t)`, `parameter_of(point)`,
    `closest_point(point)`, `distance(point)`, `contains(point, tol)`,
    `is_parallel(other, angular)` and `reversed()`, which was crate-private.
    `Axis` also gets `is_coaxial(other, tol)`: the same line, in either
    orientation.
  - `Plane::signed_distance(point)`: positive on the side the normal points
    to.
- The placement and container classes get the methods the same callers
  were missing (#198, second part).
  - `Point` and `Point2D`: `midpoint(other)`, and in Python `components()`,
    which `Vector` has. `Point` also gets `transformed(transform)`, so a
    point moves with the same spelling as a curve or a surface.
  - `Frame`: `vector_at(u, v, w)` and `local_components(vector)`, the
    direction counterparts of `point_at` and `local_coordinates`, and
    `transformed(transform)`, which keeps local coordinates.
  - `Frame2D`: `point_at`, `local_coordinates`, `vector_at`,
    `local_components`, and `x_axis()` and `y_axis()`.
  - `Transform`: `between_frames(source, target)`, the rigid motion from one
    placement to another; `scale_factor()`, which was crate-private;
    `preserves_orientation()`, false for mirrors; and `is_identity(tol)`.
  - `BoundingBox` and `BoundingBox2D`: `including(point)`, which was
    crate-private; `intersection(other)`, `None` for boxes apart; and
    `distance(point)`. In Python, `point in box` is `box.encloses(point)`.
- `start_point()` and `end_point()` on every bounded curve, through the new
  traits `BoundedCurve` and `BoundedCurve2D` (in the prelude): the circle
  and the ellipse (whose ends meet at the seam), the B-spline and Bézier
  curves, and the trimmed curves, which had them alone (#198, third part).
  In Python the `discretize` method set becomes the bounded-curve set,
  with the ends added.

### Deprecated

- Python: the modules `geomcore.curves` and `geomcore.surfaces` (#203).
  Import each shape from `geomcore` instead. Every name still resolves
  there with a `DeprecationWarning`, and the modules go in the next
  release.
- Python: `geomcore.curve_curve_extrema(first, second)` (#201). Call
  `first.extrema(second)` instead. It warns with a `DeprecationWarning`
  and goes in the next release.
- Python: the module functions `geomcore.curve_3d_to_uv_within(surface,
  curve, ...)` and `geomcore.uv_curve_to_3d_within(surface, uv_curve,
  ...)` (#185). Call the surface methods of the same names instead. The
  functions warn with a `DeprecationWarning` and go in the next release.

### Changed

- Python: the native extension module is `geomcore._native` (#203). The
  undocumented `geomcore.geomcore`, a second path to every name, is gone.
  Every class's `__module__` is `geomcore`.
- Python: the argument of a 3D curve's `extrema` is named `other`, not
  `point`, since it takes a curve too (#201). A call that passes
  `point=` by keyword needs `other=` or a positional argument.
- An axis, frame or extrusion direction that is already of unit length
  keeps its bits (#182): normalizing it again only moved its rounding, so
  a frame rebuilt from its own directions was a different frame in about
  half the cases. Rebuilding a shape from its parts now gives an equal
  shape. A closed form fed a unit direction can move by an ulp: the
  contact points of the tee of two equal cylinders are now 1 ulp off
  `±5` instead of exactly on it.
- `point_at` documents what a shape gives past its bounds (#179), in
  Rust and Python: a non-periodic B-spline or Bézier
  direction clamps to the nearer bound, point and derivatives alike;
  elementary shapes continue along their equation, and a trimmed curve
  evaluates its basis. Operations that build from a parameter, such as
  `split_at` and `u_iso`, raise on one outside the bounds. ADR 0037
  records why evaluation clamps rather than extrapolates.
- Where a result is numeric, snapped or normalized is documented (#180):
  a plane cut of a swept surface with a spline profile is made of exact
  `Circle`s or `Line`s but reports `method` `Numeric`, because the
  profile's crossings are numeric; `project_point` snaps a distance
  within `tol.confusion` to `0.0`, now said in the Python docstrings
  with the default value; and `ExtrusionSurface` drops the length of
  its direction, now said in the Python class docstring.

### Breaking

- `Plane::normal()` (Python `Plane.normal()`), the plane's
  zero-argument normal, is replaced by `Plane::axis()` (Python
  `Plane.axis()`): the axis through the frame's origin along the oriented
  normal. Write `plane.axis().direction()` where you wrote
  `plane.normal()`. The old method hid the `normal(u, v, tol)` that every
  surface has, so in Python a loop over surfaces raised `TypeError` on a
  plane, and in Rust `plane.normal(u, v, tol)` did not compile. Both now
  work on a plane and return the same direction at every point (#190).

- `IntersectionCurve` has a new public field, `deviation` (#188), so a
  struct literal must set it (`..piece.clone()` keeps working).
- Every intersection hit carries all its parameters (#177), in Rust and
  Python: `first_parameter` and `second_parameter` on `CurveCurveHit` and
  `CurveCurveHit2D`, `curve_parameter` and `surface_parameters` on
  `CurveSurfaceHit`, are now `f64` and `(f64, f64)` (Python `float` and
  `tuple[float, float]`) instead of optional. Where a closed form solves
  for one operand's parameter only, the other is recovered in closed form
  from the point, so a caller no longer calls `parameter_of(hit.point)`
  after the intersection. The hits are always sorted by the first
  operand's parameter.
- `CylinderCylinderIntersection` gains the variant `TwoEllipses` (#175):
  a `match` on it needs the arm.
- The approximation results are one type per tier (#199):
  `UvCurveApproximation` becomes `CurveApproximation`, and its error
  `UvCurveApproximationError` becomes `UvCurveTo3DApproximationError`,
  which mirrors `Curve3DToUvApproximationError`. `approximate_within` on
  `BSplineCurve`, `BSplineCurve2D` and `BSplineSurface` returns a
  `CurveApproximation`, `CurveApproximation2D` or `SurfaceApproximation`
  instead of the spline: read it with `.curve()` or `.into_curve()` in
  Rust, `.curve` in Python (`.surface()`, `.into_surface()` and
  `.surface` for a surface).
- `IntersectionCurve` gains the public fields `first_uv_curve` and
  `second_uv_curve` (#193): a struct literal building one needs them. A
  surface intersection read in the other order of the surfaces now
  differs in them, so it no longer compares equal.
- `curve_3d_to_uv_numeric` is removed (#186): use
  `curve_3d_to_uv_within`, which takes a deviation instead of a sample
  count and returns one continuous UV curve instead of `(u, v)` samples
  that wrapped at the seam.
- The Python enumerations are real enums (#202): `Side`,
  `IntersectionMethod`, `CompletenessStatus` and `IncompleteReason` are
  `enum.IntEnum` classes, so `list(IncompleteReason)` lists the members.
  Members keep their names, values and comparison with integers; their
  `repr` is now the standard `<Side.ABOVE: 1>`.
- Evaluation reads in plain words (#202), in Rust and Python, on every
  curve and surface: `eval_point` and `eval_points` become `point_at` and
  `points_at` (as `Frame::point_at` already was), `eval_derivative`
  becomes `derivative`, and `eval_d1` and `eval_d2` become
  `first_derivatives` and `second_derivatives`; the `_sided` variants
  follow. Benchmark ids keep their names.
- Methods say what they compute in plain words (#202), in Rust and
  Python: `Vector2D::perp()` becomes `perpendicular()`; `Parabola::focal()`
  and `Parabola2D::focal()` become `focal_distance()`, as do the Python
  keywords of their constructors and `from_frame`, and their `repr`;
  `Circle::circumference()` becomes `length()`; `Circle::disk_area()` and
  `Ellipse::area()` (and the Python `Ellipse2D.area()`) become
  `enclosed_area()`. `Transform::then` takes `other`, which no longer
  shadows Python's `next`.

- One accuracy argument, one name (#194). The bound on a length's error is
  `length(max_length_error)` (was `max_error`) on the B-spline and Bézier
  curves, and the fixed size of a least-squares fit is
  `approximate(points, degree, pole_count)` (was `num_poles`), as the
  getters and the surface fit already say `pole_count`. A bound on the
  distance to the true geometry stays `max_deviation`, and `discretize`
  keeps `max_chord_deviation` and `max_angle_deviation`. Only Python
  callers passing these by keyword need to change.

- A trimmed curve names its range by its ends, as it names its points
  (#195): `TrimmedCurve.first()` / `last()` become `start_parameter()` /
  `end_parameter()`, beside `start_point()` / `end_point()`, in Rust and
  Python. The Python constructor keywords follow
  (`TrimmedCurve(basis, start_parameter, end_parameter)`), as does the
  `repr`, and `extract`, `extract_u` and `extract_v` take `start` and
  `end`. `first` and `second` now only ever mean the two operands of a
  result.

- A shape's axis has one spelling (#196). `Cone::semi_angle()` and
  `ref_radius()` become `half_angle()` and `reference_radius()`, in Rust
  and Python, as do the Python keywords of `Cone(center, axis_direction,
  half_angle, reference_radius)` and `Cone.from_frame`. The Python `Torus`
  takes `axis_direction` (was `normal`), as `Cylinder` and `Cone` do; the
  circle keeps `normal`, its own word. Every elementary shape's `repr`
  now uses its constructor's keyword names (`Cone(center=…,
  axis_direction=…, …)`, `Plane(point=…, normal=…)`), which a parity test
  checks.

- Results name the same thing one way (#200).
  `CurveSurfaceIntersection::overlaps()` returns `CurveSurfaceOverlap`
  items instead of bare `(start, end)` tuples; read the range as
  `overlap.curve_range`. In Python, `CurveSurfaceHit.surface_parameters`
  becomes `u` and `v`, as on `SurfaceProjection`; the Rust field stays
  the pair `surface_parameters`.
  `BoundingBox` and `BoundingBox2D` name their corners `min_corner()` and
  `max_corner()` (were `min()` and `max()`), in Rust and Python, as every
  intersection region does.

- Rust: `TrimmedCurve::start_point` and `end_point` are now methods of
  `BoundedCurve` (`BoundedCurve2D` in the plane), so a caller brings them
  into scope with `use geomcore::prelude::*`, like every shared method
  (ADR 0015).

### Fixed

- An unbounded line, parabola or hyperbola against a revolution or
  extrusion surface missed every hit more than 10 parameter units from
  its origin (#174). A line along x from x = -10 found the near wall of
  a revolved vase and not the far one, and a revolved circle gave 2 of
  the 4 crossings that the equal torus gives. The search now seeds the
  parameters whose points lie in the revolution surface's box, or in the
  prism that an extrusion's profile box sweeps along its direction.

- Python: `surface.uv_curve_to_3d(uv_curve)` accepts a `TrimmedCurve2D`, as
  its stub and `uv_curve_to_3d_within` already did, and returns the
  `TrimmedCurve` over the same range of its basis's trace. A B-spline or
  Bézier basis is cut to the piece first, so only the piece has to lie in
  the surface's domain. The stubs of `ExtrusionSurface` and
  `RevolutionSurface` no longer promise a `TrimmedCurve` profile: the
  surface spans its profile's bounds, and `piece.to_bspline()` gives a
  whole curve for a bounded profile. A new test passes every class of
  every curve and surface `Union` in the stubs to the runtime (#176).

- Python: the docstring of `Vector.angle_with_reference` said a zero
  reference gives the unsigned angle; it gives 0 or π. `angle(other)` is
  the unsigned angle.

## [0.43.0] - 2026-11-05

A UV curve without an exact 3D curve, such as a helix on a cylinder, now
maps to one B-spline within a deviation. Unbounded marching branches are
cut where their trace leaves the seeds' box, and reported, instead of
running to a fixed cube around the origin (#147). Curve–curve tangent
contacts are complete, and close crossings stay apart (#162, #163). One
breaking change: `IncompleteReason` gains the variant `Unbounded`.

### Added

- `uv_curve_to_3d_within(surface, uv_curve, max_deviation, tol)` gives the
  3D curve of a UV curve within a deviation, for the pairs that
  `uv_curve_to_3d` refuses: a slanted line on a cylinder or cone (a helix
  or conical spiral), or a circle in the `(u, v)` plane of a sphere
  (ADR 0035). It returns one B-spline for the one UV curve, over its
  bounds and in its parameter. On an elementary surface along a UV curve
  with no knot, the composition is analytic, and the curve is of degree 9
  and C⁴, unless degree 5 takes no more poles. Everywhere else it is of
  degree 5 and C² wherever the UV curve is. Every point at `t` lies within
  the deviation of the surface's point at `c(t)`, the same-parameter
  deviation an edge needs. The result, `UvCurveApproximation`, reports the
  deviation reached and whether it is proven. It is proven by a Hermite
  remainder bound for a straight UV curve on a plane, cylinder, cone,
  sphere or torus, and sampled densely for every other pair. A deviation
  that cannot be reached is refused with
  `UvCurveApproximationError::DeviationNotReached`, not loosened. So is a
  3D curve that collapses to a point, an unbounded UV curve, or one that
  leaves the surface's domain. On a radius-10 cylinder at `1e-7`, a helix
  takes 5 spans and about 23 poles per turn.
- Python: `uv_curve_to_3d_within(surface, uv_curve, max_deviation,
  tol=None)`, taking any 2D curve or `TrimmedCurve2D`, its result
  `UvCurveApproximation` (`curve`, `deviation`, `is_proven()`), and the
  `UvCurveApproximationError` exception.
- `GenericTrimmedCurve<P, C>`: the trimmed curve of either tier, written
  once (ADR 0014, amended). `TrimmedCurve<C>` and `TrimmedCurve2D<C>` are
  now its aliases, with the same defaults, methods and paths.
- `IncompleteReason::Unbounded` (Python `IncompleteReason.UNBOUNDED`): a
  branch of two unbounded surfaces runs beyond the box where its trace was
  cut. The regions of the result hold where each branch was cut. When
  results merge, `Unbounded` ranks below `BudgetExhausted` and above
  `Unseparated`.

### Changed

- `eval_derivative` on a plane, cylinder, cone, sphere or torus evaluates
  every order `(du, dv)`, from the closed forms of the derivatives of
  `cos` and `sin`. It panicked past `du + dv = 2`. The other surfaces
  still evaluate up to second order. Python's `eval_derivative` follows,
  and still raises `ValueError` past `du + dv = 2` on the others.
- The `Debug` output of a trimmed curve names `GenericTrimmedCurve` and
  shows its tier, as the other generic curves' does.
- The trace box of an unbounded pair is the box of the seeds, grown by the
  larger surface size: a cylinder's or cone's radius, a sweep profile's
  diagonal, or for two planes the larger of the seeds' spread and the
  size of the seeding windows. It replaces the cube of
  side `Tolerance::MODEL_SIZE` (1e5) centred on the origin. A pair with a
  bounded side keeps its box of that side, grown by a tenth of its
  diagonal.
- Seeds of an unbounded side whose other side has a bounding box are read
  from that box's corners, with a margin, instead of the window `[-10, 10]`.
- A tangent contact between two curves, where they stay within
  `tol.confusion` of each other and then part, is one hit at its closest
  approach, with multiplicity 2, and the numeric curve–curve intersection
  reports the result `Complete`, as the closed forms do. Before, it was
  `Incomplete` with an `Unseparated` region. Overlaps and searches that
  run out of budget stay `Incomplete`.
- The piece searches' private parameter range is merged into the
  parameter domain the numeric searches already shared: one value, a
  period or an interval, owns the wrapping and the slack. Nothing public
  changes, and no result moves.
- CI gates the benchmarks on an instruction budget: a pull request that
  makes a budgeted benchmark more than 1% heavier fails, and the log names
  the functions that got heavier (`docs/benchmarks.md`). The intersection
  results are fingerprinted, so a change in them shows up in review.

### Breaking

- `IncompleteReason` gains the variant `Unbounded`. The enum is not
  `#[non_exhaustive]`, so an exhaustive `match` on it needs a new arm.

### Fixed

- A cone against a cylinder whose intersection has an unbounded branch ran
  to the cube and was reported as `Sampled`, as if the branch were whole.
  It is cut at its box and reported as `Incomplete` with `Unbounded`.
- An unbounded pair lying beyond the cube lost every branch beyond it. The
  same pair translated far from the origin now traces the same branches.
- The numeric curve–curve intersection no longer drops one of two
  crossings lying within `sqrt(tol.confusion)` of each other and reports
  the result `Complete`. Two hits are merged only when the curves stay
  within `tol.confusion` of each other between them.

## [0.42.2] - 2026-11-02

Two internal moves. Nothing public changes: the types, aliases, methods
and paths are the same, and so is the Python package.

### Changed

- The curve–curve distance extrema have their own module,
  `operations/extrema/` (ADR 0033): the minima search used to live inside
  the numeric intersection, which now reads its hits from it. Nothing
  public changes.
- The private spline behind the B-spline curve is folded into
  `GenericBSplineCurve`, which held it with the same type parameter and
  forwarded every method to it (ADR 0014, amended). Nothing public changes.

## [0.42.1] - 2026-11-02

The source tree, grouped by layer. Nothing public changes: every item
keeps its one path at the crate root (ADR 0022), and the Python package
is the same.

### Changed

- The private modules under `src/` are grouped into `base/` (points,
  vectors, frames, transforms, tiers, tolerances), `geometry/` (curves,
  surfaces and the shape enum macros, now `enum_variants`),
  `operations/` (discretization, intersection, projection) and `math/`
  (the B-spline engine, the closed forms of the elementary shapes and
  the numerical toolbox, now `numeric`); `bounding_box/` stays at the
  root. Every module with children uses the `foo.rs` + `foo/` layout,
  with no `mod.rs`. `ARCHITECTURE.md` maps the new tree. Unit test
  names now carry their folder (`base::frame::tests::…`).

## [0.42.0] - 2026-10-30

UV curves, both ways, on the surface. `uv_curve_to_3d` gives the exact
3D curve a parameter-space curve traces on a surface, in the same
parameter: any curve on a plane, the affine image on a line extrusion,
the iso curve for a line along `u` or `v`, and on spline surfaces the
exact composition of any line, parabola, B-spline or Bézier UV curve.
Every other pair is refused, never approximated. `curve_3d_to_uv`, the
inverse, replaces the curves' `parametrize_on`, so curves no longer
depend on surfaces (see Breaking).

### Added

- `ParametricSurface::uv_curve_to_3d(&Curve2D, tol)`, the exact 3D
  curve a UV curve (a curve of a surface's parameter space) traces on the
  surface, in the same parameter: the inverse of `curve_3d_to_uv` (#152,
  ADR 0030). Any
  curve on a plane keeps its kind; on the extrusion of a line, the
  affine image where it keeps its kind in the same parameter; a line
  along `u` or `v` on any surface gives the iso curve, reparametrized.
  Every other pair is refused with the new `UvCurveTo3DError`
  (`NotAnalytic`, `Collapsed` at a pole or apex, or `IsoCurve`), never
  approximated. The trait's default refuses, so surface types defined
  outside the crate keep compiling.
- A straight UV curve stored as a degree-1 B-spline or Bézier curve with
  two poles is read as a line: along `u` or `v` it gives the iso curve
  run at the spline's speed, exact for a line or spline iso curve at any
  speed and for a conic at speed 1, and refused otherwise.
- `uv_curve_to_3d` composes any line, parabola, B-spline or Bézier UV
  curve on a B-spline or Bézier surface, or on the extrusion of a
  B-spline or Bézier profile, exactly into a B-spline of degree
  `k (p + q)` in the UV curve's parameter (#155, ADR 0032). An unbounded
  line or parabola is composed over its part inside the surface's
  domain, and a periodic direction is unwrapped. Circle, ellipse and
  hyperbola UV curves stay refused, since their angle parameter is not
  rational; their `to_bspline` composes exactly, in the B-spline's
  parameter. New `UvCurveTo3DError` variants:
  `OutsideSurfaceDomain` and `DegreeAboveLimit { degree }`, above
  degree 25.
- Python: `uv_curve_to_3d(uv_curve, tol=None)` on every surface
  class, taking any 2D curve, and the `UvCurveTo3DError`
  exception.
- `ParametricSurface::curve_3d_to_uv(&Curve, tol)`, the exact UV curve
  of a 3D curve lying on the surface, in the same parameter (ADR 0031).
  It takes any `Curve`, so a Bézier curve in a plane now maps too, to
  the `BezierCurve2D` over its poles' plane coordinates. Python:
  `curve_3d_to_uv(curve, tol=None)` on every surface class, taking any
  3D curve.

### Breaking

- `parametrize_on` is gone from `Line`, `Circle`, `Ellipse`, `Parabola`,
  `Hyperbola` and `BSplineCurve`; the surface maps the curve instead,
  with the same closed forms and results (ADR 0031). Migrate
  `circle.parametrize_on(surface, tol)` to
  `surface.curve_3d_to_uv(&circle.into(), tol)`, or pass a `&Curve`
  directly. `ParametrizeError` is renamed `Curve3DToUvError`, with the
  same variants, and `parametrize_numeric` is renamed
  `curve_3d_to_uv_numeric`. Python: `curve.parametrize_on(surface, tol)`
  becomes `surface.curve_3d_to_uv(curve, tol)`, the exception
  `ParametrizeError` becomes `Curve3DToUvError`, and
  `geomcore.parametrize_numeric` becomes
  `geomcore.curve_3d_to_uv_numeric`.

## [0.41.0] - 2026-10-30

Deeper modules behind the same results. One closed-form table decides
which pairs `intersect` solves exactly; transform and reversal are trait
methods on every shape, with one `TransformError`; marching is split
into a surface pair, a tracer and a seeding policy, and finds small
closed loops on every pair with an elementary or patch-convertible
side; the curve–curve results and fitting are written once for both
tiers; and the Python stubs write each method set once, with
docstrings. Marching also traces through branch crossings, stays within
confusion between its points, and works far from the origin. Rust code
calling `transformed` or the reversal methods needs
`use geomcore::prelude::*`, and B-spline `transformed` returns a
`Result` (see Breaking).

### Added

- `Transformable`, `ReversibleCurve` and `ReversibleSurface`, the
  traits of moving and reversing a shape (ADR 0024), in
  `geomcore::prelude`. Every curve and surface type implements them, and
  so do the `Curve`, `ElementaryCurve`, `Curve2D`, `Surface` and
  `ElementarySurface` enums and `TrimmedCurve` / `TrimmedCurve2D`, which
  had neither. New on the way: `transformed` on `BezierCurve`,
  `BezierSurface`, `TrimmedCurve` and the enums; `reversed` and
  `reversed_parameter` on every curve of both tiers (only the B-spline
  and Bézier curves had `reversed`, and none had `reversed_parameter`).
  A circle or ellipse reverses to `2π - u`, a line, parabola or
  hyperbola to `-u`, a B-spline or Bézier curve to `first + last - u`,
  and a trimmed curve follows its basis.
- `TransformError`, the one error of a failed transform: `Collapsed`
  for a zero scale factor, `NonFinite` when the image overflows.
- Python: every 3D curve and every surface gains `transformed`, and
  every curve gains `reversed()` and `reversed_parameter(u)`, through the
  shared method sets. A failed transform raises `ConstructionError`, as
  before.
- `GenericCurveCurveIntersection<P>` and `GenericCurveCurveHit<P>`, the
  curve–curve result and hit generic over the tier (ADR 0014, amended):
  code that reads a curve–curve result can be written once for both
  tiers.
- The Python stubs carry the docstrings of the methods every curve or
  surface class shares (evaluation, local properties, projection,
  `intersect`, ...), so an IDE shows them on hover. A test keeps each one
  equal to the docstring of the method itself.

### Changed

- B-spline first derivatives are faster (#141). `eval_d1` on a
  non-periodic curve or surface of degree 1 to 3 in every direction
  computes the basis functions and their first derivatives in one
  unrolled pass, as the point evaluation already did: about 46% faster
  on a cubic curve, 31% on a bicubic surface and 18% on rational ones.
  Results are bit-for-bit unchanged.
- One closed-form pair table decides which pairs `intersect` solves in
  closed form, the plane cuts of swept surfaces included (ADR 0023).
  Results are unchanged.
- Internal: the swept surfaces reverse and move their profile through
  the traits, so `swept.rs` keeps only the profile geometry.
- The marching intersection is split into a surface pair, a tracer and a
  seeding policy (ADR 0025). Its results are unchanged.
- Marching finds loops smaller than the spacing of its iso-curves and
  grids on more pairs (ADR 0026). Once those seeds are traced, a pair
  with an elementary surface shoots an iso-curve of the other surface
  through each maximum and minimum of the elementary surface's implicit
  equation along it, where a closed loop must enclose one; other pairs
  that convert to exact B-spline patches run the patch-pair search, as
  two B-splines did. A bump of a B-spline finer than its knot spans
  against a plane, cylinder or extrusion, and a small cap a plane cuts
  off a torus, were missed and are found; a plane against a torus finds
  every such loop in closed form. Two bounded elementary surfaces, and
  pairs that only touch, are left as they were. Marching costs at
  most 1.8 % more on the benches where nothing new is found.
- A B-spline curve against an `ExtrusionSurface` or `RevolutionSurface`
  of a bounded profile goes through the proven hull search: its
  completeness is `Complete` or `Incomplete` instead of `Sampled`.
- Internal: the property tests run over one sample list per shape enum,
  checked against the enum's variants (ADR 0027).
- `CurveCurveIntersection` and `CurveCurveIntersection2D` are now the
  aliases of `GenericCurveCurveIntersection<Point>` and
  `GenericCurveCurveIntersection<Point2D>`, and
  `CurveCurveHit`/`CurveCurveHit2D` those of `GenericCurveCurveHit`;
  the fitting constructors `interpolate`,
  `approximate` and `approximate_within` are written once on
  `GenericBSplineCurve`. Every name, method and result is unchanged.
- The Python stubs write each shared method set once, as a private
  stub-only class that the curve and surface classes inherit (ADR 0029).
  Type checkers see the same methods with the same types as before.

### Deprecated

- `RevolutionSurfaceConstructionError::NonFinite`: no longer returned;
  a failed transform reports a `TransformError`.

### Breaking

- `transformed` is a `Transformable` method and returns
  `Result<Self, TransformError>` on every type: Rust code calling it on a
  concrete type needs `use geomcore::prelude::*`, and code that matched
  the per-type construction errors (`CircleConstructionError`,
  `PlaneConstructionError`, ...) matches `TransformError::Collapsed`
  instead.
- `BSplineCurve::transformed` and `BSplineSurface::transformed` return a
  `Result` instead of the curve or surface: a zero scale, which used to
  give a B-spline collapsed onto a point, now fails with
  `TransformError::Collapsed`, and an overflow with `NonFinite`. Add
  `.unwrap()` (or `?`) where the transform is known to be a similarity
  with a non-zero scale.
- `reversed` on the B-spline and Bézier curves, and `reversed_u`,
  `reversed_v`, `u_reversed_parameter` and `v_reversed_parameter` on
  every surface and on `Surface`, moved from the types to the
  `ReversibleCurve` and `ReversibleSurface` traits, with the same
  behaviour: calls on a concrete type need `use geomcore::prelude::*`.
### Fixed

- Marching found no curves for bounded surfaces lying outside the
  model-size cube centred on the origin (#138): two tori 1e5 away
  intersected in nothing. A trace now ends on the box of the smaller
  bounded side, grown by a tenth of its diagonal, and on the model-size
  cube only when neither side is bounded (ADR 0011, amended).
- Marching lost curves where branches of the intersection cross (#139).
  At a tight confusion the steps towards a crossing shrank below the
  distance the B-spline fit separates, and the fit's refusal dropped the
  whole curve: two equal tori with perpendicular axes gave no curve at
  confusion 1e-8. Those points are now thinned before the fit. And a
  stretch between two crossings had no seed of its own, the grids
  seeding only where the surfaces are closest: a tilted torus and its
  mirror image missed a third of their intersection at any confusion.
  Where a trace stalls next to a crossing, the branches crossing there
  are now seeded and traced.
- Marched curves stayed within confusion of both surfaces only at the
  traced points (#140). Between them, the interpolant strayed by up to
  the chord deviation, about a thousand times confusion: two tori in
  general position gave curves 8e-6 off at confusion 1e-7. Where the
  surfaces cross, the fit now passes cubic Hermite segments through the
  points and their exact tangents, and each step is checked against its
  segment at the midpoint. The same tori now give curves 2e-8 off, and
  the torus–torus comparison sweep stays within confusion at 1e-7 and
  1e-9. Marching costs about 12 % more.
- Projection onto a revolution surface whose profile lies off every
  plane through the axis stopped up to 3e-9 short of the closest point
  (#142). The one-parameter minimiser behind it rejected its last
  Newton steps as rounding noise: the squared distance it minimises is
  a difference of large terms, whose rounding outgrew the minimiser's
  allowance. A full Newton step that stays convex and at least halves
  the slope is now accepted too, and those projections reach 1e-15.
  B-spline curve projection and the curve–curve searches share the
  minimiser; their tests pass unchanged.

## [0.40.1] - 2026-10-30

The swept surfaces. `ExtrusionSurface` and `RevolutionSurface` join
`Surface` with exact iso curves, conversion to B-spline patches and
projection that works on the profile, and a plane cuts either in closed
form in every common position. Nothing existing changes.

### Added

- `ExtrusionSurface` and `RevolutionSurface`, the swept surfaces
  (ADR 0016): a profile `Curve` swept along a direction,
  `profile(u) + v * direction`, or turned about an `Axis`, the profile
  point at `v` rotated by the angle `u`. Each is a `Surface` variant
  with the full `ParametricSurface` interface, exact iso curves (a line
  and the moved profile; the rotated profile and a circle), reversal,
  `transformed` and `to_bspline_patch`. A revolution of a bounded
  profile has a bounding box and converts whole with `Surface::to_bspline`;
  an extrusion is unbounded in `v`. Projection works on the profile:
  onto its image across the direction for an extrusion (in closed form
  for a line or a conic), onto the meridian for a revolution (in
  closed form for a conic in a plane through the axis, a one-parameter
  search otherwise), 4 to 80 times faster than onto the converted
  B-spline patch. Intersections run the numeric search, except the plane
  cuts below. In Python, `geomcore.surfaces.ExtrusionSurface` and
  `RevolutionSurface`.
- `intersect` cuts a swept surface with a plane in closed form, in
  either order and through the `Surface` enum: an extrusion with every
  plane (generator lines through the profile's crossings when the plane
  holds the direction, the profile carried onto the plane otherwise),
  and a revolution with a plane perpendicular to its axis (a circle
  through each crossing, a contact point on the axis) or through its
  axis when the profile lies in a plane through it (the profile turned
  into the plane). Tangent crossings give tangent curves. The result is
  complete, and exact unless a spline profile's crossings came from the
  numeric curve-plane search; other positions still march. A cut takes
  under 5 µs where marching the converted patch takes 0.5 to 6 ms.

## [0.40.0] - 2026-10-30

An architecture release. Shared behaviour now lives in one place per
concept: on the traits, which `geomcore::prelude` brings into scope; in
one variant list per shape enum, which also gives Bézier curves and
surfaces their own variants; in one generic type per curve across space
and the plane; and in one step from a typed intersection result to a
family result. Every public Rust item now has one path, at the crate
root. Rust code needs a few mechanical edits, listed under Breaking.
Python only gains methods, and its `intersect` now always releases the
interpreter lock.

### Added

- `geomcore::prelude`, which re-exports the shape and operation traits
  (`ParametricCurve`, `ParametricCurve2D`, `PointExtrema`,
  `PointExtrema2D`, `ParametricSurface`, `HasBoundingBox`,
  `HasBoundingBox2D`, `Discretize`, `Discretize2D`, `Intersect`,
  `IntersectExact`). `use geomcore::prelude::*` brings every shared
  method into scope (ADR 0015).
- `project_point`, `project_points`, `contains` and `parameter_of` are
  methods of `PointExtrema` and `PointExtrema2D`, provided from
  `extrema`, so every curve answers them, `TrimmedCurve` and a curve of
  your own included. `ParametricSurface` gains `project_points`,
  `contains` and `parameters_of`. The elementary types keep their closed
  forms.
- `eval_derivative` and `eval_derivative_sided` of B-spline and Bézier
  curves evaluate any order: zero above the degree for a non-rational
  curve, every order of a rational one through the quotient rule. Orders
  1 and 2 keep their unrolled path. In Python, orders above 2 no longer
  raise `ValueError`.
- In Python, `parameter_of` and `parameters_of` take an optional `tol`
  (`Tolerance.default()` when omitted).
- `Curve::Bezier`, `Curve2D::Bezier` and `Surface::Bezier`: a Bézier
  curve or surface is a shape enum variant of its own, boxed,
  discretized, converted and parametrized as the B-spline it is
  (ADR 0017).
- `BezierCurve`, `BezierCurve2D` and `BezierSurface` are operands of
  `intersect` in both orders, against every curve and surface (and
  every 2D curve), whole or trimmed, routed to the numeric search like
  their B-spline. In Python, every function that takes a curve or
  surface accepts the Bézier classes, `TrimmedCurve` and
  `TrimmedCurve2D` included.
- `BezierSurface::u_iso`, `v_iso`, `u_reversed_parameter` and
  `v_reversed_parameter`, so `Surface::u_iso` of a Bézier surface
  returns `Curve::Bezier`.
- In Python, every class has every method its core type has. The Bézier
  classes gain `intersect`; `BezierCurve` and `BezierCurve2D` gain
  `tangent`, `curvature`, `curvature_center`, `period` and
  `is_periodic`; `BezierSurface` gains `u_iso`, `v_iso`,
  `u_reversed_parameter`, `v_reversed_parameter`, `continuity_u`,
  `continuity_v`, `min_continuity_u` and `min_continuity_v`.
  `BSplineCurve`, `TrimmedCurve` and `TrimmedCurve2D` gain
  `parameter_of`; `BSplineSurface` and `BezierSurface` gain
  `parameters_of`.
- `TierPoint` and the tier-generic curve types `GenericLine`,
  `GenericCircle`, `GenericEllipse`, `GenericParabola`,
  `GenericHyperbola`, `GenericBSplineCurve` and `GenericBezierCurve`.
  Each curve type is now one generic type over its tier, and `Circle`
  and `Circle2D` (and the other pairs) are its two aliases, so a
  function over `GenericCircle<P>` with `P: TierPoint` serves both tiers
  (ADR 0014). `CurveProjection` is generic the same way, defaulting to
  space, and `CurveProjection2D` is its planar alias.
- `TrimmedCurve::into_curve`, the new name of `into_curve3d`.
- `TrimmedCurve::to_bspline` for a B-spline, Bézier or `Curve` basis,
  as `TrimmedCurve2D` already had: the piece is extracted from the
  B-spline. `TrimmedCurve<ElementaryCurve>` implements `Discretize`.

### Fixed

- In Python, a 2D intersection with a B-spline or Bézier curve now
  releases the interpreter lock while it searches, as the 3D ones
  already did, so other Python threads keep running.

### Changed

- Each shared method is documented once, on its trait, instead of on
  every type. `ARCHITECTURE.md` placement rule 2 and ADR 0015 state the
  rule: a method that more than one type answers is a trait method, and
  a type's inherent `impl` keeps only what is specific to it.
- ADR 0012 now points to `parameter_of(point, tol)` and
  `parameters_of(point, tol)`.
- The `Curve`, `ElementaryCurve`, `Curve2D`, `Surface` and
  `ElementarySurface` enums forward every trait method to the shape they
  hold, provided methods included (`tangent`, `curvature`,
  `is_periodic`, `project_points`, `normal`, …), so a shape's own
  override is what the enum answers. Internally each enum lists its
  variants once, and every forwarding match and `From` impl is
  generated from that list (ADR 0017).
- In Python, `intersect` releases the interpreter lock around every
  call, not only around the pairs the bindings expected to search
  (ADR 0021). A closed-form call from Python costs about 47 ns more.
- In the Python bindings, each method set a family of classes shares
  (`curve`, `curve_2d`, `surface`, `discretize`, ...) is written once in
  `geomcore-py/src/method_sets.rs`, and the classes are listed once in a
  registry in `arguments.rs`, which gives the conversions to and from the
  core enums and the `TypeError` text. A shared method's docstring is
  now the same on every class; the class docstring says what its
  parameter is. The intersection results share their common members the
  same way.
- The router reads each typed exact result straight into its family
  result, in one step: no intermediate untyped result, and one bridge
  instead of four. Results are unchanged.
- The 3D and 2D copies of every curve type are merged into one generic
  type per curve: the shared methods, the trait impls of both tiers, the
  B-spline conversion and the discretization are written once. The
  public names, methods and behaviour stay, so callers compile
  unchanged; rustdoc lists the curves as the generic structs, with
  `Circle`, `Circle2D` and the others as type aliases (ADR 0014,
  amended).
- `Circle2D` projection reports parameter `0` for a point within
  `tol.confusion` of the seam, as `Circle` and `Ellipse2D` already did,
  instead of an angle just below `2π`.
- Internally, the five range traits behind trimmed curves (box of a
  range in each tier, discretization of a range, B-spline of a range in
  each tier) are one crate-private trait, generic over the tier, and the
  whole bounded curves derive their box from it over their bounds. The
  box formulas of the two tiers are written once. Results are unchanged
  (ADR 0005 and ADR 0014, amended). In Python, `TrimmedCurve.to_bspline`
  calls the core method directly.
- The documentation lists every item once, at the crate root, and the
  prelude's traits link there instead of repeating their pages. The
  folder overviews of intersections, bounding boxes and discretization
  moved onto the crate page and the `IntersectExact`, `HasBoundingBox`
  and `Discretize` traits.
- `ParametricCurve::as_bspline`, `ParametricCurve::as_elementary`,
  `ParametricCurve2D::as_bspline`, `ParametricSurface::as_bspline` and
  `ParametricSurface::as_elementary` are hidden from the documentation:
  they are routing hooks of the numeric search, with no stability
  promise. They stay on the traits, because the search is generic over
  any curve and surface (ADR 0022). A wrapper around a curve keeps
  forwarding them.
- Internal: the numeric searches on B-spline pieces (curve–curve
  intersection and extrema, curve–surface intersection, and the patch
  pairs that seed the surface–surface marching) run through one shared
  subdivision driver. Results, completeness and speed are unchanged.
- Internal: `CurveCurveRegion`, `CurveSurfaceRegion` and
  `CurveCurveExtrema` are defined beside the other intersection results.
  Their public paths are unchanged.

### Deprecated

- `TrimmedCurve::into_curve3d`: use `into_curve` (ADR 0013). The old
  name stays for one release.

### Removed

- The inherent twins of trait methods, on every curve and surface type
  and on the `Curve`, `ElementaryCurve`, `Curve2D`, `Surface` and
  `ElementarySurface` enums: `eval_point`, `eval_points`,
  `eval_derivative`, `project_point`, `project_points`, `contains` and
  `as_elementary` on both; `eval_d1`, `eval_d2`, `bounds`, `is_closed`,
  `is_periodic`, `continuity`, `min_continuity`, `extrema` and
  `parameter_of` on curves; `parameters_of`, `eval_d1`, `eval_d2`,
  `u_bounds`, `v_bounds`, `is_u_closed`, `is_v_closed`, `continuity_u`,
  `continuity_v`, `min_continuity_u` and `min_continuity_v` on surfaces.
  The trait methods of the same names remain, with the same behaviour.
  `BSplineSurface::extrema` and `BezierSurface::extrema` stay inherent.

### Breaking

- Rust code that calls a shared method on a concrete type without
  importing its trait no longer compiles: `circle.eval_point(u)` needs
  `use geomcore::prelude::*` (or the trait itself) in scope. Generic code
  and calls through the traits are unaffected.
- `parameter_of` and `parameters_of` take a tolerance:
  `curve.parameter_of(point, tol)`, `surface.parameters_of(point, tol)`.
  The elementary types still invert in closed form and ignore it; a
  B-spline projects the point within it instead of within
  `Tolerance::DEFAULT`. Pass `Tolerance::DEFAULT` for the old behaviour.
- The untyped exact results are gone (ADR 0018):
  `ExactSurfaceIntersection`, `ExactCurveSurfaceIntersection` and
  `ExactCurveCurveIntersection`, with the `IntersectExact` impls that
  returned them: `ElementarySurface` with `ElementarySurface`,
  `ElementaryCurve` with `ElementaryCurve` and with `ElementarySurface`
  (both orders), `Surface` with `Surface`, `Curve` with `Curve`, and
  `Curve` with `Surface` (both orders). `IntersectExact` itself stays,
  for pairs of concrete types. Migration: call `intersect` on the same
  operands and read the family result. `NotAnalytic` becomes
  `method() == IntersectionMethod::Numeric` (the numeric search now runs
  instead of returning nothing), `Coincident` becomes non-empty
  `overlaps()` for curves or `is_coincident()` for surfaces, `Points` /
  `Meets` become `hits()` or `curves()` and `contact_points()`, and
  `Empty` becomes `is_empty()`. A curve–surface hit's `parameter` is
  `curve_parameter`.
- Every public item has one path, at the crate root (ADR 0022). The
  folder modules `geomcore::curves`, `geomcore::surfaces`,
  `geomcore::intersect`, `geomcore::bounding_box`,
  `geomcore::discretization`, `geomcore::frame`, `geomcore::point`,
  `geomcore::vector`, `geomcore::transform`, `geomcore::tolerance`,
  `geomcore::projection` and `geomcore::tier` are private. Migration:
  replace `geomcore::curves::Line` with `geomcore::Line`, and so on;
  every item kept its name. Rust has no deprecation for a module path,
  so there is no transition release. `geomcore::prelude` stays. No item
  was made private: the audit of the 146 public items found no other
  internal one at the root. Python is unaffected.

## [0.39.1] - 2026-10-29

A layout release with no new behaviour: the source tree, in Rust and in
the Python bindings, is reorganized so that each concept has one
predictable place. Nothing changes for users of the crate or the Python
package.

### Changed

- The source tree is reorganized so that a newcomer can find things by
  folder, with no change to the public API. `ARCHITECTURE.md` maps every
  folder and states the placement rules.
  - `intersect/` is split by method into `exact/`, `numeric/` and
    `results/`. The exact solvers have left the shape files and sit in
    `intersect/exact/`, grouped by pair family.
  - Every public type has its own file, the 2D types included.
  - The shape traits sit in `parametric_curve.rs`,
    `parametric_curve_2d.rs` and `parametric_surface.rs`, apart from the
    `Curve`, `Curve2D` and `Surface` enums.
  - Unit tests have moved out of the source files into sibling `tests.rs`
    files.
  - The internal engines are grouped by concept. `bspline/` holds the
    B-spline engine of curves and surfaces, and `elementary/` holds the
    closed-form formulas. They replace `curve_math/` and `surface_math/`.
  - The Python bindings use the core crate's folder names. Each class
    has its own file under `geomcore-py/src/curves/` or `surfaces/`, and
    `numeric.rs` is split between `arguments.rs` and `intersect.rs`.
- The intersection-matrix Python example moved to `geomcore-py/examples/`.
  The scene it writes is no longer committed.

## [0.39.0] - 2026-10-29

A cleanup release with no new behaviour: the names deprecated in 0.36
and 0.37 are gone. Code that builds on 0.38 without deprecation warnings
builds on 0.39 unchanged, unless it names one of the old traits
(`ParametricCurve3D`, `PointExtrema3D`, `Bound3D`, `Discretize3D`),
which were re-exported without a warning.

### Removed

- The 0.36 names deprecated in 0.37 (breaking, ADR 0013): `Point3D`,
  `Line3D`, `Curve3D`, `ParametricCurve3D` and the other `*3D` type
  aliases and trait re-exports in Rust, and the matching names under
  `geomcore` and `geomcore.curves` in Python, which raised a
  `DeprecationWarning`. Use the unsuffixed names in the 0.37 migration
  table.
- `Vector::angle_with_ref`, deprecated in 0.36 (breaking): use
  `angle_with_reference`.

## [0.38.0] - 2026-10-29

The 2D curve tier now matches the 3D one (ADR 0014). `Ellipse2D`,
`Parabola2D`, `Hyperbola2D` and `BSplineCurve2D` join `Line2D` and
`Circle2D`. A conic or B-spline on a plane maps to its exact image in
the plane, and every pair of 2D curves intersects, exactly where a
closed form exists. Bézier curves and patches are first-class types.
This release breaks compatibility in two places: `to_bezier_pieces` /
`to_bezier_patches` return the new Bézier types, and `Curve2D` is no
longer `Copy` (see Changed).

### Added

- `BSplineCurve2D` (#90), in Rust and Python: construction (polynomial,
  rational, interpolation, approximation), evaluation and sided
  evaluation, continuity, knot insertion and refinement, degree
  elevation, knot removal, Bézier pieces, splitting, extraction,
  reversal, length, projection, discretization and intersection.
  `Curve2D` gains a `BSpline` variant, and `Circle2D`,
  `TrimmedCurve2D` and `Curve2D` gain `to_bspline`.
- `ParametricCurve2D` gains the local properties and queries of its 3D
  counterpart (#91), in Rust and Python: `tangent`, `curvature`,
  `curvature_center`, `is_closed`, `continuity`, `min_continuity` and,
  in Rust, `as_bspline`. `TrimmedCurve2D` gains `is_closed`,
  `continuity` and `min_continuity`, and `Curve2D` forwards them all.
  The planar curvature is **signed**, `(x'y'' − y'x'') / |C'|³`: positive
  where the tangent turns left, negative where it turns right (3D
  curvature stays unsigned). `curvature_center` follows the turn, so a
  `Circle2D` centres at its own centre whichever way it runs.
- `BoundingBox2D` and `HasBoundingBox2D` (#91), in Rust and Python: the
  planar counterpart of `BoundingBox`, tight for circles and arcs, from
  the control points for B-splines, `None` for lines, and always present
  on a `TrimmedCurve2D`.
- 2D curve–curve intersections with a B-spline operand run the numeric
  search, each curve seen at `z = 0`. The result reports
  `IntersectionMethod::Numeric` and is proved complete where the 3D
  search would prove it. `CurveCurveIntersection2D::method` now reports
  this method instead of always `Exact`.
- `Ellipse2D` (#92), in Rust and Python: construction from a center,
  direction and radii, from a `Frame2D`, or from a center and two
  points; evaluation and derivatives, parameter inversion, containment,
  extrema and projection; area, tight bounding boxes and proven
  discretization; exact B-spline conversion; and `TrimmedCurve2D`
  support. `Curve2D` gains an `Ellipse` variant with every operation
  forwarded, and Python gains `geomcore.curves.Ellipse2D` with stubs.
- `Parabola2D` (#93), in Rust and Python: construction from an apex,
  axis direction and focal distance or from a `Frame2D`, sharing
  `ParabolaConstructionError` with `Parabola`; evaluation and
  derivatives with the same parametrization as `Parabola`, parameter
  inversion, containment, extrema and projection; tight bounding boxes
  and proven discretization of its pieces; and `TrimmedCurve2D`
  support. `Curve2D` gains a `Parabola` variant with every operation
  forwarded, and Python gains `geomcore.curves.Parabola2D` with stubs.
  A trimmed parabola converts exactly to a B-spline (`TrimmedCurve2D::to_bspline`).
- `Hyperbola2D` (#94), in Rust and Python: construction from a center,
  major-axis direction and radii, from a `Frame2D` or from a center and
  two points, sharing `HyperbolaConstructionError` with `Hyperbola`;
  evaluation and derivatives with the same parametrization as
  `Hyperbola`, parameter inversion, containment, extrema and
  projection; tight bounding boxes, proven discretization and exact
  rational B-spline conversion of its pieces; and `TrimmedCurve2D`
  support. `Curve2D` gains a `Hyperbola` variant with every operation
  forwarded, and Python gains `geomcore.curves.Hyperbola2D` with stubs.
- Exact plane images for conics and B-splines (#95), in Rust and Python:
  `Ellipse::parametrize_on`, `Parabola::parametrize_on`,
  `Hyperbola::parametrize_on` and `BSplineCurve::parametrize_on` on a
  plane return an `Ellipse2D`, `Parabola2D`, `Hyperbola2D` or
  `BSplineCurve2D` in the plane's `(u, v)` (they always returned
  `NotAnalytic` before), with the 3D curve's parametrization:
  `plane(q(t)) == curve(t)`. A B-spline's image has the poles' plane
  coordinates with degree, knots, multiplicities, weights and periodicity
  unchanged. A curve that leaves the plane reports `CurveNotOnSurface`: a
  conic whose normal turns off the plane's by more than
  `tol.angular_for_extent(extent)` (its major radius for an ellipse, the
  model size for an unbounded parabola or hyperbola) or that sits off it, a
  B-spline with a pole farther than `tol.confusion` from the plane. The
  Python stubs return the matching 2D type.
  The other pairs were checked and stay `NotAnalytic`: an ellipse,
  parabola or hyperbola on a cylinder, cone, sphere or torus (a plane
  section of a cylinder is a sinusoid in `(u, v)`, of a cone a polar
  conic, neither a conic or a B-spline), a B-spline curve on those
  surfaces, and every curve on a B-spline surface. Golden fixtures in
  `parametrize.json` cover the plane cases, the off-plane cases and the
  curved-surface cases.
- Exact 2D intersections for every pair of lines and conics (#96), in
  Rust and Python: `Ellipse2D`, `Parabola2D` and `Hyperbola2D` gain
  `Intersect` against `Line2D`, `Circle2D`, each other, `Curve2D` and
  `TrimmedCurve2D`, in either order. Through `Curve2D` and Python these pairs
  used to run the numeric search.
  The pair is solved in closed form by the coplanar solvers of the 3D
  pair (the line-conic substitution and the conic-conic solver, up to
  four hits), the curves seen at `z = 0` with their own parameters, so
  a tangency reports multiplicity 2 and coincident pieces report
  `CurveCurveOverlap`s. The result is `IntersectionMethod::Exact` and
  complete. `IntersectExact` gives the typed result
  `Conic2DIntersection` (`Hits`, `Coincident`, `Empty`) for every pair
  of lines and conics other than line-line, line-circle and
  circle-circle, which keep their own results. `CurveCurveIntersection2D`
  holds four hits inline. A pair with a `BSplineCurve2D` is still
  searched numerically and reports its `Completeness` like the 3D
  search. Golden fixtures in `intersect_2d.json` cover every pair
  family, replayed in both orders, and the `intersect_curve_curve`
  bench times the 2D pairs.
- `BezierCurve` and `BezierCurve2D` (#98): first-class Bézier curves of
  `degree + 1` poles (polynomial or rational) with evaluation, sided
  evaluation, continuity, `split_at`, `elevate_degree`, `reversed`,
  `to_bspline`, `from_bspline`, length, projection, bounding box and
  discretization, and `BezierConstructionError`.
- `BezierSurface` (#98): the first-class Bézier patch, a grid of
  `(u_degree + 1) x (v_degree + 1)` poles (polynomial or rational) with
  evaluation to second order, sided evaluation, projection, `split_u_at`
  / `split_v_at`, `elevate_u_degree` / `elevate_v_degree`, `reversed_u`
  / `reversed_v`, `to_bspline`, `from_bspline`, bounding box and the
  `ParametricSurface` trait. It reuses `BezierConstructionError`.
- Intersections and extrema involving a `BezierCurve` or `BezierSurface`
  take the complete Bézier-piece search of the B-spline they are, as a
  `BSplineCurve` or `BSplineSurface` does, instead of the sampled search.
- `BezierCurve`, `BezierCurve2D` and `BezierSurface` in Python, with
  stubs (#98).
- Golden fixtures (`bezier.json`) check Bézier curve and surface
  evaluation and derivatives against explicit Bernstein sums to 1e-7.

### Changed

- `BSplineCurve::to_bezier_pieces` and `BSplineCurve2D::to_bezier_pieces`
  return `BezierCurve` and `BezierCurve2D` instead of single-span
  B-splines (#98), in Rust and Python; `to_bspline` gives the B-spline
  back.
- `BSplineSurface::to_bezier_patches` returns `Vec<Vec<BezierSurface>>`
  (a list of lists of `BezierSurface` in Python) instead of single-span
  B-spline surfaces (#98); `to_bspline` gives the B-spline surface back.
- `BSplineCurve` and `BSplineCurve2D` share one implementation: a
  private spline generic over its pole type (ADR 0014). The two types
  agree on every operation to 1e-12, which a new test checks over the
  B-spline fixtures.
- `BSplineCurve::project_point`, `extrema` and `contains` handle a
  query point with a non-finite coordinate as the other curves do: a
  NaN projection, no extrema, and not contained.
- `Curve2D` and `TrimmedCurve2D<Curve2D>` are no longer `Copy`: the
  enum now holds a `BSplineCurve2D`, whose poles live on the heap.
  Callers that copied one implicitly `.clone()` it.

### Fixed

- `Ellipse::from_center_and_points` (and the new
  `Ellipse2D::from_center_and_points`) returns `NonFinite` for a NaN or
  infinite input point instead of `InvertedAxis`.
- `ElementaryCurve::to_bspline` and `Curve::to_bspline` (and the new
  `Curve2D::to_bspline`) return `NotRepresentable` for a circle or
  ellipse whose poles overflow `f64` instead of panicking.
  `Circle::to_bspline` and `Ellipse::to_bspline` still panic there, and
  now document it.
- `Parabola::contains`, `extrema` and `project_point` stay correct for a
  focal distance large enough (about 1e154) that the squared gradient or
  the stationarity cubic used to overflow, which gave a wrong
  containment answer and a NaN projection.
- `BSplineSurface::project_point` now projects points inside a very narrow
  knot span correctly. Its Newton refinement minimises in each subpatch's
  own `[0, 1]` coordinates, so a 1e-9-wide span beside a full-width one no
  longer damps the wide direction or stops before the point is reached
  (#80).
- The sampled curve–surface search seeds an unbounded curve from the
  surface's bounding box, not from a fixed `[-10, 10]` parameter window.
  When the surface is a sphere, torus or B-spline (whose box is known)
  and the curve is elementary (whose reach is), the seeds cover the part
  of the curve that can reach that box, so a hit far along the curve is
  searched where the surface is instead of where the parameter happens to
  be. A parabola crossing a torus centred at its `t = 15` point now
  reports all four crossings; the old window found two. Curves and
  surfaces whose types do not report a box keep the old window, and the
  result stays `Sampled` (#74).
- The sampled curve–curve search no longer drops a transversal crossing
  that lies close to another crossing or to a near touch (#13). It runs
  in `intersect_curve_curve` and `curve_curve_extrema` (Rust and Python)
  on every pair the B-spline piece search does not take (two elementary
  curves, curves defined outside the crate), and completes a piece search
  that ends `Incomplete`, including under the routed `intersect`.
  It seeded the minimiser from the discrete local minima of a 32 × 32
  parameter grid, which sees one minimum where two lie closer together
  than the grid spacing; the minimiser then reached only one of them.
  The pairs of grid spans whose chords come within the curves' offset
  from them are now halved until each holds a single minimum, and each
  is seeded. On random coplanar circle–ellipse pairs, 184 of 20,000 in
  the XY plane and 51 of 5,000 in random frames missed a crossing; none
  do now, nor any circle–ellipse, ellipse–ellipse or line–conic pair of
  the new randomized comparison against the exact solvers. The search
  costs more where the curves meet, as each crossing is now refined: a
  circle–ellipse pair through `intersect_curve_curve` takes about 20 µs
  instead of 14 µs. Sampled extrema between curves that stay apart cost
  up to about 7 % more, near the benchmarks' noise.

## [0.37.0] - 2026-10-19

One naming release for Rust and Python, with no change in behaviour: 3D
is the unmarked dimension, so only 2D types carry a dimension suffix
(ADR 0013). `Point3D` becomes `Point`, `Line3D` becomes `Line`, and so
on, the same convention the results (`CurveProjection` /
`CurveProjection2D`) and the surfaces already followed. The old names
still work for this release and warn on use.

### Migration

Rust. The old type names are deprecated type aliases. The old trait
names are hidden re-exports, so they still compile without a warning,
and they are removed together with the aliases in the next breaking
release.

| 0.36 | 0.37 |
|---|---|
| `Point3D`, `Vector3D`, `Axis3D`, `Frame3D` | `Point`, `Vector`, `Axis`, `Frame` |
| `Line3D`, `Circle3D`, `Ellipse3D`, `Parabola3D`, `Hyperbola3D`, `BSplineCurve3D` | `Line`, `Circle`, `Ellipse`, `Parabola`, `Hyperbola`, `BSplineCurve` |
| `Curve3D`, `TrimmedCurve3D`, `ParametricCurve3D`, `PointExtrema3D` | `Curve`, `TrimmedCurve`, `ParametricCurve`, `PointExtrema` |
| `BoundingBox3D`, `Bound3D` | `BoundingBox`, `HasBoundingBox` |
| `Discretization3D`, `Discretize3D` | `Discretization`, `Discretize` |
| `LineLine3DIntersection`, `LineCircle3DIntersection`, `CircleCircle3DIntersection`, `LineConic3DIntersection` | `LineLineIntersection`, `LineCircleIntersection`, `CircleCircleIntersection`, `LineConicIntersection` |

Python. The old names resolve to the new classes and raise a
`DeprecationWarning`.

| 0.36 | 0.37 |
|---|---|
| `geomcore.Point3D`, `Vector3D`, `Axis3D`, `Frame3D`, `BoundingBox3D`, `Discretization3D` | `geomcore.Point`, `Vector`, `Axis`, `Frame`, `BoundingBox`, `Discretization` |
| `geomcore.curves.Line3D`, `Circle3D`, `Ellipse3D`, `Parabola3D`, `Hyperbola3D`, `BSplineCurve3D`, `TrimmedCurve3D` | `geomcore.curves.Line`, `Circle`, `Ellipse`, `Parabola`, `Hyperbola`, `BSplineCurve`, `TrimmedCurve` |

The 2D names (`Point2D`, `Line2D`, `Curve2D`, …) and every surface name
are unchanged.

### Added

- ADR 0014: the 2D curve tier is to mirror the 3D one (ellipse,
  parabola, hyperbola and B-spline curves in 2D, so curve-on-surface
  parametrization can return the exact image of every curve).

## [0.36.0] - 2026-10-19

One breaking release for Rust and Python. `a.intersect(&b, tol)` now
answers every supported pair of curves and surfaces: exactly when a
closed form exists, numerically otherwise, and the result says which
method ran and whether every hit was found (ADR 0011, superseding ADR
0001). Results describe themselves: every hit has its point, every
projection its closest point, and Python results are objects instead
of tuples (ADR 0012). Curves and surfaces gain their local
properties (unit tangent, curvature and centre of curvature, unit
normal, principal curvatures), closedness and continuity queries,
surface reversal and one-sided B-spline evaluation, in Rust and Python.
Mirrored elementary shapes now follow the mirror, normal included. The
migration tables below list every renamed or removed call and every
changed behaviour.

### Migration

Rust:

| 0.35 | 0.36 |
|---|---|
| `a.intersect(&b, tol)` returning the pair's exact enum (`PlaneSphereIntersection`, …) | `a.intersect_exact(&b, tol)` (`use geomcore::IntersectExact`) |
| `a.intersect(&b, tol)` on `Surface`, `Curve3D` or the elementary enums, returning `SurfaceIntersection`, `CurveSurfaceIntersection`, `CurveCurveIntersection` | `a.intersect_exact(&b, tol)`, returning `ExactSurfaceIntersection`, `ExactCurveSurfaceIntersection`, `ExactCurveCurveIntersection`; or `a.intersect(&b, tol)` for the routed result |
| `NotAnalytic`, then `marching_intersection(&a, &b, tol)` or `surface.intersect_marching(&other, tol)` | `a.intersect(&b, tol).curves()` (a `SurfaceSurfaceIntersection`; `method()` says `Numeric`) |
| `intersect_curve_surface(&c, &s, tol)` returning `CurveSurfaceHits` | `c.intersect(&s, tol)` for the crate's types; `intersect_curve_surface` (now returning `CurveSurfaceIntersection`) for types defined outside the crate; read `.hits()` and `.completeness()` |
| `intersect_curve_curve(&a, &b, tol)` returning `CurveCurveHits` | `a.intersect(&b, tol)`, or `intersect_curve_curve` returning `CurveCurveIntersection` |
| per-pair methods deprecated in 0.29 (`line.intersect_sphere(&sphere, tol)`, …) | `line.intersect_exact(&sphere, tol)` |
| `hit.first_param`, `hit.second_param`, `hit.curve_param`, `hit.surface_params` | `hit.first_parameter`, `hit.second_parameter`, `hit.curve_parameter`, `hit.surface_parameters`, each an `Option`; `None` → `curve.parameter_of(hit.point)` / `surface.parameters_of(hit.point)` |
| `LineTorusHit::line_param`, `surface_params`; `LineConicHit::line_param`, `conic_param`; `CurveCurveExtremum::first_param`, `second_param` | `line_parameter`, `surface_parameters`; `line_parameter`, `conic_parameter`; `first_parameter`, `second_parameter` (plain `f64`) |
| `projection.parameter` / `(u, v)` then `curve.eval_point(…)` for the closest point | `projection.point` |
| `project_point`, `project_points` and `PointExtrema2D::extrema` on 2D curves returning `CurveProjection` | the same calls, returning `CurveProjection2D` (with a `Point2D` `point`) |
| `v.angle_with_ref(other, reference)` (deprecated) | `v.angle_with_reference(other, reference)` |

Python:

| 0.35 | 0.36 |
|---|---|
| `kind, payload = a.intersect(b)` | `result = a.intersect(b)`: a `CurveCurveIntersection`, `CurveCurveIntersection2D`, `CurveSurfaceIntersection` or `SurfaceSurfaceIntersection` with `hits`, `overlaps`, `curves`, `contact_points`, `method`, `completeness` |
| `a.intersect_sphere(b)`, `a.intersect_line(b)`, … | `a.intersect(b)` |
| `geomcore.intersect_marching(a, b)` | `[piece.curve for piece in a.intersect(b).curves]` |
| `geomcore.intersect_curve_surface(curve, surface)` → `[(t, (u, v), point)]` | `curve.intersect(surface).hits` → `[CurveSurfaceHit]` with `curve_parameter`, `surface_parameters`, `point`, `multiplicity` |
| `geomcore.intersect_curve_curve(a, b)` → `[(t1, point, t2)]` | `a.intersect(b).hits` → `[CurveCurveHit]` with `first_parameter`, `second_parameter`, `point`, `multiplicity` |
| `curve_curve_extrema(a, b)` → `[(t1, t2, distance)]` | `curve_curve_extrema(a, b).minima` → `[CurveCurveExtremum]` with `first_parameter`, `second_parameter`, `distance` |
| `curve_curve_extrema(curve_a=a, curve_b=b)` | `curve_curve_extrema(first=a, second=b)` |
| `parametrize_numeric(curve, surface, n=64)` | `parametrize_numeric(curve, surface, sample_count=64)` |
| `parameter, distance = curve.project_point(p)` | `foot = curve.project_point(p)`: `foot.parameter`, `foot.point`, `foot.distance` |
| `u, v, distance = surface.project_point(p)` | `foot = surface.project_point(p)`: `foot.u`, `foot.v`, `foot.point`, `foot.distance` |
| `curve.extrema(p)` → `[(parameter, distance)]` | `curve.extrema(p)` → `[CurveProjection]` |
| `Vector3D.x()`, `Vector3D.y()`, `Vector3D.z()` | `Vector3D.unit_x()`, `unit_y()`, `unit_z()` (`v.x` is now the component) |
| `except ValueError` | still works; `except geomcore.ConstructionError` (and the other subclasses of `GeomcoreError`) to tell errors apart |

Behaviour, in Rust and Python alike:

| 0.35 | 0.36 |
|---|---|
| `s.transformed(t)` on a plane, cylinder, cone, sphere, torus, circle, ellipse, parabola or hyperbola re-derived a direct frame, so a mirror kept the normal on its side | every frame axis maps through `t`, so `s.transformed(t)` evaluated at `(u, v)` is `t` applied to `s` there (as B-splines already did); a mirror turns the frame indirect and the normal to the opposite of the mirrored normal. Call `reversed_u()` or `reversed_v()` on the image to get the old side back |
| `Plane::normal` and `Circle3D::normal` (Python `Plane.normal()`, `Circle3D.normal()`): the frame's z direction | the z direction on a direct frame, its negation on an indirect one (a mirrored or reversed shape). Unchanged for every plane and circle built by a constructor |
| every frame right-handed (`x × y = z`) | `Frame3D::from_axes`, reversal and mirrors make left-handed frames; test `Frame3D::is_direct()` instead of assuming `x × y = z` |
| marching judged a surface normal degenerate when `\|Su × Sv\| ≤ 1e-9 (\|Su\|² + \|Sv\|²)` | a derivative counts as zero when the distance it would cover across its parameter range is within `tol.confusion` (ADR 0010), the rule `normal`, `curvatures`, `tangent`, `curvature` and `curvature_center` also follow, so rescaling a shape's parameters changes none of them |
| no unit normal: callers normalized `Su × Sv`, which vanishes at a sphere's pole and a cone's apex | `normal(u, v, tol)` returns an `Option` (Python: `None` or a `Vector3D`): the limiting normal at a pole, `None` at a cone's apex (each meridian has its own) and wherever no normal can be recovered; `curvatures` is `None` at both |

### Added

- `Intersect<B>::intersect` routes (ADR 0011): for every pair of curves
  and surfaces, in either order, concrete types, the `Curve3D`,
  `Curve2D`, `ElementaryCurve`, `Surface` and `ElementarySurface` enums,
  B-splines and trimmed curves (`TrimmedCurve3D`, `TrimmedCurve2D`), it
  solves the pair in closed form when one exists and runs the numeric
  search otherwise: the Bézier-hull search for curve–curve and
  curve–surface, marching for surface–surface. A pair of concrete types
  calls its typed solver directly, with no run-time dispatch. Swapping
  the arguments swaps the parameters in each hit.
- Family results, structs with accessors: `CurveCurveIntersection`,
  `CurveCurveIntersection2D`, `CurveSurfaceIntersection` and
  `SurfaceSurfaceIntersection`. Each has `method()`
  (`IntersectionMethod::Exact` or `Numeric`), `completeness()`,
  `is_complete()` and `is_empty()`; the curve results have `hits()` and
  `overlaps()`, the surface result `curves()`, `contact_points()` and
  `is_coincident()`. An exact result is always complete. The curve results store up to
  four hits inline (two in 2D) and spill to the heap beyond; the exact
  line–torus and conic–surface solvers still allocate.
- `CurveCurveHit2D`, `CurveCurveOverlap` (a parameter range on each
  curve plus `same_direction`; a coincidence that crosses either curve's
  seam is split there, so each range increases in its own curve's
  parameter), `IntersectionCurve` (a `Curve3D` plus a `tangent` flag)
  and `SurfaceSurfaceRegion`. `CurveSurfaceIntersection::overlaps()`
  gives the curve's parameter ranges that lie on the surface.
- The numeric curve–curve and curve–surface search reports overlaps:
  where a curve lies on the other operand (a planar B-spline in a plane,
  a straight B-spline on a line, two B-splines along one line), the hits
  and regions left along it are joined into an overlap instead of a
  string of hits, its ends found within `tol.confusion`; the result
  stays `Incomplete`, with an `Unseparated` region over the overlap.
- `CurveSurfaceHit::multiplicity` (breaking for code that builds the
  struct with a literal): 1 where the curve crosses the surface, 2 or
  more where it is tangent.
- A trimmed curve is intersected as its basis, then its range applies:
  a hit is kept when it lies within `tol.confusion` in space of the
  range (the parameter excess beyond an end times the curve's speed
  there) and reports that end as its parameter; this holds at every
  end, the seam of a circle or ellipse included. Overlaps are clipped
  to the range, and a clipped piece no longer than `tol.confusion` in
  space is reported as a hit where the curves touch.
- The closest point on every projection: `CurveProjection::point` and
  `SurfaceProjection::point` (breaking for code that builds these
  structs with a literal), and `CurveProjection2D` (with a `Point2D`)
  for 2D curves.
- `project_point` and `project_points` as inherent methods on
  `Curve3D`, `ElementaryCurve`, `Curve2D`, `Surface` and
  `ElementarySurface`, so no trait import is needed (`Surface` and
  `ElementarySurface` had them through `ParametricSurface` before; the
  curve enums did not); `PointExtrema3D` for `ElementaryCurve`.
- `parameter_of` on `Curve3D`, `Curve2D` and `ElementaryCurve`, and
  `parameters_of` on `Surface` and `ElementarySurface`: the parameters
  of a point on a curve or surface held as an enum, by the closed form
  of an elementary variant; a B-spline projects the point with
  `Tolerance::DEFAULT`.
- `ParametricCurve2D::eval_d1` and `eval_d2`, mirroring the 3D trait.
- Unit surface normal: `ParametricSurface::normal(u, v, tol)` returns the
  unit `Su × Sv`, recovering the limiting normal at singular points (a
  sphere's poles, a B-spline row collapsing to a smooth point) from second
  derivatives, then from nearby parameters; `None` where no normal can be
  recovered, and at an apex (a cone's, or a collapsed B-spline row whose
  rows meet at an angle), whose meridians each have their own limit. A
  partial counts as zero when the distance it would cover across its
  parameter range is within `tol.confusion`, so rescaling a surface's
  knots does not change its normals. In Python,
  `normal(u, v, tol=None)` on `Cylinder`, `Cone`, `Sphere`, `Torus` and
  `BSplineSurface` (`Plane.normal()` keeps returning its constant
  normal). ADR 0010 records the design.
- Curve local properties (#69): `tangent(u, tol)` returns the unit
  direction of motion (`None` where the speed vanishes, as at a cusp),
  `curvature(u, tol)` the turn rate per unit arc length, and
  `curvature_center(u, tol)` the centre of the osculating circle
  (`None` on a straight run or where the speed vanishes). A derivative
  counts as zero when the distance it would cover across the curve's
  parameter range is within `tol.confusion`, so rescaling a curve's knots
  changes none of the three; the surface normal follows the same rule
  (ADR 0010).
- Closedness queries (#70): `is_closed(tol)` on curves and
  `is_u_closed(tol)` / `is_v_closed(tol)` on surfaces report whether
  travelling returns to its start. Periodicity implies closedness, not
  conversely: B-splines whose ends (boundary curves, for surfaces)
  coincide, clamped or not, and trims spanning a full period count as
  closed without being periodic.
- Continuity order queries (#72): `continuity(u, tol)` on curves and
  `continuity_u/v` on surfaces return the new `Continuity` enum,
  `Order(k)` for `C^k` or `Smooth` where no knot constrains the shape,
  ordered so every `Order(k)` sorts below `Smooth` (`min` picks the
  weaker), with `order()` and `is_at_least(k)` helpers; with `min_continuity` (curves) and `min_continuity_u/v` (surfaces)
  summaries. Elementary shapes are smooth everywhere; B-splines report
  the `degree − multiplicity` their knots guarantee, wrapping periodic
  parameters into the period. A trimmed curve reports its basis order
  inside its range, smooth at its ends, and only the knots inside its
  range in its minimum.
- Surface curvature properties (#73): `curvatures(u, v, tol)` returns
  principal curvatures with directions plus mean and Gaussian
  curvature (`None` where the normal degenerates). Curvatures are
  measured along the unit normal, negative bending away from it. The
  two principal directions are always a perpendicular unit pair,
  including at and near umbilical points.
- One-sided evaluation at B-spline knots and periodic seams (#78):
  `eval_d1_sided`, `eval_derivative_sided` and `eval_d2_sided` on
  B-spline curves and surfaces take a `Side` per direction.
  `Side::Below` takes the span ending at the parameter, `Side::Above`
  the span starting there, as plain evaluation does everywhere except
  at the upper bound. A periodic seam is handled like any other knot,
  whichever bound names it. Plain evaluation paths are untouched.
- Surface reversal, part 1 (#71): `reversed_u()` / `reversed_v()` on
  every elementary surface plus `Surface`, with `u/v_reversed_parameter`
  on every surface type and `Surface` mapping parameters onto the
  reversed shape. Reversal flips one frame axis (frames may now be
  indirect: `Frame3D::from_axes`, `Frame3D::is_direct`,
  `flipped_x/y/z`). Isoparametric curves, projection, conversion to
  B-splines, intersections and curve-on-surface parametrization follow
  the reversed parametrization.
- Python result classes mirroring Rust: `CurveCurveIntersection`,
  `CurveCurveIntersection2D`, `CurveSurfaceIntersection`,
  `SurfaceSurfaceIntersection`, `CurveCurveHit`, `CurveCurveHit2D`,
  `CurveSurfaceHit`, `CurveCurveOverlap`, `IntersectionCurve`,
  `Completeness` (`status`, `reason`, `regions`, `is_complete()`),
  `CompletenessStatus`, `IncompleteReason`, `IntersectionMethod`,
  `CurveCurveRegion`, `CurveSurfaceRegion`, `SurfaceSurfaceRegion`,
  `CurveCurveExtrema` and `CurveCurveExtremum`, with read-only
  properties named as in Rust, reprs and value equality. The enums are
  hashable and have a `name`; an empty result is falsy. Intersection
  results and `CurveCurveExtrema` hold their operands as `first` and
  `second`; `CurveSurfaceIntersection` also names them by role, `curve`
  and `surface`, whichever was the receiver. The stubs type the operands
  as unions of the curve or surface classes.
- Python `intersect` accepts every curve and surface class, B-splines
  and trimmed curves included, in either order, routed as in Rust.
- Python exceptions mirror the Rust error types: `GeomcoreError` (a
  `ValueError`) with `ConstructionError`, `FrameConstructionError`,
  `TrimmedCurveConstructionError`, `KnotInsertionError`, `KnotRemovalError`,
  `DegreeElevationError`, `BSplineSplitError`, `BSplineConversionError`,
  `ParametrizeError`, `SurfaceIsoError` and `DiscretizationCriteriaError`.
  Code that catches `ValueError` keeps working.
- Python binding gaps closed against the Rust API: `Axis2D` and `Frame2D`;
  `frame()` on every curve and surface that has one, `axis()` on `Line3D`,
  `Line2D` and `Cylinder`, `Torus.center()`, `Parabola3D.apex()`;
  `period()` and `is_periodic()` on curves, `u_period()` and `v_period()` on
  surfaces; `eval_d1` and `eval_d2` on curves and surfaces; `bounds()` on
  `Line2D` and `Circle2D`; `Line2D.from_two_points` and `Line2D.from_axis`;
  `extrema` on `Line3D`, `Circle3D`, `Line2D` and `Circle2D`; `point_count`
  on the discretizable curves; `+`, `-` and `*` on points and vectors;
  `Vector3D.normalized`, `square_magnitude`, `is_finite` and
  `angle_with_reference`, the members of `Vector2D` that `Vector3D` has, and
  `Frame3D.axis`, `local_coordinates` and `point_at`.
- Python bindings for the local property, closedness, continuity,
  reversal and one-sided evaluation queries. Every 3D curve class
  (`Line3D`, `Circle3D`, `Ellipse3D`, `Parabola3D`, `Hyperbola3D`,
  `BSplineCurve3D`, `TrimmedCurve3D`) gains `tangent`, `curvature`,
  `curvature_center` (each `Optional`, `None` where the Rust call
  returns `None`), `is_closed`, `continuity` and `min_continuity`. Every
  surface class gains `is_u_closed`, `is_v_closed`, `continuity_u`,
  `continuity_v`, `min_continuity_u`, `min_continuity_v`, `curvatures`
  (a frozen `SurfaceCurvatures` with `min_curvature`, `max_curvature`,
  `min_direction`, `max_direction`, `mean`, `gaussian`, or `None`) and
  `u_reversed_parameter` / `v_reversed_parameter`; the elementary
  surfaces gain `reversed_u` / `reversed_v`. `BSplineCurve3D` and
  `BSplineSurface` gain `eval_d1_sided`, `eval_d2_sided` and
  `eval_derivative_sided`, taking the new `Side` enum (`BELOW`,
  `ABOVE`; hashable, with `.name`). `Frame3D` gains `from_axes`,
  `is_direct` and `flipped_x/y/z`. Continuity queries return the `int`
  order, or `math.inf` where smooth, so `c >= 2` and `min(...)` work
  directly. Every tolerance is the keyword `tol=None`.

### Changed

- The exact, typed intersection moves to `IntersectExact<B>::intersect_exact`
  (breaking): the same per-pair enums, stack-only, never numeric. The
  untyped exact results are renamed `ExactSurfaceIntersection`,
  `ExactCurveSurfaceIntersection` and `ExactCurveCurveIntersection`,
  freeing the plain names for the family results. "Not analytic" is now
  an outcome of `intersect_exact` only.
- Hit fields are renamed `*_param` → `*_parameter` (breaking), on
  `CurveCurveHit`, `CurveSurfaceHit`, `LineTorusHit`, `LineConicHit` and
  `CurveCurveExtremum`.
- Hit parameters are optional (breaking, ADR 0012):
  `CurveCurveHit::first_parameter` / `second_parameter`, the same on
  `CurveCurveHit2D`, and `CurveSurfaceHit::curve_parameter` /
  `surface_parameters` are `Option`s, `Some` only when the solver
  computed the value on its way to the point: every parameter on the
  numeric path, the line's for a line against a quadric, both for two 2D
  lines, and always a trimmed operand's. Recover a `None` with the
  operand's `parameter_of(hit.point)` (curves) or
  `parameters_of(hit.point)` (surfaces). Hits are sorted by parameter
  when every hit has one, and otherwise keep the solver's order. Routed
  circle–circle drops from about 145 to 76 ns per pair, 2D circle–circle
  from 111 to 35 ns, line–cylinder from 184 to 127 ns, against eager
  recovery.
- `ConicConicIntersection::Hits` stores only the parameter on the conic
  the solver substituted (plus the other one where a tangency test
  computed it); the other parameter is `None`.
- `intersect_curve_surface` and `intersect_curve_curve` return the
  family results `CurveSurfaceIntersection` and `CurveCurveIntersection`
  (breaking): read the hits with `.hits()` and the search status with
  `completeness()`; `method()` is always `Numeric`. They stay generic
  over `ParametricCurve3D` / `ParametricSurface`, as the entry point for
  curve and surface types defined outside the crate. For the crate's own
  types, the routed `a.intersect(&b, tol)` solves the pairs with a
  closed form exactly and runs the same search for the others.
- 2D curve projections (`project_point` and `project_points` on
  `Line2D`, `Circle2D` and `TrimmedCurve2D`, and
  `PointExtrema2D::extrema`) return `CurveProjection2D` instead of
  `CurveProjection` (breaking).
- Marching, run by the routed `intersect` of two surfaces without a
  closed form, seeds from the points where the iso-curves of an
  elementary surface meet the other surface: rulings and circles of a
  cylinder or cone, circles of a sphere or torus (not against a torus),
  lines of a plane across the other surface's bounding box. A B-spline
  side adds its own parameter grid, and so do both sides when both
  surfaces are bounded (spheres, tori, B-splines). This replaces the
  parameter grid limited to `[-10, 10]` along unbounded directions, so an
  intersection far from the parameter origin is found. A marched
  result reports `Sampled`, or `Incomplete` where the surfaces lie on
  each other (below).
- Before seeding, marching checks whether a bounded surface lies on the
  other (a 7 × 7 grid of its points and the centres of its cells within
  confusion, normals parallel within the tilt that distance allows
  between samples): surfaces lying on each other, such as a sphere and
  its `to_bspline()` or that copy moved by up to the confusion,
  come back coincident; a patch lying on a larger surface comes back
  `Incomplete` with one `Unseparated` region over the patch. Where the
  surfaces touch along a curve, the curve is traced along the direction
  in which they keep touching and reported with `tangent` set; a touch
  at a point is a contact point. Patches overlapping in part report an
  `Unseparated` region around the overlap. A sphere against its own
  `to_bspline()` takes 0.3 ms instead of 6.8 ms (and 2 spurious curves),
  a B-spline sphere against itself 0.6 ms instead of 270 ms (and 35
  curves); a plane tangent to a cylinder patch 0.1 ms instead of 1.9 ms
  (and no curve).
- A marching trace ends where it leaves the cube of side
  `Tolerance::MODEL_SIZE` centred on the origin, instead of after 2000
  steps: an unbounded curve now runs to the cube, and a long loop closes.
  The step follows the curvature with no fixed cap of 0.1; a step is
  halved when the intersection at its midpoint strays from its chord, or
  when it crosses more than one knot of a B-spline side (so uneven knot
  spans cannot hide a feature from it). Traced curves of low curvature
  take fewer points; each step costs one more Newton solve.
- `parameter_of` on `Circle3D` and `Ellipse3D` no longer normalizes the
  frame's z direction on every call (it is unit by construction): about
  13 ns faster per call, and the parameters may differ in the last bits.
- `transformed` on every elementary surface and on circles, ellipses,
  parabolas and hyperbolas maps each frame axis through the transform,
  so the image's parametrization is the image of the original one
  (`s.transformed(t).eval_point(u, v) == t.apply_point(s.eval_point(u, v))`
  under rotations, translations and mirrors), as `BSplineSurface` and
  `BSplineCurve3D` already did. A mirror now turns the frame indirect
  and the normal to the opposite of the mirrored normal, where it used
  to re-derive a direct frame and keep the side; call `reversed_u` or
  `reversed_v` on the image to keep it.
- `Plane::normal` and `Circle3D::normal` follow orientation: the frame's
  z direction on a direct frame, its negation on an indirect one.
  Unchanged for every plane and circle built by a constructor (all
  direct).
- Python `intersect` returns a result object instead of a
  `(kind, payload)` tuple (breaking); hit parameters are `None` where
  Rust's are, recovered with `operand.parameter_of(hit.point)` or
  `operand.parameters_of(hit.point)`.
- Python `project_point`, `project_points` and `extrema` return
  `CurveProjection`, `CurveProjection2D` and `SurfaceProjection` objects
  with `parameter` (or `u`, `v`), `point` and `distance`, plus the
  `target` they were computed on and the `query_point`, instead of
  tuples (breaking).
- Python `curve_curve_extrema` returns a `CurveCurveExtrema` (`minima`,
  `completeness`, `len()`) instead of a list of tuples (breaking), and
  its arguments are named `first` and `second`, as on the result
  (breaking for callers passing `curve_a=` or `curve_b=`).
- `parametrize_numeric` names its sample count `sample_count` (breaking
  in Python for callers passing `n=`).
- `Vector3D::angle_with_ref` is renamed `angle_with_reference`; the old
  name stays, deprecated.
- Python `intersect` reads its operands by type tests that raise nothing
  on a miss, the surface classes first when the call is on a surface:
  from Python, line–cylinder drops from about 590 to 340 ns per call,
  cylinder–line from 595 to 350 ns, circle–circle from 380 to 320 ns.
  The routed calls that may run the numeric search or marching (a
  B-spline operand, a conic against a torus, two surfaces),
  `curve_curve_extrema` and `parametrize_numeric` release the
  interpreter lock while they compute, so other Python threads run.
- Python: the `Vector3D` unit constructors `x()`, `y()` and `z()` are
  renamed `unit_x()`, `unit_y()` and `unit_z()` (breaking), matching
  `Vector3D::X`, `Y` and `Z`; `x`, `y` and `z` are now the component
  properties, as on `Point3D`.
- Python reprs print floats as Python does (`1.0`, `1e-07`) and show
  coordinates below `1e-15` as `0.0`, so rounding noise such as
  `cos(pi / 2)` no longer clutters them.
- The Python type stubs no longer import `typing_extensions`, and CI
  checks them against the module with `mypy.stubtest`.
- CI runs the Python tests on Python 3.9, the floor of the wheel, as
  well as 3.12.
- CI builds the rustdoc with warnings denied, so a broken intra-doc
  link fails the build.

### Removed

- The 49 per-pair `intersect_*` methods deprecated in 0.29 (breaking):
  use `intersect_exact`, or the routed `intersect`.
- `marching_intersection` and `Surface::intersect_marching` are no longer
  public (breaking): use the routed `intersect`, which marches the pairs
  without a closed form.
- `CurveSurfaceHits` and `CurveCurveHits`, the former results of
  `intersect_curve_surface` and `intersect_curve_curve`, are no longer
  public (breaking). `ConicSurfaceHit` stays public, as the exact results
  carry it.
- Python (breaking): the `(kind, payload)` results, the per-pair
  `intersect_*` methods, and the free functions `intersect_curve_surface`,
  `intersect_curve_curve` and `intersect_marching`, all superseded by
  `intersect`.

## [0.35.0] - 2026-10-19

### Added

- Python bindings for the 0.33 isoparametric curves (#63): `u_iso(u)` and
  `v_iso(v)` on `Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus` and
  `BSplineSurface`, returning `Line3D`, `Circle3D` or `BSplineCurve3D` like
  the Rust methods, with `SurfaceIsoError` as `ValueError`.

## [0.34.0] - 2026-10-13

### Added

- Point coincidence (#65): `Point2D::coincides` and `Point3D::coincides`
  return whether another point is within `tol.confusion`, the test every
  geometric decision already made by hand; in Python `p.coincides(q)`,
  with `Tolerance.default()` when no tolerance is given.

### Changed

- Python value equality (#65): `==` on every Python class compares values
  exactly, field by field, as in Rust, so
  `Point3D(1, 2, 3) == Point3D(1, 2, 3)` is now `True` instead of an
  identity check. Comparing with a value of another class is `False`. As
  `f64` is not `Eq`, the classes are no longer hashable: code that put
  points in a set or used them as dict keys, where they were matched by
  identity, now raises `TypeError`. ADR 0009 records why equality never
  takes a tolerance.

## [0.33.0] - 2026-10-13

### Added

- Isoparametric curves (#63): `u_iso(u)` fixes `u` and returns the 3D
  curve running along `v`, `v_iso(v)` the one running along `u`, each
  parametrized like the surface so the curve at `t` is the surface point
  at `(u, t)` or `(t, v)`. A `Line3D` on a plane and for cylinder and
  cone generators, a `Circle3D` for cylinder and cone sections, sphere
  and torus meridians and parallels, a `BSplineCurve3D` on a
  `BSplineSurface` (the other direction's degree, knots and periodicity,
  rational unless its weights all come out equal), and a `Curve3D` from
  `Surface`. At a cone apex, a sphere pole or a torus pinch the circle
  has radius zero; on a collapsed B-spline edge every pole coincides.
  Errors are the new `SurfaceIsoError`: `NonFinite` for a NaN or
  infinite parameter, `OutOfRange` outside a non-periodic B-spline
  direction (a periodic one wraps the parameter into its period).

## [0.32.2] - 2026-10-13

### Added

- Python bindings for the 0.32 B-spline editing and conversion APIs
  (#58), mirroring the Rust names 1:1 with errors as `ValueError`:
  - `BSplineCurve3D`: `poles`, `weights`, `knots`, `multiplicities`,
    `insert_knot`, `refine_knots`, `elevate_degree`, `remove_knot`,
    `remove_knots`, `to_bezier_pieces`, `split_at`, `extract`,
    `reversed`, `length` and `approximate_within`.
  - `BSplineSurface`: `u_pole_count`, `v_pole_count`, `poles`,
    `weights`, `u_knots`, `u_multiplicities`, `v_knots`,
    `v_multiplicities`, `u_bounds`, `v_bounds`, per-direction
    `insert/refine/elevate/remove`, `to_bezier_patches`, `split_u_at`,
    `split_v_at`, `extract`, `extract_u`, `extract_v`, `reversed_u`,
    `reversed_v`, `approximate` and `approximate_within`.
  - Exact conversions: `to_bspline` on `Circle3D`, `Ellipse3D`,
    `TrimmedCurve3D`, `Sphere` and `Torus`, and `to_bspline_patch` on
    `Plane`, `Cylinder`, `Cone`, `Sphere` and `Torus`.
  - `bounds()` on every 3D elementary curve, `u_bounds()`/`v_bounds()`
    on every surface, the new `BoundingBox3D` (`min`, `max`, `center`,
    `diagonal`, `union`, `enlarged`, `encloses`, `overlaps`,
    `from_points`) and `bounding_box()` on every 3D curve and surface
    type (`None` when unbounded).

## [0.32.1] - 2026-10-13

### Changed

- The project moved to the `geomcore` GitHub organization: the repository is
  now <https://github.com/geomcore/geomcore> (old links redirect) and the
  documentation <https://geomcore.github.io/geomcore/>.

## [0.32.0] - 2026-10-13

### Added

- B-spline editing: knot insertion (`insert_knot`, `refine_knots`),
  degree elevation (`elevate_degree`) and splitting into Bézier pieces
  (`to_bezier_pieces`) on `BSplineCurve3D`, and the same per direction on
  `BSplineSurface` (`insert_u_knot`, `insert_v_knot`, `refine_u_knots`,
  `refine_v_knots`, `elevate_u_degree`, `elevate_v_degree`,
  `to_bezier_patches`). Every edit keeps the shape, works on clamped,
  unclamped and periodic directions, and keeps the poles it does not
  touch bit for bit, rational ones included. Rational poles are edited
  in homogeneous form, up to 128 of them on the stack, so that inserting
  64 knots allocates only the result, and an edit that rewrites every
  row in place, as degree elevation does, turns the homogeneous poles'
  storage into the points'; errors are `KnotInsertionError` and
  `DegreeElevationError`. Inserting knots into a direction of degree 1
  to 3 blends its poles on a path fixed to the degree, whose loops
  unroll: 64 knots go in 8 to 19% faster than through the general path.
- Arc length of B-spline curves: `BSplineCurve3D::length(max_error)`
  integrates the speed over each knot span by adaptive Gauss–Kronrod
  quadrature on the span's Bernstein form, each span getting a share of
  the error in proportion to its width; a sub-range's length is that of
  `extract(first, last)`.
- Knot removal within a tolerance: `remove_knot` and `remove_knots` on
  `BSplineCurve3D`, `remove_u_knot`, `remove_v_knot`, `remove_u_knots` and
  `remove_v_knots` on `BSplineSurface`. Each removal is kept only if a
  proven bound keeps every point within `max_deviation` of the original,
  the bounds of successive removals adding up; the result reports how
  many copies were removed. Rational shapes are bounded about one centre
  for every row of a surface's net, and the seam of a periodic direction
  keeps one copy. Errors are `KnotRemovalError`.
- Least-squares fitting within a deviation, and of surfaces:
  `BSplineCurve3D::approximate_within(points, degree, max_deviation,
  parametrization)` fits with few poles, every point within
  `max_deviation` of the curve at its parameter. Starting from one Bézier
  piece, each round splits the knot spans that deviate the most, so knots
  gather where the points bend; when no span is left to split, the
  fewest poles that fit on knots spread along the parameters are
  searched instead, and when nothing fits, the points are interpolated.
  `BSplineSurface::approximate(points, u_degree, v_degree, u_pole_count,
  v_pole_count, parametrization)` fits a grid of points: the fit
  separates, each row along `v` and then each column of the rows' poles
  along `u`, one factored system per direction serving every line.
  `BSplineSurface::approximate_within` fits a grid within a deviation the
  same way, the rows within 0.8 of it, each direction keeping the fewer
  poles of refined knots and of the fewest poles on spread knots.
- Splitting, extraction and reversal of B-spline curves and surfaces:
  `split_at(u)`, `extract(first, last)` and `reversed()` on
  `BSplineCurve3D`, and `split_u_at`, `split_v_at`,
  `extract(u_range, v_range)`, `extract_u`, `extract_v`, `reversed_u`
  and `reversed_v` on `BSplineSurface`. Pieces are clamped, non-periodic
  B-splines that evaluate as the original over their range, rational
  ones keeping the poles no insertion touched bit for bit; on a periodic
  direction a range may cross the seam and span up to one period. The
  cuts are raised to full multiplicity in one refinement of the spans
  between them. Reversal mirrors the parameter (`first + last - u`) and
  keeps a periodic shape periodic. A non-finite parameter, one outside
  the bounds, an empty range, a periodic range longer than one period or
  a split at an end is a `BSplineSplitError`.
- Elementary curves as exact B-splines (#37): `to_bspline` on `Circle3D`
  and `Ellipse3D` returns a periodic rational quadratic B-spline over
  `[0, 2π]` (four quarter turns, eight poles), and on a `TrimmedCurve3D`
  of a line, circle, ellipse, parabola, hyperbola or `ElementaryCurve` a
  clamped B-spline over the same range, with the same end points. A
  segment becomes degree 1 and a parabola one polynomial quadratic piece,
  both keeping the curve's parameter; a hyperbola becomes one rational
  quadratic piece, and an arc rational quadratic pieces of at most a
  quarter turn, closed exactly when it spans a whole turn. Every point
  lies on the curve; the parameter matches the curve's at the range ends,
  the knots between pieces (at the pieces' boundary angles) and each
  piece's middle, and strays at most 0.0161 rad from a circle's angle in
  between. `ElementaryCurve::to_bspline` and `Curve3D::to_bspline` (a
  B-spline is cloned) return `BSplineConversionError::Unbounded` for a
  line, parabola or hyperbola, and a trimmed curve whose B-spline would
  overflow `f64` returns `NotRepresentable`.
- Elementary surfaces as exact B-splines (#38): `to_bspline_patch(u_range,
  v_range)` on `Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus` and
  `ElementarySurface` returns the patch over those ranges as a B-spline
  surface whose domain is exactly the ranges. A plane patch is bilinear
  and keeps the plane's parameters; the other surfaces are the circle's
  rational quadratic pieces of at most a quarter turn in `u` times their
  profile in `v`, a segment for cylinders and cones (keeping their `v`)
  and circle arcs for spheres and tori. Every point lies on the surface,
  the boundary curves are the surface's at the range ends, and along a
  circular direction the parameter matches the angle at the range ends,
  the knots and each piece's middle. A range spanning a whole turn gives
  a periodic direction; a cone range ending at the apex and a sphere
  range ending at a pole give a collapsed edge, a row of poles on that
  point. `to_bspline` on `Sphere` and `Torus` converts the whole surface,
  and on `ElementarySurface` and `Surface` (a B-spline is cloned) returns
  `BSplineConversionError::Unbounded` for a plane, cylinder or cone, also
  returned for an infinite range. A NaN, empty or out-of-bounds range, or
  one longer than a turn, returns the new `InvalidRange`.
- Bounding boxes of curves and surfaces (#39): the new `BoundingBox3D`
  (`min`, `max`, `center`, `diagonal`, `union`, `enlarged`, `encloses`
  for a point, `overlaps` for another box) and the `Bound3D` trait,
  whose `bounding_box()` returns the axis-aligned box around a shape, or
  `None` for an unbounded one (lines, parabolas, hyperbolas, planes,
  cylinders, cones). Implemented for every 3D curve and surface type,
  trimmed curves, references and the run-time `Curve3D`, `Surface`,
  `ElementaryCurve` and `ElementarySurface`. Elementary boxes are the
  tightest, touching the shape on every side, in closed form and
  without allocating: a circular or elliptic arc from its ends and the
  angles where a coordinate peaks within its range, a trimmed parabola
  or hyperbola from its ends and the parameters where a coordinate is
  stationary, a torus from its radii and the tilt of its axis. A
  B-spline's box is the box of its poles; a trimmed B-spline's, of the
  poles of the knot spans its range covers and the Bézier control points
  of its pieces of those it cuts, so it holds the piece but may be loose.
- Discretization of curves with proven deviations (#17): circles,
  B-splines and trimmed circles and B-splines become points on the curve
  joined by straight segments, each segment reporting a proven upper
  bound on its chord deviation (largest distance from the curve) and
  angle deviation (largest angle to the curve's tangent). The new
  `Discretize3D` trait counts the points for a `DiscretizationCriteria`
  (`point_count`), fills caller-owned slices without allocating
  (`discretize_into`), or returns a `Discretization3D` (`discretize`).
  Circles use the exact deviations; B-splines are split per knot span
  and bounded from each span's Bézier control points. Both end points are
  the curve's points at the range ends; closed curves end on their first
  point, with at least three segments. B-splines bound the angle
  deviation from the control points of their tangent's direction (#19,
  #20; in homogeneous form for rational ones), reporting π where no bound
  was found: with an angle limit, only next to a cusp or (nearly)
  coincident poles; without one, also on pieces turning too much. Ellipses, and trimmed lines, parabolas,
  hyperbolas and run-time `Curve3D`s, use exact closed-form deviations
  of equal parameter steps (#21); a trimmed line is one segment.
  In the plane (#22), `Discretize2D` gives a `Discretization2D` of a
  `Circle2D` or a trimmed 2D curve (lines, circles, run-time `Curve2D`s,
  borrowed curves), its chord deviations being distances in the plane.
  Python: `discretize(max_chord_deviation, max_angle_deviation=None)` on
  `Circle3D`, `Ellipse3D`, `BSplineCurve3D` and `TrimmedCurve3D`,
  returning a `Discretization3D`, and on `Circle2D` and `TrimmedCurve2D`,
  returning a `Discretization2D`.
  Any positive chord deviation is accepted, also below the confusion
  tolerance; the points cannot lie closer to the curve than their
  coordinates' rounding, about 1e-15 times the curve's distance from the
  origin.
- A trimmed curve can borrow its basis (#18): the curve traits
  (`ParametricCurve3D`, `ParametricCurve2D`, `PointExtrema3D`,
  `PointExtrema2D`) are implemented for references to curves, so
  `TrimmedCurve3D<&BSplineCurve3D>` projects onto, intersects and finds
  extrema on a piece of a shared curve without cloning it.
- Points evaluated together with their derivatives (#27): `eval_d1` and
  `eval_d2` on `ParametricCurve3D` return the point and its first one or
  two derivatives, and on `ParametricSurface` the point, `[Su, Sv]` and
  `[Suu, Suv, Svv]`. B-spline curves and surfaces evaluate their basis
  functions once for all of them, 1.8× faster than separate calls on a
  cubic curve and 4× on a bicubic surface; the other types, and new
  curve and surface types, get them from `eval_point` and
  `eval_derivative`. The numeric curve–curve extrema and intersection and
  the marching intersection use them, 12–17% faster on B-splines.
- `Tolerance::MODEL_SIZE` (1e5), the side of the cube centred on the
  origin that models are assumed to fit in, and
  `Tolerance::angular_for_extent(extent)`, the angle below which two
  directions are parallel for geometry reaching `extent` from the pivot
  of the turn: `confusion / extent`, never below `angular`. The default
  angular tolerance, 1e-12, is `confusion / MODEL_SIZE`. Python:
  `Tolerance.MODEL_SIZE` and `Tolerance.angular_for_extent(extent)`.

### Changed

- `InterpParametrization` is renamed `PointParametrization`: it
  parametrizes the points of least-squares approximation as well as of
  interpolation. B-spline interpolation and least-squares approximation
  are rewritten: each data point's row holds only the basis functions
  non-zero there, and all three coordinates solve in one pass over a
  narrow band (half-bandwidth `degree - 1` for interpolation, whose
  averaged knots zero everything further out, and banded Cholesky of
  half-bandwidth `degree` for the normal equations of approximation).
  Interpolating 128 points with a cubic is about 8 times faster, and
  approximating them with 32 poles about 3.5 times faster, with a
  handful of allocations instead of hundreds. Approximation places its
  interior knots along the points' parameters, so that every knot span
  holds some, instead of evenly: unevenly spaced points no longer leave
  poles undetermined.
- B-spline curves and surfaces store their poles once, as points with a
  separate weight per pole when rational, instead of keeping a second,
  homogeneous copy. Evaluation reads them in place, folding each weight
  into its basis function. Points evaluate on a path of their own that
  skips the derivative tables, and the basis recurrence takes its
  reciprocals ahead of the chain that uses them. A cubic curve point is
  about 1.6 times faster and its second derivatives about 1.4 times; on
  surfaces, points gain about 10% and second derivatives about 15%.
  Non-periodic directions of degree 1 to 3 evaluate the basis of a point
  on a path fixed to their degree, which slices the span's knots once and
  unrolls every loop: curve points are about 1.7 times faster again, and
  surface points with both directions so about 1.6 times.
- `discretize` walks the curve once, growing its result, instead of
  counting the points twice and then filling them. A B-spline piece's
  bounds are computed once and handed to its halves, and a rational
  piece builds its tangent's direction from a table of Bernstein product
  factors shared by the whole curve. A cubic B-spline discretizes 2 to 3
  times faster, a rational one about 7 times.
- B-spline angle deviations are tighter, still proven: a step whose
  tangent turns by at most `Θ <= π/2` reports `Θ/2` (the chord lies
  within half the turning of every tangent of the step), and the turning
  rate is also bounded by `|D × D'| / |D|²` from the Bernstein
  coefficients of the tangent's direction `D`, which leaves out its change
  along itself. With an angle limit, cubic B-splines get about 40% fewer
  points, and the representative fixtures less than half as many.
- `BSplineSurface::project_point` holds fewer subpatch nets at once: the
  patches waiting their turn, but the closest, keep only their bounds,
  their nets written again when searched, and a split subpatch hands its
  storage to its first half. A typical query's peak memory drops by 20 to
  45%, at the same speed. A patch of degree 2 or 3 in both directions
  converts to its Bernstein net on a path fixed to the degree, whose loops
  unroll, for 3 to 6% faster projections.

- **Breaking:** `parametrize_on` takes a `Tolerance`, like every other
  query: `curve.parametrize_on(surface, tol)` on all six 3D curve types.
  A curve lies on the surface when the 2D image evaluates within
  `tol.confusion` of it (at least a few ulps of the coordinates far from
  the origin), where a fixed `1e-6` times the curve's distance from the
  origin was used: at the default tolerance a curve 1e-7 to 1e-6 off the
  surface is now `CurveNotOnSurface`. Its angle checks derive from `tol`
  (see the angles below), and a circle counts as a sphere's meridian or parallel when it strays
  at most `tol.confusion` from one. Python: `parametrize_on(surface,
  tol=None)`.
- Angles derive from the confusion tolerance and the extent of the
  geometry involved (ADR 0006): two directions are parallel, coplanar or
  coaxial when turning one onto the other moves that geometry by at most
  `tol.confusion`. A circle, ellipse or torus tilted below that angle now
  behaves as untilted in its intersections with planes, lines, circles,
  conics, cylinders, cones and tori, and in `parametrize_on` a cylinder,
  cone or torus; unbounded geometry (lines, planes, cylinders, cones,
  parabolas, hyperbolas) keeps `tol.angular`. Two circles or conics are
  coplanar when either lies in the other's plane.
- **Breaking:** `Frame2D::new` rejects directions off square by more
  than the angular tolerance (`|x·y| > 1e-12`), where the parametric
  tolerance (1e-9) was used: axes 1e-9 off square move a point across
  the model by 1e-4.
- Numeric curve–surface intersection of B-splines without sampling
  (#28): `intersect_curve_surface` on a B-spline curve, or a trimmed
  curve over one (owned or borrowed), against a B-spline surface
  subdivides the Bézier pieces of both where their control points' boxes
  overlap until each pair holds at most one hit, then solves it by
  Newton's method. It is about 80× faster on a bicubic surface (35 µs
  per cubic curve, from 2.8 ms), and handles rational and periodic
  B-splines. Hits are the sampled search's, except that parameters now
  always lie within the curve's range and the surface's bounds (Newton
  steps could carry a trimmed curve's hit past its ends), periodic ones
  within their period. A tangency, an overlap or a hit at a singular
  point exceeds the search's budget and falls back to the sampled search. The new
  `ParametricCurve3D::as_bspline` and `ParametricSurface::as_bspline`
  (default `None`) expose the B-spline behind a curve or surface.
- Numeric curve–curve intersection and extrema of B-splines without a
  sampling grid (#29): `intersect_curve_curve` and `curve_curve_extrema`
  on two B-spline curves, or trimmed curves over them (owned or
  borrowed), subdivide the Bézier pieces of both curves. Intersection
  pairs pieces whose control points' boxes overlap and halves them until
  each pair meets at most once, then solves it by Gauss–Newton steps:
  about 17× faster on crossing cubics (9 µs per pair, from 157 µs).
  Extrema drop the piece pairs where the distance cannot have a minimum
  and solve those where it is convex: about 3× faster on skew cubics
  (55 µs per pair, from 167 µs). Both handle rational and periodic
  B-splines and return the same hits, multiplicities and minima as
  before, plus crossings and minima the grid could miss between its
  samples. Tangencies, overlaps and other valleys of equal distance
  exceed the search's budget and fall back to the grid.
- **Breaking:** numeric curve–surface intersection reports whether its
  search was complete (#31): `intersect_curve_surface` returns
  `CurveSurfaceHits`, which reads as the slice of its hits (`hits()`,
  `into_hits()`, iteration) and carries a `Completeness`. A B-spline
  curve against a B-spline surface is `Complete` when every pair of
  pieces was proven apart or isolated to one hit. It is `Incomplete`
  when the budget ran out (`IncompleteReason::BudgetExhausted`) or
  pieces never separated at a tangent contact, an overlap or a singular
  point (`Unseparated`), with the `CurveSurfaceRegion`s where hits may
  be missing: curve and surface parameter ranges and a box in space.
  Its hits are then the proven ones together with the sampled search's
  that do not merge into them. Other pairs are `Sampled`. Code that
  only reads the hits keeps compiling unless it names the `Vec` type;
  call `into_hits()` there.
  See ADR 0007. Python is unchanged and still returns the list of hits.
- **Breaking:** numeric curve–curve intersection and extrema report
  whether their search was complete (#32), like the curve–surface case:
  `intersect_curve_curve` returns `CurveCurveHits` (`hits()`,
  `into_hits()`) and `curve_curve_extrema` returns `CurveCurveExtrema`
  (`minima()`, `into_minima()`), both reading as the slice of their
  results and carrying a `Completeness` whose regions are
  `CurveCurveRegion`s: a parameter range on each curve and a box in
  space. Two B-spline curves are `Complete` when every pair of pieces
  was proven apart (for extrema, to hold no minimum) or solved. A tangent
  contact or a straight overlap, where the pieces run along one line
  within `tol.confusion`, no longer exhausts the budget: it is kept as a
  region (`Unseparated`) and the search goes on. A curved overlap or a
  valley of equal distance, such as concentric arcs, still runs out of
  budget (`BudgetExhausted`). Either way the sampled search runs too: its
  hits that do not merge into the proven ones are added, and its seeds
  join the extrema's, one minimum kept per valley. Other pairs, two
  straight curves included, are `Sampled`. Code that only reads the
  results keeps compiling unless it names the `Vec` type; call
  `into_hits()` or `into_minima()` there. See ADR 0007. Python is
  unchanged.
- Marching intersection of two B-spline surfaces seeded from overlapping
  patch pairs (#30): `marching_intersection` on two B-spline surfaces
  (rational or periodic too) halves the pairs of their Bézier patches
  whose control points' boxes overlap until each is dropped, or flat with
  normals provably never parallel (so no closed loop of the intersection
  hides inside it), and converges one seed from each pair that no traced
  curve crosses yet. Isolated loops smaller than the spacing of the old
  seed grid are now found, and two bicubic 6×6 surfaces crossing in six
  curves take 3.4 ms, from 24 ms. Other surfaces keep the grid seeds, and
  so does a B-spline pair whose search runs out of budget (overlapping
  surfaces). A seed counts as lying on a traced curve by its distance to
  the traced points rather than to the fitted curve.
- Mixed elementary and B-spline pairs without sampling (#40):
  `intersect_curve_surface` and `intersect_curve_curve` on an elementary
  curve or surface against a B-spline curve or surface, either way
  round, convert the elementary shape to its exact rational B-spline and
  run the Bézier piece-pair search, so they report `Complete` or
  `Incomplete` instead of `Sampled`. A circle, ellipse, sphere or torus
  is converted whole and a trimmed curve over its range; a line,
  parabola, hyperbola, plane, cylinder or cone over the part inside the
  box of the B-spline's control points grown by a margin, which holds
  every hit, and a cone whose apex lies in it as two patches, one per
  nappe. Parameters are reported on the elementary shape itself. A
  circle against a bicubic surface is about 58× faster (84 µs, from
  4.9 ms), a circle against a cubic 4.4× (15 µs, from 66 µs), a cubic
  against a cylinder 1.4× (62 µs, from 84 µs). `curve_curve_extrema`
  searches a segment against a B-spline the same way, 2× faster (29 µs,
  from 58 µs); a curved elementary curve or a full line against a
  B-spline stays `Sampled`, as the grid is faster there. Two elementary
  shapes stay `Sampled`. The new `ParametricCurve3D::as_elementary` and
  `ParametricSurface::as_elementary` (default `None`) expose the
  elementary shape behind a curve or surface. See ADR 0007.

### Fixed

- Near-coplanar circles and lines are no longer missed. At the default
  tolerance, a circle tilted 1e-10 rad through a plane gave `Empty`
  instead of the whole circle; a line tilted 1e-10 rad through a circle's
  centre, or one whose origin lies far along it, gave `Empty` instead of
  two points; two circles, or a circle and an ellipse, tilted 1e-10 rad
  gave `NotAnalytic` or `Empty` instead of their crossings.
- `Circle3D::parametrize_on` a plane panicked for a circle crossing the
  plane; it returns `CurveNotOnSurface`.
- Two 2D lines nearly parallel and one of them along an axis never
  reported `Parallel` (directions 1e-14 rad apart crossed at 1e11);
  they are parallel within `tol.angular`.
- Two 3D lines up to 1e-6 rad apart were reported `Parallel`, because
  the squared sine of their angle was compared with `tol.angular`; it
  now compares the sine, and lines 1e-9 rad apart cross. The closest
  points between straight curves (`curve_curve_extrema`) agree.
- A line within 1e-6 rad of a cylinder's axis was treated as parallel to
  it for the same reason, losing hits inside the model on a thin
  cylinder; it is parallel within `tol.angular`.

- A line or circle meeting a circle in its plane is found at a loose
  tolerance (#26). With `tol.confusion` above about `1e-6`, a line or
  second circle accepted as lying in the circle's plane was then rejected
  by `parametrize_on`'s fixed check, and the intersection returned
  `Empty` or `NotAnalytic` instead of the crossings.
- Parameters on a circle or ellipse are reported within its `bounds()`,
  `[0, 2π)`, everywhere (#15). Hits against a plane or sphere came back
  anywhere in `(-2π, 2π)`, hits against a cylinder or cone and the
  ellipse's `project_point` and `extrema` in `(-π, π]`, so a caller
  checking them against `bounds()` or an arc's range dropped valid hits.
  The conic–surface hits are now also sorted by parameter, and a hit or
  projection within `tol.confusion` of the seam (along the curve) reports
  `0` whichever side it was found from, in every intersection and
  projection of a circle or ellipse.
- `marching_intersection` keeps the curves that run into a parameter
  bound (#30). Their traces crept along the bound in steps too short to
  fit, and the whole curve was dropped: on two B-spline surfaces, open
  curves crossing the domain went missing. A trace now ends on the bound.

## [0.31.0] - 2026-10-12

### Added

- Exact intersections between any two of a circle, ellipse, parabola and
  hyperbola (#14), returning the new `ConicConicIntersection`: hits with
  the parameter on each curve and a multiplicity (2 or more marks
  tangency), `Coincident` or `Empty`. Coplanar conics substitute one
  curve's parametrization into the other's implicit equation and meet at
  up to four points; ellipse–parabola, ellipse–hyperbola,
  parabola–parabola, parabola–hyperbola and hyperbola–hyperbola were
  `NotAnalytic` before. Conics in crossing planes meet on the line the
  planes share, at up to two points, instead of reporting `NotAnalytic`;
  this covers two circles in crossing planes through `ElementaryCurve`
  and `Curve3D`. A hyperbola is one branch: hits on its mirror branch are
  not reported.
- Python: `intersect` between any two of `Circle3D`, `Ellipse3D`,
  `Parabola3D` and `Hyperbola3D` (except two circles) returns
  `("hits", [(s, Point3D, t, multiplicity)])`, `("coincident", None)` or
  `("empty", None)`, with `s` on the receiver, where unsupported pairs
  raised `TypeError`.

### Changed

- **Breaking:** a circle against an ellipse, parabola or hyperbola, and
  an ellipse against an ellipse, return `ConicConicIntersection`, which
  replaces `CircleConic3DIntersection`, `CircleConicHit`,
  `EllipseEllipse3DIntersection` and `EllipseEllipseHit`. Its hits carry
  a multiplicity, where a tangency was a separate `Tangent` variant, and
  `circle_param`/`conic_param` become `first_param`/`second_param`. In
  Python the same pairs return `"hits"` with the multiplicity instead of
  `"points"`/`"tangent"`, and the receiver's parameter now comes first,
  so `ellipse.intersect(circle)` puts the ellipse parameter first.
- `contains` on an ellipse, parabola or hyperbola decides most points from
  the curve's implicit equation, with rigorous bounds, and only projects
  points within a few tolerances of the curve: a point off the curve is
  rejected 15–250 times faster (2 µs to 8 ns off an ellipse) and one on
  it accepted 1.5–6 times faster. Results are unchanged.

- The 16 exact intersections of a circle, ellipse, parabola or hyperbola
  with a plane, sphere, cylinder or cone now share one implementation.
  Each surface's section by the curve's plane is written once and each
  curve's parametrization substituted into it once, instead of one
  derivation per pair. Results are unchanged up to rounding, with two
  exceptions:
  - a circle lies on a plane when the normals are parallel and its center
    is on the plane, as for the other curves, instead of when five points
    of it are, which is also about 10% faster; its tangency band against
    the plane scales with its diameter, like the ellipse's, instead of its
    radius plus one;
  - a hit of a circle far from the origin is verified with the same
    coordinate-scaled tolerance as the other curves.
  Hyperbola–cone and ellipse–cone are about 5% faster; no pair is slower.
- Tangencies found through a quartic are cheaper to detect: the critical
  points tested for contact now come from the roots already found (one
  bounded Newton search between neighbouring roots and near each complex
  pair) instead of solving the derivative cubic, and the candidates no
  longer allocate. The conic–quadric pairs that solve a quartic are 12–28%
  faster; line–torus and a circle against a coplanar ellipse, parabola or
  hyperbola gain too. Results are unchanged up to rounding; a hit at the
  ellipse parameter ±π may now be reported as −π.

## [0.30.0] - 2026-10-12

### Fixed

- `Circle3D` against a cylinder or cone no longer reports `Circle` (the
  whole circle on the surface) for a circle that only meets it at three
  points (#10). Whether the circle lies on the surface was decided from
  three samples, which fixes a circle on a plane or sphere but not on a
  cylinder or cone, whose section by the circle's plane is a conic. It
  now takes five, like the ellipse.
- Tangent contacts of a circle with an ellipse, parabola, hyperbola,
  cylinder or cone, and of an ellipse with a sphere, cylinder or cone, are
  no longer lost or split when they fall near half a turn from the
  curve's angular origin (#11), and a tangent contact exactly there
  reports multiplicity 2 instead of 1 (#9). These pairs solve a quartic in
  `tan(t/2)`, where such a contact is a root near infinity: rounding
  turned its double root into a complex pair, giving `Empty`, or into two
  hits about 1e-8 apart, and the contact at exactly half a turn was
  always added with multiplicity 1. The quartic is now set up so that no
  root lies near infinity, and tangency is decided by distance, as for the
  parabola and hyperbola. Parameters stay in `(-π, π]`. These
  intersections are about 60% slower (circle–cylinder from ~0.96 to
  ~1.55 µs).
- `intersect_curve_surface`, the numeric path behind `Intersect` for
  B-spline curves and surfaces, no longer misses a crossing at a point
  where the surface's parametrization is singular: a sphere's pole, a
  cone's apex, a B-spline edge collapsed to a point (#12). Newton's
  method stalls there because the surface's partial derivatives vanish or
  lose rank. A seed at a local minimum of the sampled distance that
  Newton cannot converge from is now refined by a one-dimensional
  minimization of the distance from the curve to the surface. Misses are
  about 10% slower.

## [0.29.0] - 2026-10-12

### Added

- Python: `geomcore.curves.TrimmedCurve3D` and `TrimmedCurve2D`, built
  from any curve of the matching dimension and a range, or with
  `segment(start, end)`. They have `basis()`, `first()`, `last()`,
  `bounds()`, the end points, evaluation, `extrema`, `project_point(s)`
  and `contains`. `curve_curve_extrema`, `intersect_curve_curve`,
  `intersect_curve_surface` and `parametrize_numeric` accept a
  `TrimmedCurve3D` and only search its range.
- `Intersect` trait: one entry point for exact intersections,
  `a.intersect(&b, tol)`, implemented in both orders for every pair of
  elementary surfaces; lines against the five elementary surfaces and
  circles and conics against planes, spheres, cylinders and cones; and
  every curve pair with an exact solver: lines with lines, circles and
  conics, circles with circles and conics, ellipses with ellipses, and
  the 2D line and circle pairs. The output is the pair's own result
  type, so dispatch is static and the call costs the same as the
  per-pair method it replaces (`benches/intersect_dispatch.rs`).
- `ElementaryCurve` and `ElementarySurface`: `Copy` enums of the
  closed-form curves and surfaces, reached from `Curve3D` and `Surface`
  with `as_elementary()`.
- `SurfaceIntersection`, the general result of intersecting two surfaces
  of unknown type: `Surface` and `ElementarySurface` implement
  `Intersect` with it. It reports `Empty`, `Coincident`, `NotAnalytic`
  (also for a B-spline), or `Meets` with up to four curves, each flagged
  when tangent, and up to four contact points, held inline in an
  `InlineList`. Every typed surface result converts into it with
  `From`. It never runs a numeric intersection; call
  `intersect_marching` for that.
- `CurveSurfaceIntersection`, the general result of intersecting a curve
  and a surface of unknown type, from `Curve3D` with `Surface` and
  `ElementaryCurve` with `ElementarySurface`, in either order: `Empty`,
  `Coincident`, `NotAnalytic` (a B-spline, or a circle or conic against
  a torus), or `Points` with up to four `ConicSurfaceHit`s ordered by
  curve parameter, a multiplicity of 2 or more marking a tangency.
- `CurveCurveIntersection`, the general result of intersecting two 3D
  curves of unknown type, from `Curve3D` and `ElementaryCurve`: `Empty`
  (also skew or parallel lines), `Coincident`, `NotAnalytic` (a B-spline,
  crossing curve planes, or a conic pair other than two ellipses or a
  circle with a conic), or `Points` with up to four `CurveCurveHit`s
  ordered by first parameter, with both parameters even where the typed
  result omits them. A lone line hit counts as tangent only when the
  directions are parallel: the typed line-circle and line-conic results
  report `Tangent` for any single hit.

- Python: `intersect(other, tol=None)` on `Plane`, `Cylinder`, `Cone`,
  `Sphere`, `Torus`, `Line3D`, `Circle3D`, `Ellipse3D`, `Parabola3D`,
  `Hyperbola3D`, `Line2D` and `Circle2D`. It picks the exact solver by
  the type of `other`, works in either order, and returns the same
  `(kind, payload)` tuple as the matching per-pair method, whose
  parameters keep their roles: in `circle.intersect(line)` each `t` is
  a line parameter. A pair without an exact solver, such as a circle
  and a torus, raises `TypeError`.

### Changed

- **Breaking:** `CurveCurveHit` has a `multiplicity` field, 1 for a
  crossing and 2 for a tangency. `intersect_curve_curve` sets it to 2
  where the tangents are parallel within `tol.angular`. Code that builds
  a `CurveCurveHit` literal must add the field.

### Deprecated

- The per-pair surface methods `Plane::intersect_plane`,
  `intersect_sphere`, `intersect_cylinder`, `intersect_cone`, and the
  same on `Sphere`, `Cylinder`, `Cone` and `Torus`; and the curve
  methods against surfaces, `intersect_plane`, `intersect_sphere`,
  `intersect_cylinder`, `intersect_cone` and `Line3D::intersect_torus`
  on `Line3D`, `Circle3D`, `Ellipse3D`, `Parabola3D` and `Hyperbola3D`;
  and the curve pairs, `intersect_line`, `intersect_circle`,
  `intersect_ellipse`, `intersect_parabola` and `intersect_hyperbola` on
  `Line3D`, `Circle3D` and `Ellipse3D`, and `intersect_line` and
  `intersect_circle` on `Line2D` and `Circle2D`: use
  `Intersect::intersect`.
- Python: the per-pair `intersect_*` methods emit a `DeprecationWarning`
  and are marked deprecated in the stubs; use `intersect`.

## [0.28.0] - 2026-10-12

### Added

- `TrimmedCurve3D` and `TrimmedCurve2D`: a basis curve restricted to a
  parameter range `[first, last]`, such as a line segment, an arc of a
  circle or conic, or a piece of a B-spline. The basis keeps its own
  parametrisation and the range becomes `bounds()`, so
  `intersect_curve_curve`, `curve_curve_extrema` and
  `intersect_curve_surface` only report results inside it. On a periodic
  basis the range may cross the seam, up to one period. A trimmed curve
  is never periodic. The basis type defaults to `Curve3D`/`Curve2D`.
  - `segment(start, end)` builds a line segment parametrised by arc
    length.
  - `project_point` and `contains` are exact: the closest point is either
    one of the basis's distance extrema inside the range or one of the two
    ends.
- `PointExtrema3D` and `PointExtrema2D` traits: the local minima of the
  distance from a point to a curve, implemented by every curve type and by
  the trimmed curves themselves.

### Fixed

- `Plane::intersect_cone` decides tangency on the angle between the
  plane and the nearest generator, within `tol.angular`, for planes
  through the apex too (#8). It used to threshold the cosine of the
  generator angle, which is quadratic in the tilt, so the `TangentLine`
  band was only a third of `tol.angular` wide (4e-13 rad by default), and
  planes just short of tangency came back as `ApexPoint`. The parabola
  band for planes off the apex already used the angle, and both now share
  it.
- `marching_intersection` stops a trace where the intersection crosses
  itself (#7). At such a point the surfaces touch, and a trace used to
  turn onto the other branch and report an arc a second time. At the
  Villarceau tilt of a torus, the section now comes back as the four arcs
  between the two contact points, with no overlap.

## [0.27.0] - 2026-10-11

### Added

- `is_finite()` on `Point2D`, `Point3D`, `Vector2D` and `Vector3D`.
- `Torus::intersect_plane` returns the two Villarceau circles for a
  plane through the center that is tangent to the tube twice (tilted
  `asin(r / R)` from the equatorial plane), as `TwoCircles`. It used to
  report `NotAnalytic`.

### Changed

- **Breaking:** cones are double. `v` was already unbounded and
  evaluation, containment and projection covered both nappes; every
  analytic intersection now does too, instead of keeping only the
  nappe in front of the apex.
  - Line–cone returns the hit on the opposite nappe.
  - `PlaneConeIntersection::Hyperbola` carries both branches, one per
    nappe. A plane beyond the apex yields the circle, ellipse or
    parabola on the other nappe. `PlaneConeIntersection::Empty` is gone.
  - Sphere–cone finds rings on both nappes, so a single ring is now
    reported as two. `SphereConeIntersection::Circle` is gone.
  - `ConeCylinderIntersection::Circle` became `TwoCircles`, one ring
    per nappe.
  - Coaxial cone–cone gains `TwoCircles` for unequal angles.
    `ConeConeIntersection::Empty` is gone. A flipped axis with the
    opposite semi-angle is `Coincident`.
  - `TorusConeIntersection` is now `Rings(Vec<LatitudeRing>)`, with up
    to four rings sorted along the torus axis, each flagged when
    tangent, plus `Empty` and `NotAnalytic`.
  - Python: plane–cone `"hyperbola"` returns `(Hyperbola3D, Hyperbola3D)`;
    cone–cylinder returns `"two_circles"`; cone–cone can return
    `"two_circles"`; sphere–cone no longer returns `"circle"`;
    torus–cone returns `("rings", [(Circle3D, tangent), ...])`.
- `Circle3D::intersect_circle` and `Ellipse3D::intersect_ellipse`
  return `Empty` instead of `NotAnalytic` for curves in distinct
  parallel planes, matching circle–conic.
- **Breaking:** constructors validate their inputs.
  - Every constructor rejects NaN and infinite inputs: radii, focal
    distances, semi-angles, point and direction coordinates, and plane
    coefficients. They used to be accepted and turned every later
    result into NaN. Each construction error enum gains a `NonFinite`
    variant, `FrameConstructionError` included.
  - Conics that a zero size collapses onto a line or segment are
    rejected: parabola focal `<= 0`, ellipse minor radius `<= 0`,
    hyperbola radii `<= 0`. Their evaluation and inversion produced
    NaN or infinity. The variants are renamed
    `ParabolaConstructionError::NonPositiveFocal`,
    `EllipseConstructionError::NonPositiveRadius` and
    `HyperbolaConstructionError::NonPositiveRadius`.
  - Zero radii stay valid where the shape degenerates cleanly: circle,
    sphere, cylinder, torus and cone reference radius.
  - `Plane::from_coefficients` rejects finite coefficients whose plane
    lies beyond `f64` range, with `BadEquation`.
  - Python: the same inputs raise `ValueError`.

### Fixed

- `Plane::new` with a NaN normal and `Sphere::new` or `Circle2D::new`
  with a non-finite center no longer panic; they return `NonFinite`.

- `marching_intersection` (#6):
  - Finds curves that run between its seed grid lines, e.g. a sphere
    against a plane at most heights, or a cylinder against an offset
    sphere. It used to return nothing for them. Grid samples where the
    distance to the other surface is a local minimum are now converged
    onto the intersection by Newton's method.
  - Closed loops are traced once and closed on their seed. They used to
    wind around until the step limit, about 30 times for a radius-2 circle.
  - Open curves, such as a plane–plane line, are no longer fitted as
    loops back to their start.
  - Each step is corrected by Newton's method on both surfaces in the
    plane normal to the tangent, instead of alternating projections, and
    the step length follows the curvature, so small loops and shallow
    crossings trace too.

- Cones with a negative semi-angle no longer produce negative ring
  radii (a panic) in the cone–cylinder, cone–cone, sphere–cone and
  torus–cone intersections.

## [0.26.0] - 2026-10-11

### Added

- `transformed(transform)` on every 3D curve and surface (`Line3D`,
  `Circle3D`, `Ellipse3D`, `Parabola3D`, `Hyperbola3D`,
  `BSplineCurve3D`, `Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus`,
  `BSplineSurface`). Analytic types map their frame and scale their
  radii by the uniform scale factor (a mirror re-derives a right-handed
  frame, reversing the angular direction) and fail only for a
  collapsing transform; B-splines map every control pole. Python:
  same-named methods (analytic types raise `ValueError` on collapse).
- Second-order derivatives for `BSplineSurface::eval_derivative`
  (`(2,0)`, `(1,1)`, `(0,2)`, rational via the quotient rule); the
  fixture goldens now carry them.

### Changed

- B-spline curve and surface evaluation rewritten around the local
  basis functions (The NURBS Book A2.2/A2.3) over the `p+1` (resp.
  `(p+1)(q+1)`) affected poles with fixed-size stack tables, instead of
  copying all poles per call: 10-30x faster on large control nets,
  results unchanged to 1e-14.
- `gauss_newton_1d`/`gauss_newton_2d` evaluate each trial residual once.
- `BSplineSurface::project_point` and `extrema` are exact. Each
  knot-span patch (from a bounding-box hierarchy over the patch poles,
  built on the first projection and visited by increasing box distance)
  is converted to Bernstein form and halved, closest control point
  first. A piece is dropped once its control points are provably
  farther than the best distance so far, and solved by a box-bounded
  Newton minimisation with exact second derivatives once interval
  bounds on the derivatives of the squared distance show it is convex
  or that its minimum lies on one edge or corner; otherwise it is
  halved again. Minima on crease and boundary edges are found on the
  edge curve itself. `extrema` runs the same search on every patch
  without pruning, so it returns every local minimum, its first entry
  being `project_point`'s result. Against dense brute force on 6x6 to
  80x80 nets the previous version missed the global minimum for 7-15%
  of query points (by up to 106 units); on 5,400 random surfaces
  (rational, periodic, creased) the new one is never farther than a
  reference kernel or brute force, at 9-52 us per projection. `extrema`
  searches every patch, so it costs more on large nets (3.4 ms on
  80x80).
- `BSplineCurve3D::project_point` is exact: a branch and bound over the
  knot spans (bounding boxes of the span poles, visited nearest first)
  in which each visited span is scanned exactly. On a span the
  derivative of the squared distance has the sign of a polynomial
  (degree `2p - 1`, `3p - 1` if rational); its rising roots are
  bracketed in Bernstein form (Descartes' rule of signs with de
  Casteljau halving) and refined by a Newton minimisation bounded to
  the bracket, so every minimum is found however close basins sit, and
  the foot point is exact to rounding. On 9,900 random curves (degree
  1-5, rational, periodic, interior knots of any multiplicity) none is
  farther than a dense brute-force search; the previous sampling-based
  search missed 3.9% of them. Speed is unchanged (1-8 us at degree <= 5).
- `BSplineCurve3D::extrema` returns the local minima of the distance
  (plus the ends of an open curve), found by the same exact scan over
  every span; a knot is a minimum when the distance falls into it and
  rises out of it. Its documentation (and `BSplineSurface::extrema`'s,
  in Python) promised all stationary points, which neither returns.
- The bounded Newton minimisation behind B-spline projection accepts a
  step that shrinks the gradient once the squared distance is flat to
  rounding, so the foot converges to full precision (it stalled near
  `sqrt(eps)` in the parameter).
- `Circle3D` intersection with a conic in a parallel, offset plane
  returns `Empty` (it returned `NotAnalytic`).
- `real_roots` (behind every analytic intersection and conic
  projection) scales the polynomial exactly by powers of two, splits it
  into real linear and quadratic factors refined by Newton, and solves
  roots of very different magnitudes each at its own scale. It is about
  2x slower (150-630 ns) for far fewer failures (below).
- Line intersections with spheres, cylinders, cones, tori and in-plane
  conics form their equation about the line's closest approach to the
  geometry, so the tangency band is `tol.confusion` around the chord's
  midpoint wherever the line's origin lies.
- Tangency in curve-surface intersections is decided by distance: a
  line, parabola or hyperbola whose closest approach to a plane, sphere,
  cylinder, cone or torus is within `tol.confusion` of it reports a
  single tangent hit at the touching point. It was decided by how far
  apart the roots fell in parameter space, which construction rounding
  splits into two hits or none (a line tangent to a torus equator was
  missed in 43 of 300 randomized cases).
- `ConicSurfaceIntersection::Coincident` also covers an ellipse,
  parabola or hyperbola lying entirely on a sphere, cylinder or cone
  (a plane section): these reported one to three isolated hits, or
  none.
- `curve_curve_extrema` returns the local minima of the distance over
  the curves' parameter domains, minima at the ends of bounded curves
  included. Two straight curves are solved in closed form, a full line
  against another curve reduces to that curve's distance to the line,
  and other pairs are minimised over the parameter box by a bounded
  Newton with exact second derivatives, each result checked to be a
  true minimum on the box. Unbounded parameters are seeded from the
  curves' own length scale instead of a fixed `[-10, 10]` window.
  Parallel lines, concentric circles and overlapping curves report the
  distance once, with a representative pair (they returned nothing).
  3-16x faster on the pairs measured.
- Python: `eval_derivative` raises `ValueError` for an unsupported
  order (order 0, or above 2 for B-splines and surfaces) instead of a
  `PanicException`.

### Fixed

- B-spline `extrema`/`project_point` silently dropped seeds whose
  Gauss-Newton refinement did not converge, which could return a far
  local minimum (e.g. 44.4 instead of 2.35); neither search seeds any
  more (see Changed).
- Periodic B-spline derivatives just below the last parameter came from
  the first span: wrapping the parameter into the period could round it
  onto the first knot. A parameter already inside the period is left
  unchanged.
- `Circle3D` intersection with an ellipse, parabola or hyperbola whose
  frame's x-axis differs from the circle's lost every hit: the angles
  were solved in the conic's frame but evaluated in the circle's.
- `Ellipse3D`/`Parabola3D`/`Hyperbola3D::contains` agree with
  `project_point`: a point within tolerance of the curve but where
  parameter inversion lands elsewhere (far out on an arm, say)
  was rejected (38, 285 and 54 wrong answers in 1,800 randomized checks
  per curve type).
- `real_roots`: a small nonzero leading coefficient was treated as zero,
  dropping the large root it implies; small roots next to large ones lost
  their precision; huge coefficients produced NaN and panicked; tiny ones
  underflowed into false double roots. Exact triple and quadruple roots
  now report their full multiplicity. On about 7,200 stress polynomials
  per seed, failures against 80-digit arithmetic drop from about 2,000
  (and about 115 panics) to 2, both errors of the reference itself.
- Removed a leftover debug `eprintln!` from the quartic solver that
  printed to stderr on every call since 0.24.0.
- Periodic B-spline curves and surfaces with a seam multiplicity above
  1 evaluated against the wrong poles: pole 0 now starts at the span
  opened by the last copy of the first knot. Pinned by new fixture
  cases verified in exact arithmetic.
- B-spline derivatives exactly on an interior knot now come from the
  span to its right (they came from the left one, which differs at
  C0/C1 knots); on a periodic curve only exactly the last parameter
  closes the last span, other seam copies wrap to the first.
- `Ellipse3D` projection/`extrema` near the parameter `pi`: the
  `tan(t/2)` quartic degenerates there, so on-curve points could
  project up to 10.6 away on a flat ellipse. The reversed quartic now
  covers that half, and every root is Newton-polished on the
  trigonometric stationarity equation.
- `Hyperbola3D` projection/`extrema` roots are Newton-polished on the
  hyperbolic stationarity equation (on-curve points of a thin
  hyperbola projected 3.4e-6 away).
- `Ellipse3D`/`Parabola3D`/`Hyperbola3D::project_point` panicked on a
  query point with a NaN or infinite coordinate; it now returns a NaN
  parameter and distance (and `extrema` an empty list). A hyperbola
  with both radii zero panicked on every projection.
- `curve_curve_extrema` on trimmed curves reported points stuck at a
  domain end as extrema and could miss the true minimum (5.00637
  instead of 5 for concentric arcs); on random trimmed arc pairs its
  minimum was worse than a 400-sample grid in 17 of 60 general and 32
  of 60 concentric pairs. It also missed minima of lines, parabolas and
  hyperbolas beyond parameter 10, and returned nothing for some pairs
  of skew ellipses.
- `Ellipse3D::intersect_ellipse` reported a tangency on an axis vertex
  `sqrt(eps)` off the axis (4.2e-8 for unit-sized ellipses); it is now
  exact.
- `Parabola3D` intersections with a plane missed roots far out on the
  arms (around `u = 3.4e9`), where checking the hit on the surface
  rounds by more than `tol.confusion`.

## [0.25.0] - 2026-09-30

### Added

- Analytic line/circle-vs-conic intersections: `intersect_ellipse/
  parabola/hyperbola` on `Line3D` (plane-piercing graze/miss, coplanar
  quadratic substitution) and on `Circle3D` (coplanar trigonometric
  solve, circular-ellipse coincidence, `NotAnalytic` otherwise),
  sharing `LineConic3DIntersection`/`CircleConic3DIntersection` with
  both curve parameters per hit. Hits on the unmodeled hyperbola
  branch are dropped. Python: per-type methods returning
  `(kind, payload)` tuples.

## [0.24.0] - 2026-09-30

### Added

- Analytic conic-vs-quadric intersections: `intersect_plane/sphere/
  cylinder/cone` on `Ellipse3D`, `Parabola3D`, and `Hyperbola3D`
  (linear trigonometry, quartics in `tan(t/2)`/`t`/`e^t`, tangency
  multiplicities, whole-curve coincidence, nappe filtering), sharing
  `ConicSurfaceIntersection`. Cross-validated against the generic
  solver. Python: per-type methods returning `(kind, payload)` tuples.

## [0.23.0] - 2026-09-30

### Added

- Analytic circle-vs-quadric intersections: `Circle3D::intersect_plane/
  sphere/cylinder/cone` (linear trigonometry for planes/spheres,
  double-angle quartics otherwise, whole-circle detection by
  three-sample containment, tangency multiplicities). Python returns
  `(kind, payload)` tuples.
- Closed-geometry measures: `Circle3D::circumference/disk_area`,
  `Ellipse3D::area`, `Sphere::area/volume` (exact; ellipse perimeter
  deliberately absent — no closed form).

## [0.22.0] - 2026-09-30

### Added

- B-spline least-squares approximation: `BSplineCurve3D::approximate`
  (uniform interior knots, normal equations by band LU with an explicit
  support guard). Python: `BSplineCurve3D.approximate(points, degree,
  num_poles)`.

## [0.21.0] - 2026-09-30

### Added

- Generic 3D curve–curve extrema and intersection
  (`curve_curve_extrema`, `intersect_curve_curve`, over
  `ParametricCurve3D`): coarse grid seeds refined by Gauss-Newton on
  `|C1(s) - C2(t)|^2`, meetings filtered by tolerance. Python:
  `geomcore.curve_curve_extrema` and `geomcore.intersect_curve_curve`.

## [0.20.0] - 2026-09-30

### Added

- Numeric curve-on-surface parametrization fallback
  (`parametrize_numeric`, generic over both traits): uniform samples
  inverted through surface projection, validated against every analytic
  pcurve. Python: `geomcore.parametrize_numeric(curve, surface, n=64)`.

## [0.19.0] - 2026-09-30

### Added

- Generic marching surface-surface intersection
  (`marching_intersection`, also `Surface::intersect_marching`): grid
  seeds projected both ways, adaptive tangent tracing with alternating
  projections, loop/boundary termination, `BSplineCurve3D` fits.
  Approximate by design (fits, no tangency classification, no singular
  points). Python: `geomcore.intersect_marching(a, b)` returning a list
  of `BSplineCurve3D`.

## [0.18.0] - 2026-09-30

### Added

- Points-to-B-spline interpolation: `BSplineCurve3D::interpolate`
  (chordal/centripetal parametrization, averaged knots, banded LU —
  the tracing-to-curve primitive marching SSI will build on). Python:
  `BSplineCurve3D.interpolate(points, degree)`.

## [0.17.0] - 2026-09-30

### Added

- Generic curve–surface intersection (`intersect_curve_surface`, free
  function over `ParametricCurve3D` × `ParametricSurface`): curve samples
  projected onto the surface seed Newton iteration on `C(t) - S(u,v) =
  0`, reporting ordered transversal hits. `project_point` joined the
  `ParametricSurface` trait so every surface works uniformly. Also
  available as top-level `geomcore.intersect_curve_surface` in Python.

## [0.16.0] - 2026-09-30

### Added

- Analytic curve intersections, batch 2: `Line3D::intersect_circle`
  (transversal graze or in-plane 2D reduction) and
  `Circle3D::intersect_circle` (coplanar pairs via shared-plane 2D
  images, `NotAnalytic` otherwise). Python returns `(kind, payload)`
  tuples.

## [0.15.0] - 2026-09-30

### Added

- Analytic curve intersections, batch 1: `Line2D::intersect_line/
  intersect_circle`, `Circle2D::intersect_circle`, and
  `Line3D::intersect_line` (transversal hits with both parameters,
  tangencies, parallel/coincident/skew classification with distances).
  Python returns `(kind, payload)` tuples.

## [0.14.0] - 2026-09-30

### Added

- Point projection, extrema, and containment on `BSplineSurface`:
  dense grid seeds (boundaries included) refined by 2D Gauss-Newton
  (first derivatives only, clamped into non-periodic directions).
  Python returns `(u, v, distance)` tuples. Surface containment is now
  complete across every surface type.

## [0.13.0] - 2026-09-30

### Added

- Point projection, extrema, and containment on `BSplineCurve3D`:
  dense seeds refined by Gauss-Newton (first derivatives only),
  endpoint candidates on open curves, batch variants. Python returns
  `(parameter, distance)` tuples. Curve containment is now complete
  across every curve type.

## [0.12.0] - 2026-09-30

### Added

- Point projection and extrema on conics: `extrema(point, tol)` (all
  stationary points by ascending distance, via quartic/cubic stationarity
  solved with the `math` module) plus `project_point`/`project_points`
  on `Ellipse3D`, `Parabola3D`, and `Hyperbola3D`. Python returns
  `(parameter, distance)` tuples and lists thereof.

## [0.11.0] - 2026-09-30

### Added

- Internal numerical toolbox (`src/math`, zero dependencies): real roots
  of degree ≤ 4 polynomials (stable quadratic, Cardano/trigonometric
  cubic, Ferrari quartic with biquadratic shortcut, Newton polish,
  tolerance clustering with multiplicity so tangencies survive),
  1D Newton, Gauss-Newton drivers for curve/surface projection, Brent
  minimum, and 2x2 solves. Powers projection and extrema next.

## [0.10.0] - 2026-09-29

### Added

- Analytic intersections: `Cone::intersect_cone` for coaxial cones
  (latitude ring, apex point, coincident, empty, `NotAnalytic`).
- Analytic intersections, torus-symmetric batch:
  `Torus::intersect_plane` (meridian/latitude rings),
  `intersect_sphere/cylinder/cone` (coaxial latitude rings via axial
  quadratics with nappe filtering) and `intersect_torus` (coaxial
  rings, coincident). Non-symmetric pairs report `NotAnalytic`.
  Python returns `(kind, payload)` tuples throughout.

## [0.9.0] - 2026-09-29

### Added

- Analytic intersections: `Cone::intersect_cylinder` for coaxial
  cone–cylinder pairs (single latitude ring, apex point for degenerate
  tubes, `NotAnalytic` otherwise). Python returns `(kind, payload)`
  tuples.

## [0.8.0] - 2026-09-29

### Added

- Analytic intersections, Stream 2 batch 3 (symmetric quadrics):
  `Sphere::intersect_cylinder` and `Sphere::intersect_cone` (latitude
  rings, single/tangent/empty) and `Cylinder::intersect_cylinder`
  (generators, tangent, coincident, empty). Each checks its analytic
  criterion first (axis through center / on axis / parallel axes) and
  reports `NotAnalytic` otherwise — the reserved hook for the future
  numeric surface-surface path, never a silent miss. Python returns
  `(kind, payload)` tuples throughout.

## [0.7.0] - 2026-09-29

### Added

- Analytic intersections, Stream 2 batch 2: `Plane::intersect_cylinder`
  (→ circle / ellipse / two lines / tangent line / empty),
  `Plane::intersect_cone` (→ circle / ellipse / parabola / hyperbola
  branch / two apex lines / tangent generator / apex point / empty, with
  single-nappe filtering), and `Line3D::intersect_plane/sphere/cylinder/
  cone` (→ ordered hits with parameters, tangent, single transversal
  hit, empty, coincident). A shared tolerance-aware quadratic solver
  classifies by root separation, so grazing contact is robust across
  scales. Python returns `(kind, payload)` tuples throughout.

## [0.6.0] - 2026-09-29

### Added

- Analytic intersections, Stream 2 batch 1: `Plane::intersect_plane`
  (→ `Line3D`, direction = n1 × n2, with `Parallel`/`Coincident`
  classification), `Plane::intersect_sphere` and
  `Sphere::intersect_sphere` (→ `Circle3D`, with
  `TangentPoint`/`Empty`/`Coincident` classification). Near-degenerate
  configurations classify by `Tolerance` instead of collapsing to noise.
  Python returns `(kind, payload)` tuples, e.g. `("line", Line3D)`,
  `("tangent_point", Point3D)`, `("empty", None)`.

## [0.5.0] - 2026-09-29

### Added

- Closed-form point projection: `project_point(point, tol)` and batch
  `project_points(points, tol)` on `Line3D`, `Circle3D`, `Line2D`,
  `Circle2D` (returning `CurveProjection { parameter, distance }`) and on
  `Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus` (returning
  `SurfaceProjection { u, v, distance }`). Sub-tolerance distances snap
  to `0.0`, so projection agrees with `contains`. Python returns
  `(parameter, distance)` / `(u, v, distance)` tuples.
  (Conic and B-spline projection arrive with the `math` foundations.)

## [0.4.0] - 2026-09-29

### Added

- Point containment: `contains(point, tol)` on all analytic surfaces
  (`Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus`), verified by
  inverse-parameter re-evaluation against `tol.confusion`.
  Python: `contains(point, tol=None)`.
  (`BSplineSurface` containment arrives with numeric projection.)

## [0.3.0] - 2026-09-29

### Added

- Public `Tolerance` API (`confusion`, `angular`, `parametric`, with
  `Tolerance::DEFAULT` matching the kernel's historical internal values);
  internal tolerances now derive from it. Python: `geomcore.Tolerance`.
- Point containment: `contains(point, tol)` on all analytic curves
  (`Line3D`, `Circle3D`, `Ellipse3D`, `Parabola3D`, `Hyperbola3D`,
  `Line2D`, `Circle2D`), verified by inverse-parameter re-evaluation
  against `tol.confusion`. Python: `contains(point, tol=None)`.
  (`BSplineCurve3D` containment arrives with numeric projection.)

## [0.2.0] - 2026-09-29

### Changed

- **Breaking:** renamed value and placement types for `2D`/`3D` consistency
  with the curve names: `Point2` → `Point2D`, `Point3` → `Point3D`,
  `Vector2` → `Vector2D`, `Vector3` → `Vector3D`, `Axis2` → `Axis2D`,
  `Axis3` → `Axis3D`, `Frame2` → `Frame2D`, `Frame3` → `Frame3D`
  (Rust and Python).
- **Breaking (Python):** geometric objects are now constructed with
  `__init__` only (e.g. `Circle3D(Point3D.origin(), Vector3D.z(), 2.0)`);
  the duplicate `.new()` staticmethods were removed.

## [0.1.0] - 2026-07-02

### Added

- Parametric evaluation (point, first and second derivatives, bulk
  `eval_points`) for lines, circles, ellipses, parabolas, hyperbolas, and
  B-spline curves (rational and periodic) in 3D, plus 2D lines and circles.
- Parametric evaluation for planes, cylinders, cones, spheres, tori, and
  B-spline surfaces.
- Rigid transformations (translation, rotation, mirroring, scaling) via a
  single `Transform` type, applicable to all geometry.
- Analytic curve-on-surface parametrization for lines and circles on the five
  elementary surfaces (plane, cylinder, cone, sphere, torus).
- Python bindings (`pip install geomcore`) with native `geomcore.curves` and
  `geomcore.surfaces` submodules, for Python 3.9+ (abi3).
- Golden-fixture validation of all numeric results at 1e-7 tolerance, and
  criterion throughput benchmarks (`docs/benchmarks.md`).

[Unreleased]: https://github.com/geomcore/geomcore/compare/v0.46.0...HEAD
[0.46.0]: https://github.com/geomcore/geomcore/compare/v0.45.0...v0.46.0
[0.45.0]: https://github.com/geomcore/geomcore/compare/v0.44.0...v0.45.0
[0.44.0]: https://github.com/geomcore/geomcore/compare/v0.43.0...v0.44.0
[0.43.0]: https://github.com/geomcore/geomcore/compare/v0.42.2...v0.43.0
[0.42.2]: https://github.com/geomcore/geomcore/compare/v0.42.1...v0.42.2
[0.42.1]: https://github.com/geomcore/geomcore/compare/v0.42.0...v0.42.1
[0.42.0]: https://github.com/geomcore/geomcore/compare/v0.41.0...v0.42.0
[0.41.0]: https://github.com/geomcore/geomcore/compare/v0.40.1...v0.41.0
[0.40.1]: https://github.com/geomcore/geomcore/compare/v0.40.0...v0.40.1
[0.40.0]: https://github.com/geomcore/geomcore/compare/v0.39.1...v0.40.0
[0.39.1]: https://github.com/geomcore/geomcore/compare/v0.39.0...v0.39.1
[0.39.0]: https://github.com/geomcore/geomcore/compare/v0.38.0...v0.39.0
[0.38.0]: https://github.com/geomcore/geomcore/compare/v0.37.0...v0.38.0
[0.37.0]: https://github.com/geomcore/geomcore/compare/v0.36.0...v0.37.0
[0.36.0]: https://github.com/geomcore/geomcore/compare/v0.35.0...v0.36.0
[0.35.0]: https://github.com/geomcore/geomcore/releases/tag/v0.35.0
[0.34.0]: https://github.com/geomcore/geomcore/releases/tag/v0.34.0
[0.33.0]: https://github.com/geomcore/geomcore/releases/tag/v0.33.0
[0.32.2]: https://github.com/geomcore/geomcore/releases/tag/v0.32.2
[0.32.1]: https://github.com/geomcore/geomcore/releases/tag/v0.32.1
[0.32.0]: https://github.com/geomcore/geomcore/releases/tag/v0.32.0
[0.31.0]: https://github.com/geomcore/geomcore/releases/tag/v0.31.0
[0.30.0]: https://github.com/geomcore/geomcore/releases/tag/v0.30.0
[0.29.0]: https://github.com/geomcore/geomcore/releases/tag/v0.29.0
[0.28.0]: https://github.com/geomcore/geomcore/releases/tag/v0.28.0
[0.27.0]: https://github.com/geomcore/geomcore/releases/tag/v0.27.0
[0.26.0]: https://github.com/geomcore/geomcore/releases/tag/v0.26.0
[0.25.0]: https://github.com/geomcore/geomcore/releases/tag/v0.25.0
[0.24.0]: https://github.com/geomcore/geomcore/releases/tag/v0.24.0
[0.23.0]: https://github.com/geomcore/geomcore/releases/tag/v0.23.0
[0.22.0]: https://github.com/geomcore/geomcore/releases/tag/v0.22.0
[0.21.0]: https://github.com/geomcore/geomcore/releases/tag/v0.21.0
[0.20.0]: https://github.com/geomcore/geomcore/releases/tag/v0.20.0
[0.19.0]: https://github.com/geomcore/geomcore/releases/tag/v0.19.0
[0.18.0]: https://github.com/geomcore/geomcore/releases/tag/v0.18.0
[0.17.0]: https://github.com/geomcore/geomcore/releases/tag/v0.17.0
[0.16.0]: https://github.com/geomcore/geomcore/releases/tag/v0.16.0
[0.15.0]: https://github.com/geomcore/geomcore/releases/tag/v0.15.0
[0.14.0]: https://github.com/geomcore/geomcore/releases/tag/v0.14.0
[0.13.0]: https://github.com/geomcore/geomcore/releases/tag/v0.13.0
[0.12.0]: https://github.com/geomcore/geomcore/releases/tag/v0.12.0
[0.11.0]: https://github.com/geomcore/geomcore/releases/tag/v0.11.0
[0.10.0]: https://github.com/geomcore/geomcore/releases/tag/v0.10.0
[0.9.0]: https://github.com/geomcore/geomcore/releases/tag/v0.9.0
[0.8.0]: https://github.com/geomcore/geomcore/releases/tag/v0.8.0
[0.7.0]: https://github.com/geomcore/geomcore/releases/tag/v0.7.0
[0.6.0]: https://github.com/geomcore/geomcore/releases/tag/v0.6.0
[0.5.0]: https://github.com/geomcore/geomcore/releases/tag/v0.5.0
[0.4.0]: https://github.com/geomcore/geomcore/releases/tag/v0.4.0
[0.3.0]: https://github.com/geomcore/geomcore/releases/tag/v0.3.0
[0.2.0]: https://github.com/geomcore/geomcore/releases/tag/v0.2.0
[0.1.0]: https://github.com/geomcore/geomcore/releases/tag/v0.1.0
