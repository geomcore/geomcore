# geomcore

A geometry-only kernel: evaluation, transformations, projection, extrema and intersection of curves and surfaces. Topology (faces, edges, solids) belongs to the caller.

## Shapes

**Elementary curve**:
A line, circle, ellipse, parabola or hyperbola: a curve with a closed-form definition.
_Avoid_: analytic curve, conic (a line is not a conic), primitive

**Elementary surface**:
A plane, cylinder, cone, sphere or torus: a surface with a closed-form definition.
_Avoid_: analytic surface, quadric (a torus is not a quadric), primitive

**Conic**:
A circle, ellipse, parabola or hyperbola: an elementary curve other than the line.

**Quadric**:
A plane, sphere, cylinder or cone: an elementary surface whose implicit equation has degree at most two, the plane being the degree-one case.

**Closed**:
A curve whose two ends meet, or a surface direction whose boundary curves at its two bounds coincide, so travelling along it returns to its start. Periodicity implies closedness, not conversely: a clamped shape with coincident ends is closed without a period.
_Avoid_: periodic for merely closed shapes, seamless

**Lies on**:
A point lies on a curve or surface when it is within `tol.confusion` of it, tested by `contains`.
_Avoid_: inside (reserved for polygons), containment for anything other than this test

**Normal**:
The unit direction of `Su × Sv` at a surface's parameters. A partial counts as zero when the distance it would cover across its parameter range is within `tol.confusion`, so the normal does not depend on how the surface is parametrized. Singular points resolve their limiting normal only if it is unique along the parameter line that collapses there: a sphere's pole reports its axis, a cone's apex reports nothing. A point beside the apex, even within `tol.confusion` of it, still has its meridian's normal.
_Avoid_: raw cross product for the normal, tangent for either partial derivative

**B-spline**:
A curve or surface defined by poles, knots and optional weights. Every curve is either elementary or a B-spline; a surface may also be swept.
A Bézier curve or surface is a B-spline of a single span, and every operation reads it as one.

**Swept surface**:
An extrusion surface or a revolution surface: a surface traced by moving a **profile** curve. An extrusion sweeps it along a direction, `S(u, v) = C(u) + v d`; a revolution turns it about an axis, `S(u, v)` being `C(v)` rotated by `u`. Some profiles sweep an elementary surface (a circle revolved in a plane through the axis sweeps a torus); the swept surface keeps its own parametrization all the same.
_Avoid_: sweep for a profile moved along an arbitrary path (not supported), surface of linear extrusion

**Profile**:
The curve a swept surface moves.
_Avoid_: generatrix, basis curve, directrix

**Trimmed surface**:
A surface restricted to a rectangle of its parameters, a `u` range by a `v` range: the patch a face of a model lies on, without its boundary curves or orientation. It evaluates as its basis, past the ranges too; a periodic direction's range may cross the seam. Its intersection with another surface is the intersection of the bases, cut to the pieces inside the ranges (ADR 0039).
_Avoid_: face (which adds boundary curves and topology), patch for the type name, bounded surface

**UV curve**:
A 2D curve in a surface's `(u, v)` parameter space; its 3D curve is the surface composed with it, `S∘c`, which `uv_curve_to_3d` gives exactly or refuses (ADR 0030; on a spline surface by composing the two, ADR 0032), and `uv_curve_to_3d_within` gives within a deviation as one B-spline in the same parameter (ADR 0035). `curve_3d_to_uv` goes the other way: the UV curve `c` of a 3D curve `C` on the surface, with `S∘c = C` (ADR 0031).
Synonyms: pcurve (the ISO 10303-42 entity name), parameter-space curve.
_Avoid_: `UV` in Rust names (the casing is `Uv`, one word), parametric curve (`ParametricCurve` already means any curve parametrized by `t`), curve on surface, lift, drop

**Same-parameter deviation**:
The largest distance between a 3D curve and the surface over its UV curve at the same parameter, `max_t |Q(t) − S(c(t))|`: what an edge's tolerance bounds, and what `uv_curve_to_3d_within` holds to the asked deviation (ADR 0035). It is stronger than the Hausdorff distance between the two curves, which ignores the parameter.
_Avoid_: same-parameter error, SameParameter

**Offset**:
The surface at a signed distance along another surface's normal. It is an operation, not a stored shape: an elementary surface's offset is elementary, and other surfaces have none yet (ADR 0016).
_Avoid_: shell, thickening

**Tier**:
Space or the plane: where a curve lives. Each curve type is one generic type over its tier's point, `TierPoint` (`Point` in space, `Point2D` in the plane), named `GenericCircle`, `GenericBSplineCurve` and so on, and used through its two aliases: `Circle` in space, `Circle2D` in the plane (ADR 0014). Only the plane's name carries a suffix (ADR 0013).
_Avoid_: dimension for the type parameter, `Circle3D`, 3D tier as a name

**Shape enum**:
`Curve`, `ElementaryCurve`, `Curve2D`, `Surface` or `ElementarySurface`: one variant per shape type, answering every trait method by forwarding to the shape it holds. Its variants are written once, in the enum's variant list (ADR 0017).
_Avoid_: any-curve, dynamic curve, union

**Continuity**:
The order `C^n` up to which derivatives are guaranteed to agree at a parameter, typed `Continuity`: `Order(degree − multiplicity)` at a B-spline knot, `Smooth` everywhere else and on elementary shapes. Every `Order(k)` ranks below `Smooth`, so the minimum over a range is the weakest point; Python reports the order as an `int` and smooth as `math.inf`. It is read from the knots, so the shape may be smoother than reported. The ends of a non-periodic range report smooth, seen from the inside; a periodic seam is a knot like any other. Trimming never lowers it, and a piece that avoids a knot does not inherit it. Derivatives of higher order than the continuity may differ on the two sides of the parameter; one-sided evaluation (see **Side**) reads each.
_Avoid_: smoothness as a yes-or-no answer, regularity

**Reversal**:
Running a curve, or a surface in its `u` or `v` direction, the other way over the same point set; a reversed surface's normal reverses. The **reversed parameter** is where the reversed shape has the point the shape has at `u`: `first + last - u` over bounded parameters (B-spline and Bézier shapes, which reorder their poles and mirror their knots; a circle or ellipse, which flips its frame's y axis, so `2π - u`; a periodic elementary surface direction, `period - u`), and `-u` over unbounded ones (a line negates its direction, a parabola or hyperbola flips its frame's y axis, an unbounded surface direction flips its frame axis), keeping the point at parameter 0. A trimmed curve follows its basis's rule, and a swept surface its profile's. The map is its own inverse; an elementary shape's flipped frame turns indirect (left-handed).
_Avoid_: mirror (a transform of space, not of traversal)

**Transformed shape**:
The image of a curve or surface under a transform, parametrized by the image of its parametrization: `transformed(t)` evaluated at a parameter is `t` applied to the original point there (length parameters scale with a uniform scale). Every frame axis maps through the transform, so a mirror turns a direct frame indirect and the image's normal is the opposite of the mirrored normal, for elementary shapes and B-splines alike. Reversing the image restores the mirrored side. Every shape fails a transform the same two ways: collapsed by a zero scale, or overflowing.

**Bounding box**:
The axis-aligned box around a curve or surface (`BoundingBox` in space, `BoundingBox2D` in the plane), used to cull pairs before intersecting them. An elementary shape's box is the tightest, touching the shape on every side; a B-spline's comes from its control points, holding the shape but possibly loose. An unbounded shape has none. A box encloses a point and overlaps another box.
_Avoid_: AABB, bbox, extent, contains or inside (reserved) for a point in a box

## Precision

**Model size**:
The side of the cube that models are assumed to fit in: 1e5 model units at the default tolerance, `tol.confusion / tol.angular` in general. The cube sits on the model, wherever the model is (ADR 0046). It bounds the model, not its carriers: an infinite line or plane matters only inside it. Within it, `f64` coordinates are exact to far below the confusion tolerance.

**Similarity invariance**:
Moving every operand by one similarity and scaling the confusion tolerance alike gives the same result, moved. A result therefore does not depend on the model's size nor on where it sits (ADR 0046).
_Avoid_: scale-free or unit-independent for this property

**Confusion tolerance**:
The distance below which two points coincide (`tol.confusion`), tested by `coincides`. Every geometric decision is a distance decision against it.
_Avoid_: precision, resolution, epsilon, close or near for points that coincide

**Equal**:
Two values are equal when they are the same value, field by field (`==`). Equality never takes a tolerance: two points that are not equal may still coincide, within the confusion tolerance.
_Avoid_: equal for points that merely coincide, fuzzy or approximate equality

**Parallel**:
Two directions are parallel for a configuration when turning one onto the other moves the geometry by at most the confusion tolerance: an angle below the confusion tolerance over the geometry's extent, or over the model size for unbounded geometry (`tol.angular_for_extent`).

## Intersection

**Exact intersection**:
An intersection computed in closed form, for a pair the closed-form pair table lists: elementary shapes, or a plane with a swept surface. Its result is itself made of elementary curves and points, or of curves read off a swept surface's profile, except a quartic curve, which is fitted within confusion from its closed form (ADR 0036), its pieces that run to infinity exact quartic branches (ADR 0042).
_Avoid_: analytic intersection

**Closed-form pair table**:
The one place that says which pairs of shapes `intersect` solves in closed form, and which solver answers each. Every other pair, and every pair with a B-spline or Bézier shape, goes to the numeric intersection.
_Avoid_: exact routes, analytic pairs

**Quartic curve**:
The quartic where a cylinder or a cone meets another quadric away from the configurations that split it into conics and lines: two cylinders whose axes cross (the bicylindric curve), a sphere and a cylinder or a cone whose axis misses its centre, a cone and a cylinder off its axis (ADR 0040), or two cones in general position (ADR 0041). Also where a ring torus meets a plane (a spiric section) or a sphere in no special position, read along the torus's tube circles instead of generators (ADR 0048). One or two closed loops, read along the generators of the cylinder or cone: each loop runs around it, over an arc of its generators, twice around through the point where the surfaces touch and the loop crosses itself, over an arc that crosses itself there or at a cone's apex on the other surface, or through the apex of a cone that lies on the sphere. Its closed form gives each point exactly; `intersect` fits it within confusion and reports it exact (ADR 0036). Between two cones a loop may run off to infinity along a generator of one parallel to one of the other's and come back from the other end: closed only through its points at infinity, it comes back as quartic branches. So may the one loop where a cylinder whose axis runs along a generator of a cone meets it, each generator of the cylinder meeting the cone once.
_Avoid_: bicylindric curve for anything but two cylinders, Steinmetz curve (the case of equal radii and meeting axes only)

**Spiric section**:
The curve where a plane cuts a torus, a quartic curve. A plane parallel to the axis cuts two ovals while it is nearer the axis than the inner equator, Cassini's ovals at the distance of the minor radius; touching the inner equator, one loop crossing itself there; and one loop between the equators. Through the axis it is two meridian circles, across the axis latitude circles, and in a bitangent plane through the centre the two Villarceau circles.
_Avoid_: toric section (also used for any section of a torus by a surface), Cassini oval for any section but at the distance of the minor radius

**Tube circle**:
The circle of a torus at one angle around its axis, in a half-plane through the axis: a meridian circle of the torus, its `u` iso-curve. A plane or a sphere meets each where a sine wave in the angle around the tube vanishes, which reads their quartic curve. A cylinder or a cone meets each in up to four points, the roots of a quartic, which reads their octic curve.
_Avoid_: poloidal circle

**Octic curve**:
The curve of degree eight where a cylinder or a cone off its axis meets a ring torus (ADR 0049). Read along the torus's tube circles, each holding up to four of its points, between the tube circles where a pair of points appears or vanishes and the points where the surfaces touch. One to four loops: two where a pipe passes through the tube, four where it runs along a diameter, one crossing itself where the surfaces touch or where a cone's apex lies on the torus. A point where they touch and meet nowhere else is a tangent point. Its closed form gives each point exactly; `intersect` fits it within confusion and reports it exact.
_Avoid_: degree-8 curve, octic

**Quartic branch**:
A piece of a quartic curve between two of its points at infinity, where a generator of one cone is parallel to one of the other's, or a generator of a cylinder runs along one of a cone's: exact at every point, unbounded both ways, its parameter running over the real line (ADR 0042). `intersect` returns it as the curve, not a piece cut at a box.
_Avoid_: unbounded component, infinite loop, asymptotic arc

**Common sphere conics**:
The two plane conics where a cone and a cylinder, or two cones, meet when their axes meet at a point as far from both surfaces, so that both touch one sphere centred there (Monge's theorem on quadrics sharing an inscribed sphere). Each conic is a plane section of either surface; the planes share the line through the points where the surfaces touch along the sphere, two of them, one or none.
_Avoid_: reducer curves (one application), degenerate quartic

**Shared generator**:
A line through the common apex of two cones with axes apart that lies on both. There are up to four: on the unit sphere around the apex, each cone's directions are two circles, and each line is where they meet. The cones touch along a line where their circles touch, and cross along the others.
_Avoid_: common ruling

**Viviani's curve**:
The quartic curve of a sphere and a cylinder that touches it from inside, the cylinder's axis and its far side as far from the centre as the sphere's radius: one loop crossing itself where they touch. Classically the sphere of radius `2r` and the cylinder of radius `r` through its centre.
_Avoid_: hippopede (its projection), Viviani's window (the surface it bounds)

**Numeric intersection**:
An intersection traced or refined iteratively, whose curves are fitted approximations. `intersect` runs it only for a pair or configuration without a closed form, and the result says so (its intersection method).
_Avoid_: fallback, approximate intersection

**Intersection method**:
Which solver produced an intersection result: exact (a closed form, always complete) or numeric (a search, with its own completeness). Every result of `intersect` states it.
_Avoid_: route, path, strategy

**Typed result**:
The result `intersect_exact` returns for one pair of concrete types, such as `PlaneSphereIntersection`: an enum of the outcomes that pair can have, on the stack.

**Family result**:
The result `intersect` returns for a family of pairs: `CurveCurveIntersection`, `CurveCurveIntersection2D`, `CurveSurfaceIntersection` or `SurfaceSurfaceIntersection`. The two curve–curve results are the instances of one `GenericCurveCurveIntersection` over the tier (ADR 0014). It holds hits or curves, overlaps or coincidence, the intersection method and the completeness, whatever the operands. The router reads an exact pair's typed result straight into it.
_Avoid_: general result, untyped result

**Not analytic**:
An outcome of `intersect_exact` only: the typed solver for the pair has no closed form for this configuration. `intersect` never reports it: it runs a numeric intersection instead, or a closed form the typed solver lacks (two circles in crossing planes are solved exactly as conics, on the line the planes share).

**Completeness**:
Whether a numeric intersection's hits, or the distance minima of numeric curve–curve extrema, are all of them. Complete: the search proved that no other hit exists. Incomplete: it did not settle everywhere, and says why (its budget ran out, pieces never separated along an overlap or at a tangent contact of a curve with a surface, or a branch of two unbounded surfaces runs beyond the box where its trace was cut) and the regions where hits may be missing. Sampled: the hits were found from samples, so nothing is proven. A tangent contact between two curves, where they stay within confusion of each other and then part, is one hit at its closest approach and stays Complete.
_Avoid_: converged, exact (reserved for exact intersections), fallback

**Coincident**:
The outcome when two shapes overlap entirely, so their intersection is a whole shape rather than curves or points.

**Contact point**:
An isolated point where two surfaces meet, either tangentially or at a singular point such as a cone apex. It has a position but no surface parameters, which a singular point does not define.

**Hit**:
A point where a curve meets another curve or a surface, with its multiplicity (1 for a crossing, 2 or more for a tangency) and its parameters on both operands (each curve's parameter, the surface's `(u, v)`). The point is always present. A parameter is stored only when the solver computed it on its way to the point, and is otherwise absent (`None`): recover it with the operand's `parameter_of(hit.point, tol)` or `parameters_of(hit.point, tol)`.
_Avoid_: intersection point, contact (reserved for surfaces)

**Overlap**:
A parameter range where a curve lies on another curve or on a surface, instead of meeting it at isolated hits. Between two curves it is a range on each curve, each increasing in its own curve's parameter (so a coincidence across a seam is split there), plus whether the two curves run the same way along it (`same_direction`); with opposite directions the start of one range meets the end of the other. On a surface it is a range on the curve.
_Avoid_: coincidence (reserved for two shapes that overlap entirely), segment

## Projection and extrema

**Projection**:
The closest point of a curve or surface to a point: its parameter, the point itself and the distance, as a `CurveProjection` or `SurfaceProjection` from `project_point`. On a curve, `extrema` lists the local minima of the distance from the point, closest first; the projection is the first of them, or a search for that point alone where it is faster (a line in closed form, a B-spline searching its spans closest first). On an elementary surface it is the closed-form inversion `parameters_of` and one evaluation. `contains` agrees with the projection's distance being `0.0`.
_Avoid_: foot (in names), closest point as a type name

**Extrema**:
The local minima of the distance between two curves, closest first, each with the parameters on both curves and the distance: `curve_curve_extrema`. A valley of equal minima (parallel lines, concentric arcs) is represented by one of its points. The result carries its **completeness**, as a numeric intersection does: the curve–curve intersection is read from the minima at zero distance. Distances within the confusion tolerance snap to `0.0`.
_Avoid_: extremum for the level-function extrema that seed marching (those are level extrema), distance for the whole result

## Discretization

**Discretization**:
Points on a curve, in order of parameter, joined by straight segments. Each point keeps its parameter; each segment keeps its chord deviation and angle deviation.
_Avoid_: polyline, polygon (reserved for a closed 2D outline), tessellation (reserved for surfaces), sampling

**Chord deviation**:
The largest distance from a piece of curve to its chord, the straight segment joining the piece's ends. Each segment of a discretized curve carries a proven upper bound on it, never an estimate.
_Avoid_: sag, deflection, chordal error, tolerance (reserved for `Tolerance`)

**Angle deviation**:
The largest angle between a curve's tangent and the chord, over one segment. Proven like the chord deviation; where the curve's speed reaches zero no bound exists, and the segment reports π.
_Avoid_: angular deflection, turning angle

## Polygons

**Polygon**:
A closed 2D outline of straight segments. It encloses a region.

**Point-in-polygon classification**:
Deciding whether a point is inside the region a closed polygon encloses, outside it, or on its boundary: one of Inside, Outside, OnBoundary. A polygon never offers a `contains` method.
_Avoid_: containment, contains (reserved for a point that lies on a curve or surface), point on polygon

## Local properties

**Tangent**:
The unit direction of a curve's motion at a parameter. It is undefined where the speed vanishes, as at a cusp. The speed counts as zero when the distance it would cover across the curve's parameter range is within `tol.confusion`, so the tangent does not depend on how the curve is parametrized.
_Avoid_: derivative (the un-normalized rate), tangent for a surface partial derivative

**Curvature**:
How fast a curve's tangent turns per unit of arc length. Its reciprocal is the radius of the osculating circle, whose centre is the centre of curvature; both are undefined on a straight run and where the speed vanishes. A run is straight where the part of the second derivative that turns away from the tangent is negligible by the same distance rule, or is parallel to the tangent within `tol.angular`.
A curve in space has unsigned curvature `|C' × C''| / |C'|³`: there is no side to turn toward. A curve in the plane has signed curvature `(x'y'' − y'x'') / |C'|³`, positive where the tangent turns left (counterclockwise in a direct plane), negative where it turns right, zero on a straight run. The centre of curvature lies on the left of a left-turning curve and on the right of a right-turning one, so a circle's centre is its own whichever way it runs: a counterclockwise `Circle2D` has curvature `+1/r`, one with an indirect frame `−1/r`.
_Avoid_: bend, curvature for a surface (which has two principal values, not one), unsigned curvature for a 2D curve

**Principal curvature**:
One of the two extremal normal curvatures at a surface point, with its principal direction; the other direction is perpendicular to it. At an umbilical point every direction is principal.
_Avoid_: curvature without qualification on a surface

**Mean curvature**:
Half the sum of the principal curvatures, measured along the unit normal (negative bends away from it).

**Gaussian curvature**:
The product of the principal curvatures. Its sign classifies the point (elliptic, hyperbolic, parabolic); it vanishes where the surface is developable.
_Avoid_: mean curvature for the product

## Evaluation

**Side**:
Which span to evaluate at a knot or periodic seam: `Below` approaches from the span ending there, `Above` from the span starting there (the side plain evaluation takes, except at the upper bound). Both sides agree away from knots, and at a knot they differ only in derivatives of higher order than its **Continuity**.
_Avoid_: left/right (which flip with traversal direction), sidedness
