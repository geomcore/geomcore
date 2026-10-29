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
A curve or surface defined by poles, knots and optional weights. Every curve or surface is either elementary or a B-spline.

**Continuity**:
The order `C^n` up to which derivatives are guaranteed to agree at a parameter, typed `Continuity`: `Order(degree − multiplicity)` at a B-spline knot, `Smooth` everywhere else and on elementary shapes. Every `Order(k)` ranks below `Smooth`, so the minimum over a range is the weakest point; Python reports the order as an `int` and smooth as `math.inf`. It is read from the knots, so the shape may be smoother than reported. The ends of a non-periodic range report smooth, seen from the inside; a periodic seam is a knot like any other. Trimming never lowers it, and a piece that avoids a knot does not inherit it. Derivatives of higher order than the continuity may differ on the two sides of the parameter; one-sided evaluation (see **Side**) reads each.
_Avoid_: smoothness as a yes-or-no answer, regularity

**Reversal**:
Running a surface in the opposite `u` or `v` direction over the same point set, so the normal reverses. An elementary surface flips one frame axis and its frame turns indirect (left-handed); a B-spline reorders its poles and mirrors its knots.
_Avoid_: mirror (a transform of space, not of traversal)

**Transformed shape**:
The image of a curve or surface under a transform, parametrized by the image of its parametrization: `transformed(t)` evaluated at a parameter is `t` applied to the original point there (length parameters scale with a uniform scale). Every frame axis maps through the transform, so a mirror turns a direct frame indirect and the image's normal is the opposite of the mirrored normal, for elementary shapes and B-splines alike. Reversing the image restores the mirrored side.

**Bounding box**:
The axis-aligned box around a curve or surface (`BoundingBox` in space, `BoundingBox2D` in the plane), used to cull pairs before intersecting them. An elementary shape's box is the tightest, touching the shape on every side; a B-spline's comes from its control points, holding the shape but possibly loose. An unbounded shape has none. A box encloses a point and overlaps another box.
_Avoid_: AABB, bbox, extent, contains or inside (reserved) for a point in a box

## Precision

**Model size**:
The side of the cube, centred on the origin, that models are assumed to fit in: 1e5 model units. It bounds the model, not its carriers: an infinite line or plane matters only inside it. Within it, `f64` coordinates are exact to far below the confusion tolerance.

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
An intersection computed in closed form between elementary shapes, whose result is itself made of elementary curves and points.
_Avoid_: analytic intersection

**Numeric intersection**:
An intersection traced or refined iteratively, whose curves are fitted approximations. `intersect` runs it only for a pair or configuration without a closed form, and the result says so (its intersection method).
_Avoid_: fallback, approximate intersection

**Intersection method**:
Which solver produced an intersection result: exact (a closed form, always complete) or numeric (a search, with its own completeness). Every result of `intersect` states it.
_Avoid_: route, path, strategy

**Not analytic**:
An outcome of `intersect_exact` only: the typed solver for the pair has no closed form for this configuration. `intersect` never reports it: it runs a numeric intersection instead, or a closed form the typed solver lacks (two circles in crossing planes are solved exactly as conics, on the line the planes share).

**Completeness**:
Whether a numeric intersection's hits, or the distance minima of numeric curve–curve extrema, are all of them. Complete: the search proved that no other hit exists. Incomplete: it did not settle everywhere, and says why (its budget ran out, or pieces never separated at a tangent contact or an overlap) and the regions where hits may be missing. Sampled: the hits were found from samples, so nothing is proven.
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
