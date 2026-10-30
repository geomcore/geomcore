# The 2D curve tier mirrors the 3D one

Every 3D curve type has a 2D counterpart with the same capabilities: `Line2D`, `Circle2D`, `Ellipse2D`, `Parabola2D`, `Hyperbola2D`, `BSplineCurve2D` and `TrimmedCurve2D`. `Curve2D` has the same variants as `Curve`, and `ParametricCurve2D` offers what `ParametricCurve` does where it makes sense in the plane: tangent, curvature, continuity, closedness, sided evaluation, and conversion to a B-spline. Discretization, bounding boxes, projection and intersection extend to the new types.

Curve-on-surface parametrization is one of the crate's stated purposes, and its result is a 2D curve. Before this decision `Curve2D` held only lines and circles, so `parametrize_on` could not return the image of an ellipse in its own plane, the most common planar edge after lines and circles, and reported `NotAnalytic`, although that image is an ellipse in `(u, v)`. Parabolas and hyperbolas on a plane had the same gap, and a B-spline curve had no 2D curve to become: `parametrize_numeric` returned bare points. A kernel builds edges on faces from these images, so a 2D tier narrower than the 3D one leaves those edges without a curve.

`BSplineCurve2D` reuses the B-spline engine that already serves 3D curves and surfaces: the engine is generic over its poles, so 2D poles are one more pole type, not a second implementation. Above the engine, both curve types wrap one private spline generic over its pole type, which holds the validation, editing, evaluation, length and projection.

## Amendment: the mirror is one generic type per curve

The mirror was first built as two copies of every curve type, a 3D file and a 2D file that kept the same methods by hand. It is now one generic type per curve, over the tier's point: `GenericLine`, `GenericCircle`, `GenericEllipse`, `GenericParabola`, `GenericHyperbola`, `GenericBSplineCurve` and `GenericBezierCurve`, each `<P: TierPoint>`. The public names are type aliases of its two instances, `Circle = GenericCircle<Point>` and `Circle2D = GenericCircle<Point2D>`, so callers name the curves as before. `TierPoint`, implemented by `Point` and `Point2D` only, names the tier's vector, frame and axis. A method both tiers answer is written once on the generic type; constructors and what only one tier has (`Circle::from_axis`, `transformed`, `parametrize_on`) sit in an `impl` of one alias. The capability promise above stays: the two instances answer the same methods because they are the same code.

The trait families stay two: `ParametricCurve` and `ParametricCurve2D`, `PointExtrema` and `PointExtrema2D`, and the operation traits likewise. A single trait generic over the tier would put `<P as TierPoint>::Vector` in every signature a caller reads, and the families really differ: `as_elementary` and unsigned curvature exist in space only, signed curvature in the plane only. Their impls are stamped for both tiers from one body by a crate-private macro, so the code is still written once. `Curve` and `Curve2D` stay two enums, since their variant sets differ. Bounding boxes stay per tier: `BoundingBox` and `BoundingBox2D` are distinct values with their own coordinates. Their formulas are written once all the same, over the coordinates the tier's points already expose to the B-spline engine and a crate-private interface to the tier's box (ADR 0005, amended).

rustdoc shows each curve as its generic struct, with the two aliases as separate pages; the docs and examples name the aliases.

### Considered options

- **Keep two types, generate the second with a macro**: rejected. Each curve's methods would still be documented twice, and no caller could write code once for both tiers.
- **A defaulted type parameter, `Circle<P = Point>`, with `Circle2D = Circle<Point2D>`**: rejected. A default is not used in expressions, so `Circle::new(...)` would be ambiguous between the constructors of the two tiers.
- **One trait generic over the tier**: rejected, for the signatures and the tier-only methods above. It can follow if a consumer needs a bound over both tiers that the generic types do not already give.

## Amendment (2026-10-30): curve–curve intersection and fitting are written once

The curve–curve intersection still had two copies of everything above the curve types: two result types and two hit types that kept the same methods by hand, two sets of readers turning each exact solver's result into them, two descriptions of a routed curve, two copies of the overlap mapping's curve interface and two macro families stamping the `Intersect` impls, while `BSplineCurve` and `BSplineCurve2D` repeated the three fitting constructors. They now follow the rule of this ADR. `GenericCurveCurveIntersection<P>` and `GenericCurveCurveHit<P>` are the result and the hit, and `CurveCurveIntersection`/`CurveCurveIntersection2D`, `CurveCurveHit`/`CurveCurveHit2D` are their aliases, so callers name them as before. `interpolate`, `approximate` and `approximate_within` sit on `GenericBSplineCurve<P>`. Inside the router, one crate-private `CurveBasis<P>`, implemented for every curve that converts into its tier's elementary enum (`ElementaryCurve` or `Curve2D`), carries what the readers and the overlap mapping ask of a curve, so each reader is written once and a coincidence is mapped on the tier's enum. One macro stamps the curve–curve `Intersect` impls of both tiers, and the 2D numeric path maps the hits of the search in space into plane hits and hands them to the same constructor.

Two things stay per tier. The operand traits `CurveOperand` and `CurveOperand2D` differ in what the tier's numeric search needs, `ParametricCurve` in space and an embedding in space in the plane, and `TrimmedCurve`/`TrimmedCurve2D`, outside this change, are two types. One exact constructor serves both tiers, the 3D form: the readers, inlined where the solver's hits are read, build the list, and the constructor takes it whole.

### Considered options

- **One operand trait generic over the tier**: rejected. Its supertrait would depend on the tier, which a trait cannot say without one more trait per tier naming the search's bound; the two traits are a few lines each, and the routing they feed is the shared part.
- **The 2D constructor for both, pushing the hits into the result**: rejected on the instruction counts. It cost 7 to 16 % on routed line–circle and circle–circle pairs in space, while the list form keeps every 2D pair within 2 % of its own former path.

## Considered options

- **Add 2D types only when a parametrization needs them**: rejected. Every 3D curve type appears as the image of some curve on a plane, so this ends at the full mirror anyway, one exception at a time.
- **Return a 3D curve lying in the plane `z = 0`**: rejected. It puts the parameter space in model coordinates, loses the 2D type safety, and every caller would have to drop `z`.
- **Approximate every non-line, non-circle image with a 2D B-spline**: rejected. An exact image exists, and approximating it throws away the exactness that the exact intersection and parametrization paths keep.
