# The 2D curve tier mirrors the 3D one

Every 3D curve type has a 2D counterpart with the same capabilities: `Line2D`, `Circle2D`, `Ellipse2D`, `Parabola2D`, `Hyperbola2D`, `BSplineCurve2D` and `TrimmedCurve2D`. `Curve2D` has the same variants as `Curve`, and `ParametricCurve2D` offers what `ParametricCurve` does where it makes sense in the plane: tangent, curvature, continuity, closedness, sided evaluation, and conversion to a B-spline. Discretization, bounding boxes, projection and intersection extend to the new types.

Curve-on-surface parametrization is one of the crate's stated purposes, and its result is a 2D curve. Before this decision `Curve2D` held only lines and circles, so `parametrize_on` could not return the image of an ellipse in its own plane, the most common planar edge after lines and circles, and reported `NotAnalytic`, although that image is an ellipse in `(u, v)`. Parabolas and hyperbolas on a plane had the same gap, and a B-spline curve had no 2D curve to become: `parametrize_numeric` returned bare points. A kernel builds edges on faces from these images, so a 2D tier narrower than the 3D one leaves those edges without a curve.

`BSplineCurve2D` reuses the B-spline engine that already serves 3D curves and surfaces: the engine is generic over its poles, so 2D poles are one more pole type, not a second implementation. Above the engine, both curve types wrap one private spline generic over its pole type, which holds the validation, editing, evaluation, length and projection, so each public type is a thin layer of one-line delegations.

## Considered options

- **Add 2D types only when a parametrization needs them**: rejected. Every 3D curve type appears as the image of some curve on a plane, so this ends at the full mirror anyway, one exception at a time.
- **Return a 3D curve lying in the plane `z = 0`**: rejected. It puts the parameter space in model coordinates, loses the 2D type safety, and every caller would have to drop `z`.
- **Approximate every non-line, non-circle image with a 2D B-spline**: rejected. An exact image exists, and approximating it throws away the exactness that the exact intersection and parametrization paths keep.
