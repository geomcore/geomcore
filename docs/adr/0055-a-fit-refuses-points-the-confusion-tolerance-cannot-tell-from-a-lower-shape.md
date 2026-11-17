# A fit refuses points the confusion tolerance cannot tell from a lower shape

`Plane::fit`, `Circle::fit`, `Circle2D::fit`, `Sphere::fit` and `Cylinder::fit` (#253) take the points and a `Tolerance`, and return the shape with its residuals, the root mean square and the largest distance from the points to it, in a `Fit`. They refuse, with a `FitError`, fewer points than the shape's degrees of freedom, and points that lie within `tol.confusion` of one point, of their least-squares line, or, for a sphere or a cylinder, of their least-squares plane.

Points on a line determine no plane and no circle; points on a plane determine no sphere (a circle lies on a family of them) and no cylinder (an ellipse lies on two). The question is when points are on a line, and that is a distance decision like every other one in the crate (ADR 0006): against the caller's confusion tolerance, which the caller scales with the model, so the refusal does not depend on the model's size or placement (ADR 0051). Points just outside it are fitted, and a fit that is barely determined shows it in its residuals and its shape, a radius far larger than the points' spread.

The shapes take their orientation from the points' order, the only orientation the points carry: a circle runs through them in order from the first, and the normal of a plane or the axis of a cylinder is the one they turn counterclockwise about. A fit of the same points in the same order is the same shape, moved with them.

## Considered options

- **No tolerance, refusing only a numerically singular system**: rejected. Whether a system is singular in floating point depends on rounding, so points a model cannot tell from a line would be fitted with an arbitrary circle, its centre set by noise, or refused at one placement and not at another.
- **`Tolerance::DEFAULT` inside the fit**, as `Plane::from_three_points` does: rejected. It makes the refusal depend on the model's size, which ADR 0051 rules out for new code.
- **A threshold on the ratio of the covariance's eigenvalues**: rejected. It is dimensionless and needs a constant of its own, and it measures the spread of the points rather than their distance to the line, which is the question.
- **An unoriented shape, its normal or axis either way**: rejected. The sign would follow the eigenvector's, which rounding sets, so the same points moved would fit a shape facing the other way.

## Consequences

A caller with metrology points passes the tolerance of its model, as for every other query. A sphere or a cylinder fit to points on one plane is refused even where the circle they lie on would give a unique cylinder perpendicular to it: the points say nothing about the length along the axis, and a caller wanting that circle asks `Circle::fit`.
