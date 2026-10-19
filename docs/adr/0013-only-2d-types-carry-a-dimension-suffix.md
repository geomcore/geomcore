# Only 2D types carry a dimension suffix

geomcore is a kernel for 3D models: every surface is 3D, and 2D exists for curves in a surface's parameter space and for planar work. So 3D is the unmarked dimension. A type with no suffix is 3D (`Point`, `Vector`, `Axis`, `Frame`, `Line`, `Circle`, `BSplineCurve`, `Curve`, `ParametricCurve`, `TrimmedCurve`, `BoundingBox`, `Plane`, `Surface`, `Transform`), and only its 2D counterpart says `2D` (`Point2D`, `Line2D`, `Curve2D`). Types that combine operands follow the same rule: `LineLineIntersection` and `LineLine2DIntersection`.

The bounding-box trait could not simply drop its suffix. `Bound` would collide with `std::ops::Bound` and with the Python bindings' `pyo3::Bound`, and it says nothing about what the trait provides. It is named `HasBoundingBox` after its one method, `bounding_box`. `Bounded` was rejected because unbounded shapes such as a plane implement the trait too, and return `None`.

This is the rule that the results (`CurveProjection` / `CurveProjection2D`, `CurveCurveHit` / `CurveCurveHit2D`), `ElementaryCurve`, `Transform` and every surface already followed. The 24 public types that carried `3D` were the exception. With two conventions, a reader could not tell whether a name without a suffix was 3D or had simply never been given a 2D version.

The rename ships in one breaking release. The old names remain for one release as deprecated type aliases (`#[deprecated] pub type Point3D = Point;`), the same way `angle_with_ref` was retired, and the Python bindings keep matching aliases for the same period.

## Considered options

- **Suffix on every type that has, or could have, a 2D counterpart** (the previous state, applied consistently): rejected. It means `Plane3D` and `Sphere3D`, or a rule that depends on guessing which types might one day gain a 2D version.
- **Keep the mixed state**: rejected. It is the friction a new reader hits first, and it gets worse as the 2D tier grows (ADR 0014).
- **Collisions with other crates' `Point` / `Vector`**: accepted. A caller who mixes crates renames on import (`use geomcore::Point as GeomPoint`). That is a one-off cost for the caller, while the suffix would be a cost on every name in this crate.
