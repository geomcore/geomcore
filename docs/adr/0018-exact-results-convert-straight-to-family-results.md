# Exact results convert straight to family results

The router reads each pair's typed exact result (`PlaneSphereIntersection`, `LineConicIntersection`, …) straight into its family result (`SurfaceSurfaceIntersection`, `CurveSurfaceIntersection`, `CurveCurveIntersection`, `CurveCurveIntersection2D`), in one step. One crate-private trait, `ExactFamily<B>` in `intersect/routing/exact_family.rs`, does it for every pair the router solves in closed form: the concrete curves and surfaces, and the `ElementaryCurve`, `ElementarySurface` and `Curve2D` enums, which match on their variant first. It receives both operands and their parameter ranges, because a lone line hit on a circle or conic needs the operands to decide its multiplicity, and a coincidence of two curves needs them to map its overlaps. It also receives the router's numeric search as a fallback, and calls it where the configuration has no closed form.

The untyped exact results are removed: `ExactSurfaceIntersection`, `ExactCurveSurfaceIntersection` and `ExactCurveCurveIntersection`, the crate-private `ExactCurveCurve2DIntersection`, and the `IntersectExact` impls on the enums that returned them. `IntersectExact` stays for pairs of concrete types, the nanosecond path of ADR 0011.

ADR 0011 already made `intersect` the answer for shapes whose type is only known at run time: it gives the whole result, exact when a closed form exists and numeric otherwise. The untyped results were what remained of ADR 0001's untyped path. They sat between the typed results and the family results, so every exact pair crossed two conversions and the router held four bridges, one per family, each shaped differently: a blanket impl for surfaces, a macro for curve–surface, a closure table for curve–curve, and hand-written impls for 2D. A change to a typed result had to be followed through the middle tier and its bridge. With one bridge, each typed result has one reading, next to the others, and adding a pair means adding one `into_family` and one line to a macro.

## Considered options

- **Keep the untyped results public**: rejected. They answer the question `intersect` already answers, with less: no numeric fallback, no method, no parameters of a coincidence. Two public answers for run-time shapes would leave callers to pick between them, and the middle tier would stay in every routed call.
- **Keep them, crate-private, as the router's internal step**: rejected. It keeps the two conversions and the four bridges, which is the cost this decision removes; nothing outside the router reads them.
- **Return `None` for no closed form, and let the router fall back**: rejected after measuring. A curve family result holds its hits inline, a few hundred bytes, and moving it out of an `Option` costs a copy on every routed call: up to a third more instructions on the cheapest 2D pairs. With the fallback passed in, the result is built where the caller reads it.

## Consequences

- Breaking for callers of the removed enums and impls. They call `intersect` on the same operands and read the family result (CHANGELOG lists the mapping). A `NotAnalytic` pair now gets the numeric result instead of nothing: a caller that must not pay for a search checks the operand types and calls `intersect_exact` on the concrete pair.
- `exact/conic_2d.rs` builds the typed 2D conic result from the 3D bridge's family result, so the exact module depends on the router's bridge for that one pair family.
- Routed calls do no more work than before, and the enum and surface pairs do less: measured in instructions per call, plane–sphere is a fifth cheaper and line–circle through `ElementaryCurve` a tenth, while every pair stays within a few percent of before. The 2D readers and the curve–curve hit builders are marked always-inline, because the compiler otherwise keeps them out of line and copies the result on return.
- ADR 0001's untyped path is fully superseded.
