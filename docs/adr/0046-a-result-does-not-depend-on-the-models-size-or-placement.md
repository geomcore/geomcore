# A result does not depend on the model's size or placement

Moving every operand by one similarity (a uniform scale, a rotation, a translation) and scaling `tol.confusion` by the same factor gives the same result, moved: the same hits, curves, overlaps and extrema, the same method and the same completeness. The caller owns the one confusion tolerance of its model (ADR 0006) and scales it with the model; the library adds no length of its own to it.

Three rules follow, each of which a search or a closed form can break on its own:

1. **A length threshold comes from `tol` or from the operands.** No comparison uses a length literal, and no window, step or budget is a fixed number of model units. A dimensionless guard (a normalized coefficient, a sine) may use `tol.angular`.
2. **A rounding floor is measured about the operands, not about the origin.** A quantity built from differences of coordinates is rounded within `EPSILON` times the coordinates and is as large as the differences; a floor of `EPSILON` times the coordinates squared erases it far from the origin. The floor is the rounding of one term times the size of the other.
3. **A window on unbounded geometry is placed on the operands.** A line or a plane matters across the model (ADR 0006), and the model is where the other operand is, not a cube centred on the origin. The window's side stays `tol.confusion / tol.angular`, which scales with the tolerance.

This amends ADR 0006 and the glossary's model size on one point: the cube that bounds the model is no longer centred on the origin. A model may sit anywhere an `f64` coordinate still resolves well below the confusion tolerance; at the default tolerance of 1e-7 that is beyond 1e7 model units, where a coordinate is exact to about 2e-9.

**The judge is `tests/similarity_invariance.rs`.** It takes every sample of `tests/common/samples.rs` (ADR 0027) and every pair of them, moves them by a scale of 2^-20 and 2^20 and by a translation of 1e3 and 1e7, and compares projections, curve–curve extrema and the three families of intersections with the unmoved result. A power of two scales coordinates without rounding, so the scaled comparisons separate a length literal from rounding. On 0.45.0 it failed 166 of 2,044 comparisons, 10 of them for a translation of 1e3, inside the model cube of ADR 0006. The comparisons that do not hold yet are listed in `tests/fixtures/similarity_invariance_known.txt`; a new mismatch fails the test, and so does a listed comparison that now holds, so the list only shrinks.

## Considered options

- **Keep the model cube centred on the origin**: rejected. Results already change for a translation well inside the cube, so the assumption does not protect them, and a model placed in site coordinates sits millions of units from the origin.
- **A tolerance the library derives from the operands' extent**: rejected, as in ADR 0006. Two shapes of different sizes share points that must coincide by one distance, so the tolerance belongs to the model, and only the caller knows the model. A caller that wants a tolerance relative to its model computes it once and passes it.
- **Decisions on interval arithmetic, with a third answer when the tolerance cannot settle them**: not now. `Completeness` already reports a search that did not settle, and the mismatches measured are length literals, floors and windows, which exact arithmetic would not remove. It is reopened when the judge shows a decision that floating rounding alone gets wrong.
- **A script that greps for float literals**: not now. It flags dimensionless constants and misses a window or a budget written as an integer. The judge measures the property itself; a literal that never changes a result is not a defect.
- **Fix every mismatch before landing the judge**: rejected. The list makes each fix a small change that shortens it, and stops a new dependence on size or placement from entering meanwhile.

## Consequences

A fix to a search or a closed form that removes a dependence on size or placement shortens the list in the same change. A new operation is added to the judge with its family. The numeric searches on unbounded geometry and on free-form surface pairs hold most of the list today.
