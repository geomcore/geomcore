# Parameter ranges are given by trimmed curves

A parameter range on a curve is expressed in one way only: a trimmed curve, whose `bounds()` are the range. Algorithms take no `first, last` arguments; they work over `bounds()`, so every generic numeric algorithm (intersection, extrema, projection, discretization) restricts itself to a trimmed curve's range without knowing trimming exists, and the range is validated once, at construction (seam crossings, ranges longer than a period, ranges outside the basis). The curve traits are implemented for references, so a trimmed curve can borrow its basis: a kernel whose edges share one B-spline builds a `TrimmedCurve<&BSplineCurve>` per edge, a pointer and two numbers, instead of cloning the curve. Operations that need a finite range, such as discretization, are not offered on unbounded curves (lines, parabolas, hyperbolas), so the compiler rejects them until they are trimmed.

## Amendment: one range trait behind the operations

What a trimmed curve needs from its basis is one crate-private trait, `RangeOperations` in `curves/range.rs`: the box, the discretization and the exact B-spline of a range `[first, last]`. Each curve type implements it once for both tiers, its methods calling the math in the operation folders, and the curve enums forward it to the curve they hold. The public operations work over a range: a trimmed curve's `bounds()`, or a bounded curve's own, so the box and the discretization of a whole circle are those of its range over `[0, 2π]`, and a curve type writes its range math once. The B-spline of a whole curve is the exception: a whole circle or ellipse converts to a periodic B-spline, which no range gives, since a trimmed curve is never periodic.

The trait stays crate-private and sealed. Its methods take the `first, last` arguments this decision keeps out of the public interface. And unsealing it would let a curve type from outside the crate have boxes and discretizations whose bounds are its own claims rather than bounds geomcore proved from the curve's definition (ADR 0003), which point classification and the pruning of numeric seeds are to rely on. A user curve that needs a box or a discretization converts to a geomcore curve, a B-spline if nothing else fits.

### Considered options

- **One range trait per operation and tier** (five of them, the state this replaces): rejected. Each curve type implemented every one and each enum forwarded every one, and a missing impl went unnoticed: a trimmed 3D B-spline had no `to_bspline`, and a trimmed `ElementaryCurve` no discretization.
- **Public range methods, or an unsealed trait**: rejected for the reasons above.

## Considered options

- **`first, last` arguments on every range-aware method**: rejected. A symmetric intersection would take two ranges, each method would re-validate them, and it would add a second way to say what a trimmed curve already says.
- **Dropping trimmed curves for range arguments**: rejected for the same reasons, and it would undo the range-restricted numeric algorithms of 0.28 and their Python API.
- **Trimmed curves that own their basis only**: rejected. It forces a kernel to clone a shared basis curve per edge.
