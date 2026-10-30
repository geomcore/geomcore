# The surface maps a 3D curve to its UV curve, and curves do not depend on surfaces

`ParametricSurface::curve_3d_to_uv(&self, curve: &Curve, tol) -> Result<Curve2D, Curve3DToUvError>` gives the exact UV curve `c` of a 3D curve `C` lying on the surface `S`, in the same parameter: `S∘c = C`, so `self.eval_point(c(t))` is `curve.eval_point(t)` for every `t`. It is the mirror of `uv_curve_to_3d` (ADR 0030), and it replaces the six inherent `parametrize_on` methods of `Line`, `Circle`, `Ellipse`, `Parabola`, `Hyperbola` and `BSplineCurve`. `Curve3DToUvError` replaces `ParametrizeError`, with the same variants: `NotAnalytic` and `CurveNotOnSurface`. The free sampled fallback `parametrize_numeric` becomes `curve_3d_to_uv_numeric`, re-exported at the crate root next to `intersect_curve_surface`.

The method has the same shape as `uv_curve_to_3d`. It is a provided method whose default refuses with `NotAnalytic`, so a surface type defined outside the crate keeps compiling. Each crate surface overrides it with a one-line forward to `surfaces/uv_curve/from_curve.rs`, and `Surface` and `ElementarySurface` forward through their variant lists (ADR 0017). It takes `&Curve`, so the trait stays dyn-compatible.

**The closed forms live on the elementary surfaces.** The table is unchanged: a line or circle on the plane, cylinder, cone, sphere and torus where its image is a 2D line or circle, and every conic or spline on a plane. The forward reads the surface through `as_elementary`, and every curve on an extrusion, revolution, B-spline or Bézier surface is `NotAnalytic`, as it was. `curve_3d_to_uv_numeric` samples those pairs. Mapping a curve onto the extrusion of a line, the mirror of ADR 0030's affine rows, has no consumer yet and is left out.

**Taking `&Curve` adds `Curve::Bezier`.** A Bézier curve in a plane maps its poles exactly as a B-spline does, since both are the same spline underneath. Its image is a `BezierCurve2D` with the same weights, and a pole off the plane gives `CurveNotOnSurface`. On every other surface it is `NotAnalytic`, as a B-spline is.

**Curves do not depend on surfaces.** With `parametrize_on` gone, no file under `src/curves/` names `crate::surfaces`. `tests/module_dependencies.rs` reads the sources and fails if one does, so the rule stays enforced.

## Considered options

- **Keep `parametrize_on` as a deprecated shim**: rejected. The crate is pre-1.0, and the shim would keep the curves importing surfaces until it is removed. The migration is one line, given in the changelog.
- **A free function `curve_3d_to_uv(curve, surface, tol)`**: rejected. The surface owns both directions, as it owns `eval_point` and `parameters_of` at point level.
- **Closed forms on the swept surfaces now**: deferred until a consumer needs them. The numeric fallback answers in the meantime.

## Consequences

Rust callers write `surface.curve_3d_to_uv(&circle.into(), tol)` for `circle.parametrize_on(surface, tol)`, and match on `Curve3DToUvError`. A caller holding a `Curve` no longer dispatches on its variant. Python moves the method from the 3D curve classes to every surface class: `surface.curve_3d_to_uv(curve, tol=None)` accepts any 3D curve class, raises `Curve3DToUvError`, and `geomcore.parametrize_numeric` becomes `geomcore.curve_3d_to_uv_numeric`. The stub-only `_ParametrizeOnMethods` class is gone. The method is now part of `_SurfaceMethods`.
