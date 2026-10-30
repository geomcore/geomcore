# A UV curve maps to an exact curve or is refused

`ParametricSurface::uv_curve_to_3d(&self, uv_curve: &Curve2D, tol) -> Result<Curve, UvCurveTo3DError>` gives the 3D curve that a **UV curve** `c`, a curve in the surface's `(u, v)` parameter space, traces on the surface `S`, in the same parameter: `S∘c`, whose point at `t` is `self.eval_point(c(t))` for every `t`. It is the inverse of `parametrize_on`. It answers with a curve of `Curve` only when that curve is exact; every other pair is refused, never approximated, so a caller that needs an approximation chooses it.

The surface owns the mapping, as it owns `eval_point` and `parameters_of` at point level. It is a provided method of `ParametricSurface` whose default refuses with `NotAnalytic`, so a surface type defined outside the crate keeps compiling. Each crate surface overrides it with a one-line forward to the pair table in `surfaces/uv_curve/to_curve.rs`, and `Surface` and `ElementarySurface` forward through their variant lists (ADR 0017). The method takes `&Curve2D`, not `impl Into<Curve2D>`, so the trait stays dyn-compatible: marching works on `&dyn ParametricSurface`. Curves do not depend on surfaces for it: nothing about it lives under `curves/`.

**Naming.** Both ends are explicit: "uv curve" is the 2D side, named by the parameters the crate already uses everywhere (`eval_point(u, v)`, `u_iso`, `v_iso`), so it is neither jargon nor an abbreviation; "3d" is the 3D side. The mirror swaps the halves: `uv_curve_to_3d` returns `S∘c`, and the follow-up `curve_3d_to_uv` returns the `c` with `S∘c = C`, a valid identifier where `3d_curve_to_uv` is not. The Rust casing treats "Uv" as one word and writes "3D" as in `Curve2D` and `Vector2D`: `UvCurveTo3DError`. The B-rep literature calls a UV curve a pcurve (the ISO 10303-42 entity), and the glossary keeps that name as a synonym. Rejected: `uv_to_curve`, since "uv → curve" can read as "uv points to a curve" and relies on knowing that a bare "curve" means 3D; `pcurve_to_curve`, which misreads at a glance; `parametric_to_3d`, which has no valid mirror (`3d_to_uv` starts with a digit) and collides with the `ParametricCurve` trait; `parameter_space_to_curve`, too long; and lift / drop, which are metaphors.

## The pair table

| Curve | Surface | 3D curve |
|---|---|---|
| any | plane | the same kind, placed by the plane's frame |
| any | extrusion of a line | its affine image, where it keeps its kind in the same parameter |
| line along `v`, `u = u0` | any other | `u_iso(u0)`, run from `v0` at rate `±1` |
| line along `u`, `v = v0` | any other | `v_iso(v0)`, run from `u0` at rate `±1` |
| straight spline along `u` or `v` | any other | the iso curve, run at the spline's signed speed `σ`: a spline at any speed, a line at speed 1 or else the degree-1 B-spline over the UV curve's bounds, a conic at speed 1 only |

A line is along `u` or `v` when its direction is within `tol.angular` of the axis. Its fixed parameter is then taken from its origin. `u_iso` and `v_iso` are the special case `c = 0`, rate `1`, and `uv_curve_to_3d` is built on them, not the other way round.

**A straight spline** is a degree-1 B-spline or Bézier curve with two poles `p0`, `p1` over `[a, b]` and equal or no weights: `p0 + (t - a) / (b - a) (p1 - p0)`. Files exchanged between kernels often store a straight UV curve that way, so it is read as a line along `u` or `v` instead of being refused. Its speed `|p1 - p0| / (b - a)` need not be `1`, so the iso curve runs at `t ↦ c + σ t` with `σ` any nonzero real. A B-spline or Bézier iso curve maps its knots affinely, exact at any `σ`. A `Line` runs at speed `1`, so at any other `σ` the line iso curve is given as the degree-1 B-spline through its points at the UV curve's bounds, exact over them. A circle, ellipse, parabola or hyperbola has no other speed, so it needs `|σ| = 1` within `tol.angular`, or the pair is refused. Unequal weights run the segment at a varying speed, which is refused. Coincident poles, within `tol.parametric`, are refused as `Collapsed`.

**Reparametrizing an iso curve** by `t ↦ c + s t`, `s = ±1`, keeps every kind it can. A line moves its origin and flips its direction, a circle turns or mirrors its frame, and a B-spline shifts and mirrors its knots. A Bézier curve keeps arbitrary bounds, so the shifted one is a Bézier curve too. Every elementary curve is rebuilt from its point, `s` times its first derivative and its second derivative at `c`. It is a curve of its kind in `t` only when those vectors fit the kind's own parametrization. An ellipse fits from either end of its major axis (`c` a multiple of `π`) and a parabola or hyperbola from its vertex only (`c = 0`). Elsewhere the curve in that parameter is not a conic of `Curve`, so the pair is refused.

**A degenerate iso curve**, at a sphere's pole, a cone's apex, where a torus or revolution surface meets its axis, or on a collapsed row of a spline surface, is refused with `Collapsed`. The iso curves themselves accept it as a zero-radius circle, and `SurfaceIsoError` has no degenerate case, so the check is the mapping's own.

**The extrusion of a line** is the affine map `(u, v) ↦ origin + u e + v d`, oblique when `e · d = k ≠ 0`. Its principal axes are `(1, 1)/√2`, stretched by `√(1 + k)`, and `(1, -1)/√2`, shrunk by `√(1 - k)`. The property test confirms each of these outcomes:

- A B-spline or Bézier curve maps its poles and is always kept.
- A line keeps unit speed only along `u` or `v`.
- A circle becomes an ellipse only when its `x` axis lies along the stretched principal axis, either way. Its cosine vector must stay the longer, as the ellipse's major axis.
- An ellipse is kept along either principal axis while its cosine image stays the longer one: radii 2 and 1 along the shrunk axis give an ellipse, 2 and 1.9 are refused.
- A hyperbola has no order between its radii and is kept along either principal axis.
- A parabola needs a unit tangent at its vertex perpendicular to its axis, which no shear keeps.

On a plane, or an extrusion perpendicular to its line, every kind is kept.

## Considered options

- **A method on every 2D curve type, as issue #152 asked**: rejected. The curves would depend on surfaces. The pair table reads iso curves and affine maps, which are surface knowledge. On the surface, the pair of point-level mappings (`eval_point`, `parameters_of`) gets its curve-level pair in one place.
- **A composed curve type that evaluates any curve through any surface**: deferred. It answers every pair, but with a curve that is neither elementary nor a B-spline, which no operation reads yet. It would be added when a consumer needs it, and the exact table stays its fast path.
- **Approximating the refused pairs**: rejected. A tolerance-dependent fit under a name that promises exactness would hide the choice from the caller.
- **`impl Into<Curve2D>`**: rejected, as it would make the trait not dyn-compatible.
- **`TrimmedCurve2D`**: not accepted, since it is not a `Curve2D` variant. A caller maps its basis and trims the result in the same parameter.

## Consequences

A surface type defined outside the crate gets `NotAnalytic` until it overrides the method. A new crate surface adds its iso curves to the table (`IsoCurves` in `to_curve.rs`) and its one-line forward. Python exposes the method on every surface class, through the `surface` method set, and raises `UvCurveTo3DError`. A follow-up replaces `parametrize_on` with `ParametricSurface::curve_3d_to_uv(&self, curve: &Curve, tol) -> Result<Curve2D, Curve3DToUvError>` in `surfaces/uv_curve/from_curve.rs`, so the surface owns both directions.
