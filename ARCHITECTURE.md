# Architecture

A map of the repository: where each concept lives, the rules that decide
where new code goes, and where to start for the common changes. The domain
terms (elementary, exact, numeric, completeness, and so on) are defined in
[`CONTEXT.md`](CONTEXT.md). The decisions behind the rules are recorded in
[`docs/adr/`](docs/adr/).

## The crate at a glance

```text
src/
├── lib.rs               crate docs and the flat re-exports (geomcore::Line, …)
├── prelude.rs           the shape and operation traits, for `use geomcore::prelude::*`
├── point.rs vector.rs   Point/Point2D, Vector/Vector2D
├── frame.rs             Axis, Frame and their 2D forms: where a shape is placed
├── transform.rs         rigid transforms
├── tolerance.rs         Tolerance: every geometric decision takes one
├── projection.rs        results of point projection (CurveProjection, …)
│
├── curves/              the public curve types, one file per type
│   ├── parametric_curve.rs      ParametricCurve, PointExtrema, Side, Continuity
│   ├── parametric_curve_2d.rs   ParametricCurve2D, PointExtrema2D
│   ├── curve.rs curve_2d.rs     the Curve, ElementaryCurve, Curve2D enums
│   ├── line.rs … hyperbola.rs   elementary 3D curves
│   ├── line_2d.rs … hyperbola_2d.rs
│   ├── bspline.rs bspline_2d.rs spline.rs   the two B-spline tiers and their shared type
│   ├── bezier.rs bezier_2d.rs
│   ├── trimmed.rs               TrimmedCurve and TrimmedCurve2D
│   ├── conversion.rs fitting.rs parametrize.rs   to B-spline, from points, onto a surface
├── surfaces/            the public surface types, one file per type
│   ├── parametric_surface.rs    ParametricSurface, SurfaceCurvatures, normal recovery
│   ├── surface.rs               the Surface and ElementarySurface enums
│   ├── plane.rs … torus.rs bspline.rs bezier.rs
│   ├── iso.rs                   helpers shared by the iso curves (u_iso, v_iso)
│   └── conversion.rs fitting.rs
│
├── intersect/           intersections, by method (ADR 0011)
│   ├── routing.rs       Intersect: exact when a closed form exists, numeric otherwise
│   ├── exact/           IntersectExact and the closed-form solvers, by family:
│   │                    curve_curve.rs, curve_surface.rs, surface_surface.rs
│   ├── numeric/         the searches for pairs with no closed form
│   └── results/         the results Intersect returns, Completeness, InlineList
├── bounding_box/        HasBoundingBox: the box type and one file per family
├── discretization/      Discretize: proven point sequences on bounded curves
│
├── elementary/          internal: closed-form formulas of the elementary curves and surfaces
├── bspline/             internal: the B-spline engine (poles, curve evaluation, editing,
│                        tensor-product surfaces, exact conic nets)
└── math/                internal: polynomial roots, Bernstein form, solvers, quadrature
```

Outside `src/`:

| Path | Holds |
|---|---|
| `tests/` | Integration tests through the public API, one file per area (see below) |
| `tests/fixtures/` | Golden JSON data; [its README](tests/fixtures/README.md) says where each file comes from |
| `benches/` | Criterion benchmarks; results are discussed in `docs/benchmarks.md` |
| `examples/` | Small Rust programs (`cargo run --example sample_conics`) |
| `geomcore-py/` | The Python bindings: `src/` (PyO3), `geomcore/` (stubs), `tests/`, `examples/` |
| `docs/` | The user documentation site (`mkdocs.yml`), with `guides/`, `geometry/` and `adr/` |
| `CONTEXT.md` | The glossary |

## Placement rules

1. **One public type per file, named after the type.** `Ellipse` lives in
   `curves/ellipse.rs`, `Ellipse2D` in `curves/ellipse_2d.rs`. Only 2D
   types carry a dimension suffix (ADR 0013), and the 2D files mirror the 3D
   ones (ADR 0014). A type's construction error lives beside it.
2. **Shared behaviour lives on the traits, implemented in the type's
   file.** `ParametricCurve`, `PointExtrema`, `ParametricCurve2D`,
   `PointExtrema2D` and `ParametricSurface` are defined in the
   `parametric_*.rs` files, and each `impl` sits with its type. A method
   that more than one type answers is a trait method, documented once on
   the trait; a type overrides a provided method only for a faster closed
   form, and its inherent `impl` keeps only what is specific to it
   (`Circle::radius`, `BSplineCurve::insert_knot`). Callers bring the
   traits in with `use geomcore::prelude::*` (ADR 0015).
3. **Operation traits are implemented in their operation's folder.**
   `HasBoundingBox` impls live in `bounding_box/`, `Discretize` impls in
   `discretization/`, and `IntersectExact` impls in `intersect/exact/`, all
   grouped by family. To see everything an operation does, open one folder.
4. **Exact intersection solvers are grouped by pair family.**
   `Line`–`Plane` is in `intersect/exact/curve_surface.rs`, next to its
   result type `LinePlaneIntersection`. The general results that `Intersect`
   returns live in `intersect/results/`.
5. **Internal math sits under the shapes.** `elementary/`, `bspline/` and
   `math/` are `pub(crate)`. They work on poles, knots and numbers,
   not on the public shape types; a shape calls into them, never the
   other way round.
6. **Unit tests live in a sibling `tests.rs`.** `curves/line.rs` ends with
   `#[cfg(test)] mod tests;` and its tests are in `curves/line/tests.rs`
   (for a `mod.rs`, the tests are in `tests.rs` beside it). Tests that only
   use the public API go in `tests/` instead.

## Where to start

**A new curve or surface type.**
- Add `curves/<name>.rs` (or `surfaces/<name>.rs`) holding the type, its
  constructor error and its shape-trait impls.
- Re-export it from the folder's `mod.rs` and from `lib.rs`.
- Add it as a variant of the `Curve`/`ElementaryCurve` (or `Surface`) enum.
- Give it a box in `bounding_box/` and a discretization in `discretization/`.
- Bind it in `geomcore-py/src/curves/<name>.rs` (or `surfaces/`), register
  the class in `geomcore-py/src/lib.rs` and add its stub to
  `geomcore-py/geomcore/*.pyi`.

**A new exact intersection pair.**
- Write the solver as `impl IntersectExact<B> for A` in the family file
  under `intersect/exact/`, with its result type, in both argument orders.
- Register it in the dispatch in `intersect/exact/mod.rs`, so that
  `Intersect` routes the pair there instead of to the numeric search.

**A new operation.**
- Give it its own folder at the top of `src/`, like `bounding_box/` and
  `discretization/`: the trait in `mod.rs`, the impls by family.
- Record the term in `CONTEXT.md`.

## The Python bindings

`geomcore-py/src/` uses the core crate's folder names. Each class sits in
the file its core type has, under `curves/` and `surfaces/`. The other
files are:

| File | Holds |
|---|---|
| `lib.rs` | the module wiring (`#[pymodule]`) and small shared argument checks |
| `values.rs` | `Point`, `Vector`, `Axis`, `Frame`, `Tolerance`, `Transform`, `Side` |
| `bounding_box.rs`, `discretization.rs` | the classes of those operations |
| `intersect.rs` | the routed `intersect` every class shares, and `curve_curve_extrema` |
| `results.rs` | the intersection, projection, extrema and curvature results |
| `arguments.rs` | conversions between Python objects and the core `Curve`/`Surface` enums |
| `errors.rs` | the Python exception hierarchy |

## Integration tests

Each file in `tests/` is its own test crate. Its first line says what it
checks.

| Area | Files |
|---|---|
| Evaluation and construction | `combined_evaluation.rs`, `construction_validation.rs`, `precision.rs`, `surface_reversal.rs` |
| Golden fixtures | `vs_fixtures.rs`, `bezier_fixtures.rs`, `discretization_fixtures.rs`, `intersect_2d_fixtures.rs` |
| B-splines and Bézier | `bspline_tiers.rs`, `degree_elevation.rs`, `bezier_surface.rs` |
| Projection | `projection_point.rs` |
| Intersection | `routed_intersections.rs`, `routed_surface_surface.rs`, `bezier_intersections.rs` |
| Bounding boxes | `bounding_box.rs`, `bounding_box_allocation.rs` |
| Discretization | `discretization.rs`, `discretization_allocation.rs` |

## Checks before a pull request

CI runs these; run them locally first:

```sh
cargo fmt --check
cargo clippy --all-targets -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
cargo test
cargo check -p geomcore-py
```

When you change the bindings, also build them with maturin and run
`pytest geomcore-py/tests`.
