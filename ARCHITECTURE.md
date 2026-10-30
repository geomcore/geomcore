# Architecture

A map of the repository: where each concept lives, the rules that decide
where new code goes, and where to start for the common changes. The domain
terms (elementary, exact, numeric, completeness, and so on) are defined in
[`CONTEXT.md`](CONTEXT.md). The decisions behind the rules are recorded in
[`docs/adr/`](docs/adr/).

## The crate at a glance

```text
src/
├── lib.rs               crate docs and the flat re-exports: each item's one public path
├── prelude.rs           the shape and operation traits, for `use geomcore::prelude::*`
├── point.rs vector.rs   Point/Point2D, Vector/Vector2D
├── frame.rs             Axis, Frame and their 2D forms: where a shape is placed
├── transform.rs         Transform, and the Transformable trait of the shapes it moves
├── tolerance.rs         Tolerance: every geometric decision takes one
├── projection.rs        results of point projection (CurveProjection, …)
├── tier.rs              TierPoint: the two tiers (space, plane) the curve types are generic over
├── shape_enum.rs        the macros that expand a shape enum's variant list (ADR 0017)
│
├── curves/              the public curve types, one file per type
│   ├── parametric_curve.rs      ParametricCurve, PointExtrema, ReversibleCurve, Side, Continuity
│   ├── parametric_curve_2d.rs   ParametricCurve2D, PointExtrema2D
│   ├── curve.rs curve_2d.rs     the Curve, ElementaryCurve, Curve2D enums
│   ├── line.rs … hyperbola.rs   elementary curves, each one generic type for both tiers
│   ├── bspline.rs spline.rs     B-spline curves of both tiers and the spline they wrap
│   ├── bezier.rs                Bézier curves of both tiers
│   ├── trimmed.rs               TrimmedCurve and TrimmedCurve2D
│   ├── range.rs                 internal: RangeOperations, the box, discretization and
│   │                            B-spline of a range of each curve (ADR 0005)
│   ├── conversion.rs fitting.rs parametrize.rs   to B-spline, from points, onto a surface
├── surfaces/            the public surface types, one file per type
│   ├── parametric_surface.rs    ParametricSurface, ReversibleSurface, SurfaceCurvatures,
│   │                            normal recovery
│   ├── surface.rs               the Surface and ElementarySurface enums
│   ├── plane.rs … torus.rs bspline.rs bezier.rs
│   ├── extrusion.rs revolution.rs   the swept surfaces; swept.rs holds the
│   │                            profile geometry they share (ADR 0016)
│   ├── iso.rs                   helpers shared by the iso curves (u_iso, v_iso)
│   ├── pcurve.rs                pcurves to and from 3D curves (ADR 0030);
│   │                            pcurve/to_curve.rs is the pair table of
│   │                            pcurve_to_curve
│   └── conversion.rs fitting.rs
│
├── intersect/           intersections, by method (ADR 0011)
│   ├── routing.rs       Intersect: exact when a closed form exists, numeric otherwise
│   ├── routing/closed_form.rs  the closed-form pair table: which pairs have a closed
│   │                    form and which solver answers each (ADR 0023)
│   ├── routing/exact_family.rs  each typed exact result read as its family
│   │                    result (ADR 0018)
│   ├── routing/swept_plane.rs  plane cuts of the swept surfaces in closed form,
│   │                    entries of the table
│   ├── exact/           IntersectExact and the closed-form solvers, by family:
│   │                    curve_curve.rs, curve_surface.rs, surface_surface.rs
│   ├── numeric/         the searches for pairs with no closed form. subdivision.rs is
│   │                    the one driver of the searches on Bézier pieces (the queue of
│   │                    pairs, the budget of halvings, the pairs left unsettled);
│   │                    *_hull.rs and patch_pairs.rs bring their tests on a pair (ADR 0019).
│   │                    surface_patches.rs turns any surface into Bézier patches with a
│   │                    map back to its own parameters (ADR 0026).
│   │                    marching.rs composes marching/: the surface pair (pair.rs), the
│   │                    tracer's step loop (tracer.rs), the seeding policy (seeding.rs)
│   │                    and the coincidence check (coincidence.rs) (ADR 0025)
│   └── results/         the results Intersect returns, their regions (regions.rs),
│                        Completeness, InlineList
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

Every public item has one path, at the crate root. The folders are
private modules; `lib.rs` re-exports each public item (`geomcore::Line`,
never `geomcore::curves::Line`), and `prelude` re-exports the traits.
Moving a file between folders is therefore never a public change
(ADR 0022). Doc examples, tests, benches and the bindings use the flat
paths.

1. **One public type per file, named after the type.** A curve type is one
   generic type over its tier (ADR 0014): `GenericEllipse`, with its aliases
   `Ellipse` and `Ellipse2D`, lives in `curves/ellipse.rs`, and its 2D unit
   tests in `curves/ellipse/tests_2d.rs`. Only 2D types carry a dimension
   suffix (ADR 0013). A type's construction error lives beside it.
2. **Shared behaviour lives on the traits, implemented in the type's
   file.** `ParametricCurve`, `PointExtrema`, `ParametricCurve2D`,
   `PointExtrema2D`, `ParametricSurface`, `ReversibleCurve` and
   `ReversibleSurface` are defined in the `parametric_*.rs` files and
   `Transformable` in `transform.rs`, and each `impl` sits with its type.
   The enums forward them through their variant lists (ADR 0017), so no
   caller matches on a variant to move or reverse a shape (ADR 0024). A method
   that more than one type answers is a trait method, documented once on
   the trait; a type overrides a provided method only for a faster closed
   form, and its inherent `impl` keeps only what is specific to it
   (`Circle::radius`, `BSplineCurve::insert_knot`). Callers bring the
   traits in with `use geomcore::prelude::*` (ADR 0015). A generic curve
   type implements the space and plane traits from one body, with the
   crate-private `impl_for_both_tiers!` of `tier.rs`.
3. **Operation traits are implemented in their operation's folder.**
   `HasBoundingBox` impls live in `bounding_box/`, `Discretize` impls in
   `discretization/`, and `IntersectExact` impls in `intersect/exact/`, all
   grouped by family. To see everything an operation does, open one folder.
   The one exception is the table of the crate-private `RangeOperations`
   in `curves/range.rs`: each curve's box, discretization and B-spline
   over a range, in one impl per type that calls the formulas of those
   folders. Trimmed curves and the whole bounded curves read their
   operations from it (ADR 0005).
4. **Exact intersection solvers are grouped by pair family.**
   `Line`–`Plane` is in `intersect/exact/curve_surface.rs`, next to its
   result type `LinePlaneIntersection`. The general results that `Intersect`
   returns live in `intersect/results/`.
5. **Internal math sits under the shapes.** `elementary/`, `bspline/` and
   `math/` re-export nothing from `lib.rs`. They work on poles, knots and numbers,
   not on the public shape types; a shape calls into them, never the
   other way round.
6. **Unit tests live in a sibling `tests.rs`.** `curves/line.rs` ends with
   `#[cfg(test)] mod tests;` and its tests are in `curves/line/tests.rs`
   (for a `mod.rs`, the tests are in `tests.rs` beside it). Tests that only
   use the public API go in `tests/` instead.

## Where to start

**A new curve or surface type.**
- Add `curves/<name>.rs` (or `surfaces/<name>.rs`) holding the type, its
  constructor error and its shape-trait impls. A curve that exists in both
  tiers is one `Generic<Name><P: TierPoint>` with the aliases `<Name>` and
  `<Name>2D`.
- Re-export it from the folder's `mod.rs` and from `lib.rs`, its one
  public path.
- Add it as a variant of the `Curve`/`ElementaryCurve` (or `Surface`) enum,
  and as one line in that enum's `*_variants!` list beside it: every
  forwarding method and `From` impl follows from the list (ADR 0017), and
  the compiler points at the matches that do real per-variant work.
- Make it an operand of `Intersect`: a `CurveOperand`/`SurfaceOperand`
  impl in `intersect/routing.rs` saying which basis it is read as. The
  `Intersect` pair lists follow the variant list, so they need no entry
  (ADR 0023).
- Implement `RangeOperations` for it in `curves/range.rs`, calling a box
  formula from `bounding_box/`, a discretization from `discretization/`
  and a B-spline conversion from `curves/conversion.rs`. If it is
  bounded, add it to `box_over_bounds!` in `bounding_box/curves.rs` and
  to `discretize_bounded!` in `discretization/mod.rs`; otherwise its
  `bounding_box` is `None`.
- Add a sample of it to `tests/common/samples.rs`; `cargo test` names
  the variant until there is one, and every property test then covers it
  (ADR 0027).
- Bind it in Python: see [The Python bindings](#the-python-bindings).

**A new exact intersection pair.**
- Write the solver as `impl IntersectExact<B> for A` in the family file
  under `intersect/exact/`, with its result type, in both argument orders.
- Read its result into the family result with an `into_family` method
  in `intersect/routing/exact_family.rs`, and list the pair in its
  family's entry of the closed-form pair table,
  `intersect/routing/closed_form.rs`, so that `Intersect` routes the pair
  there instead of to the numeric search.

**A change to marching.** A new seed source, or a change to when one
runs, goes in `intersect/numeric/marching/seeding.rs`; a change to the
step loop, such as tracing through a point where the surfaces touch, in
`tracer.rs`; a new question asked of the two surfaces at a point on
both, in `pair.rs` (ADR 0025). A new surface kind joins the patch-pair
searches of marching and the curve searches by converting to patches in
`intersect/numeric/surface_patches.rs` (ADR 0026).

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
| `method_sets.rs` | the method sets the curve and surface classes share |
| `arguments.rs` | the class registry: the Python class of each core enum variant |
| `intersect.rs` | the routed `intersect` every class shares, and `curve_curve_extrema` |
| `results.rs` | the intersection, projection, extrema and curvature results |
| `errors.rs` | the Python exception hierarchy |

**Method sets.** A class writes its `#[pymethods]` block through
`shape_methods!`, naming the sets it takes (`curve`, `discretize`,
`surface`, ...) and listing only its own methods. Each set is the Python
face of a core trait (`curve` is `ParametricCurve`, `PointExtrema`,
`ReversibleCurve`, `Transformable`, `HasBoundingBox` and `Intersect`), with its docstrings and signatures
written once, so a class has a method exactly when its core type has it.
The table at the top of `method_sets.rs` says which set follows from
which trait.

The stubs mirror the sets (ADR 0029). Each set is a stub-only private
class in `curves.pyi` or `surfaces.pyi` (`_CurveMethods`,
`_DiscretizeMethods`, `_SurfaceMethods`, ...) holding the set's typed
signatures and its docstrings. A class stub inherits the sets its
`shape_methods!` line names and lists only its own methods. What differs
by class is a type argument: `_DiscretizeMethods[Discretization2D]`,
`_ParametrizeOnMethods[Ellipse2D]`, `_SurfaceMethods[Line, Circle]` (the
classes of the `u` and `v` iso curves).

**The class registry.** `shape_classes!` in `arguments.rs` lists, once per
core enum, the Python class of each variant (`Line => PyLine`). From that
list come the conversions both ways (`curve_classes::extract`,
`curve_classes::to_py`), the `TypeError` text naming the accepted
classes, and the classes' registration in their module.

**The interpreter lock.** `intersect` releases it around every core call,
so other Python threads run during a search (ADR 0021).

**Binding a new curve or surface type.**
- Add `geomcore-py/src/curves/<name>.rs` (or `surfaces/`) with the
  `#[pyclass]` holding the core type as field `0`, and its methods through
  `shape_methods!`: the sets its core traits give it, then its
  constructors and inherent methods.
- Add one `Variant => PyClass` line to the registry in `arguments.rs`.
  `intersect`, every argument that takes a curve or surface, and every
  method that returns one then accept the class.
- Add its stub to `geomcore-py/geomcore/curves.pyi` (or `surfaces.pyi`):
  a class inheriting the stub classes of the sets named in its
  `shape_methods!` line, with its own methods only, and add it to the
  `_WholeCurve`, `_Curve2D` or `_Surface` union there. stubtest checks
  the stub's methods against the class.
- `tests/test_parity.py` fails until the new class has a fixture. It then
  checks that the class has every method its stub inherits, and runs it
  through `intersect` against every other class.

## Integration tests

Each file in `tests/` is its own test crate. Its first line says what it
checks.

| Area | Files |
|---|---|
| Evaluation and construction | `combined_evaluation.rs`, `construction_validation.rs`, `precision.rs`, `surface_reversal.rs`, `transform_and_reversal.rs` |
| Golden fixtures | `vs_fixtures.rs`, `bezier_fixtures.rs`, `discretization_fixtures.rs`, `intersect_2d_fixtures.rs`, `swept_fixtures.rs` |
| B-splines and Bézier | `bspline_tiers.rs`, `degree_elevation.rs`, `bezier_surface.rs` |
| Projection | `projection_point.rs` |
| Intersection | `routed_intersections.rs`, `routed_surface_surface.rs`, `bezier_intersections.rs`, `bezier_through_enums.rs` |
| Bounding boxes | `bounding_box.rs`, `bounding_box_allocation.rs` |
| Discretization | `discretization.rs`, `discretization_allocation.rs` |
| Every variant | `shape_conformance.rs` |

`tests/common/samples.rs` holds samples of every variant of every shape
enum; the tests that check a property on every shape take their shapes
from it. `src/shape_enum/tests.rs` expands each enum's variant list over
it, so a new variant fails `cargo test`, naming itself, until it has a
sample there (ADR 0027). `tests/common/shapes.rs` decodes a fixture's
curve or surface from its `type` tag and replays its recorded samples,
inversions and projections through the traits.

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
