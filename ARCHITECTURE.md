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
│
├── base/                the values every shape is built from
│   ├── point.rs vector.rs   Point/Point2D, Vector/Vector2D, with the predicates the
│   │                        solvers ask (is_parallel, is_perpendicular,
│   │                        component_perpendicular_to)
│   ├── direction.rs         Direction/Direction2D: unit vectors by type (ADR 0054),
│   │                        whose is_parallel and is_perpendicular skip normalizing
│   ├── frame.rs             Axis, Frame and their 2D forms: where a shape is placed;
│   │                        Axis::is_coaxial
│   ├── transform.rs         Transform, and the Transformable trait of the shapes it moves
│   ├── tolerance.rs         Tolerance: every geometric decision takes one
│   ├── tier.rs              TierPoint: the two tiers (space, plane) the curve types are
│   │                        generic over
│   ├── pole.rs              Pole, HomogeneousPole, PolePoint: the poles a B-spline
│   │                        blends, points or weighted points, for both tiers
│   ├── side.rs continuity.rs   Side and Continuity: the one-sided queries at a knot
│   │                        or a seam, and how smooth a shape is there
│   └── inline_list.rs       InlineList: up to a proven bound of items, without allocating
│
├── geometry/            the curve and surface types
│   ├── curves/              the public curve types, one file per type
│   │   ├── parametric_curve.rs      ParametricCurve, PointExtrema, ReversibleCurve
│   │   ├── parametric_curve_2d.rs   ParametricCurve2D, PointExtrema2D
│   │   ├── curve.rs curve_2d.rs     the Curve, ElementaryCurve, Curve2D enums
│   │   ├── line.rs … hyperbola.rs   elementary curves, each one generic type for both tiers
│   │   ├── bspline.rs               B-spline curves of both tiers
│   │   ├── bezier.rs                Bézier curves of both tiers
│   │   ├── trimmed.rs               GenericTrimmedCurve: TrimmedCurve and TrimmedCurve2D
│   │   ├── range.rs                 internal: RangeOperations, the box, discretization and
│   │   │                            B-spline of a range of each curve (ADR 0005)
│   │   ├── conversion.rs fitting.rs   to B-spline, from points
│   │   └── (QuarticBranch)          a Curve variant that lives beside the quartic it is
│   │                                a piece of: operations/intersect/exact/quartic/branch.rs
│   ├── surfaces/            the public surface types, one file per type
│   │   ├── parametric_surface.rs    ParametricSurface, ReversibleSurface, SurfaceCurvatures,
│   │   │                            normal recovery
│   │   ├── surface.rs               the Surface and ElementarySurface enums
│   │   ├── plane.rs … torus.rs bspline.rs bezier.rs
│   │   ├── extrusion.rs revolution.rs   the swept surfaces; swept.rs holds the
│   │   │                            profile geometry they share (ADR 0016)
│   │   ├── trimmed.rs               TrimmedSurface: a basis surface over a u range and
│   │   │                            a v range (ADR 0039)
│   │   ├── iso.rs                   helpers shared by the iso curves (u_iso, v_iso)
│   │   ├── uv_curve.rs              UV curves to and from 3D curves (ADR 0030,
│   │   │                            0031); uv_curve/to_curve.rs is the pair table
│   │   │                            of uv_curve_to_3d, uv_curve/composition.rs its
│   │   │                            exact composition on a B-spline surface (ADR
│   │   │                            0032), uv_curve/from_curve.rs the closed forms
│   │   │                            of curve_3d_to_uv (from_curve/swept.rs those
│   │   │                            along the iso-curves of the swept surfaces),
│   │   │                            uv_curve/approximation.rs
│   │   │                            uv_curve_to_3d_within (ADR 0035) and
│   │   │                            uv_curve/from_curve_within.rs its mirror
│   │   │                            curve_3d_to_uv_within
│   │   └── conversion.rs fitting.rs
│   ├── approximation.rs         the results of the operations within a deviation:
│   │                            CurveApproximation, SurfaceApproximation
│   ├── fit.rs                   the least-squares fits of planes, circles, spheres
│   │                            and cylinders to points: Fit and FitError
│   └── enum_variants.rs         the macros that expand a shape enum's variant list
│                                (ADR 0017)
│
├── operations/          the operations on curves and surfaces, one folder each
│   ├── discretization/      Discretize: proven point sequences on bounded curves;
│   │                        sampled.rs measures the curves with no closed-form
│   │                        deviations (a quartic branch)
│   ├── extrema/             distance extrema between two curves (ADR 0033)
│   │   ├── curve_curve.rs   curve_curve_extrema: seeds, refinement and its results
│   │   └── curve_curve_pieces.rs  the minima search on the Bézier pieces of two
│   │                        B-splines
│   ├── intersect/           intersections, by method (ADR 0011)
│   │   ├── routing.rs       Intersect: exact when a closed form exists, numeric otherwise
│   │   ├── routing/closed_form.rs  the closed-form pair table: which pairs have a closed
│   │   │                    form and which solver answers each (ADR 0023)
│   │   ├── routing/exact_family.rs  each typed exact result read as its family
│   │   │                    result (ADR 0018)
│   │   ├── routing/swept_plane.rs  plane cuts of the swept surfaces in closed form,
│   │   │                    entries of the table
│   │   ├── routing/trimmed_surface.rs  the intersection of trimmed surfaces, with a
│   │   │                    surface or a curve: the bases intersected, the result
│   │   │                    clipped to the ranges (ADR 0039)
│   │   ├── routing/trimmed_surface/shared_boundary.rs  what two trimmed surfaces of
│   │   │                    one surface share: an overlap with area, or the edges
│   │   │                    and corners they touch along (#406)
│   │   ├── exact/           IntersectExact and the closed-form solvers, by family:
│   │   │                    curve_curve.rs, curve_surface.rs, surface_surface.rs.
│   │   │                    conic.rs and conic_2d.rs solve a conic against a quadric
│   │   │                    or another conic for them. quartic.rs is the quartic
│   │   │                    curve where two quadrics meet, read along the generators
│   │   │                    of the ruled one and fitted within confusion (ADR 0036,
│   │   │                    0040, 0041); quartic/branch.rs is QuarticBranch, its
│   │   │                    pieces that run to infinity (ADR 0042)
│   │   ├── numeric/         the searches for pairs with no closed form. subdivision.rs is
│   │   │                    the one driver of the searches on Bézier pieces (the queue of
│   │   │                    pairs, the budget of halvings, the pairs left unsettled);
│   │   │                    *_hull.rs and patch_pairs.rs bring their tests on a pair (ADR 0019);
│   │   │                    bezier_pieces.rs holds the Bézier pieces and the tests on
│   │   │                    pairs of them that the searches share;
│   │   │                    pairs_with_bspline.rs runs the searches on the pairs holding
│   │   │                    a B-spline, the other shape converted to its exact B-spline.
│   │   │                    curve_curve.rs reads its hits from the extrema at zero distance;
│   │   │                    curve_surface.rs is the curve-surface search (sampled seeds
│   │   │                    refined, or piece pairs on B-splines); overlap_runs.rs joins
│   │   │                    the runs along which a curve lies on the other operand;
│   │   │                    embedding_2d.rs meets 2D curves through the 3D search, at z = 0.
│   │   │                    parameter_domain.rs is the bounds or period of a curve or a
│   │   │                    surface direction, with the searches' wrapping and slack.
│   │   │                    surface_patches.rs turns any surface into Bézier patches with a
│   │   │                    map back to its own parameters (ADR 0026).
│   │   │                    marching.rs composes marching/: the surface pair (pair.rs), the
│   │   │                    tracer's step loop (tracer.rs), the seeding policy (seeding.rs,
│   │   │                    with the seeds on closed loops of level_extrema.rs)
│   │   │                    and the coincidence check (coincidence.rs) (ADR 0025)
│   │   └── results/         the results Intersect returns, their regions (regions.rs),
│   │                        Completeness, InlineOrHeapList
│   └── projection.rs    results of point projection (CurveProjection, …)
│
├── bounding_box/        HasBoundingBox: the box type and one file per family; also the
│                        culling box of the numeric searches and the tiers. tight.rs
│                        holds the tight boxes of B-spline shapes
│
└── math/                internal math under the curves and surfaces
    ├── elementary/      closed-form formulas of the elementary curves and surfaces,
    │                    and fitting.rs, the numbers under their least-squares fits
    ├── bspline/         the B-spline engine (curve evaluation, editing,
    │                    tensor-product surfaces, exact conic nets) and errors.rs,
    │                    the errors its checks raise
    ├── numeric/         polynomial roots, Bernstein form, solvers, quadrature,
    │                    symmetric eigenvectors, Levenberg–Marquardt steps
    └── parameter_range.rs  the rules of a parameter range: the shift by whole
                         periods, the slack at the ends, the nearest end across
                         the seam; the trimmed curves and surfaces and the
                         intersection router read them here
```

Outside `src/`:

| Path | Holds |
|---|---|
| `tests/` | Integration tests through the public API, one file per area (see below) |
| `tests/fixtures/` | Golden JSON data; [its README](tests/fixtures/README.md) says where each file comes from |
| `benches/` | Criterion benchmarks; results are discussed in `docs/benchmarks.md` |
| `examples/` | Small Rust programs (`cargo run --example sample_conics`) |
| `geomcore-py/` | The Python bindings: `src/` (PyO3), `geomcore/` (stubs), `tests/`, `examples/` |
| `docs/` | The user documentation site (`mkdocs.yml`), organised by Diátaxis mode: tutorials, how-to guides, reference and explanation, plus the `adr/` records |
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
   `Ellipse` and `Ellipse2D`, lives in `geometry/curves/ellipse.rs`, and its 2D unit
   tests in `geometry/curves/ellipse/tests_2d.rs`. Only 2D types carry a dimension
   suffix (ADR 0013). An error lives with the code that checks it: a type's
   construction error lives beside the type, and the errors of the B-spline
   engine's checks (the knot vector, knot insertion and removal, degree
   elevation, splits) live in `math/bspline/errors.rs`, beside those checks.
   One curve type lives elsewhere. `QuarticBranch`, a `Curve` variant, is
   built from the two quadrics it lies on and read through the quartic's
   walk, so it sits beside that quartic in
   `operations/intersect/exact/quartic/branch.rs`, and `geometry/curves.rs`
   re-exports it (ADR 0042). Through it the `Curve` enum reaches the cone
   and the cylinder. The check that curves do not depend on surfaces
   (`tests/module_dependencies.rs`, ADR 0031) reads the files under
   `geometry/curves/` only, so it does not see this.
2. **Shared behaviour lives on the traits, implemented in the type's
   file.** `ParametricCurve`, `PointExtrema`, `ParametricCurve2D`,
   `PointExtrema2D`, `ParametricSurface`, `ReversibleCurve` and
   `ReversibleSurface` are defined in the `parametric_*.rs` files and
   `Transformable` in `base/transform.rs`, and each `impl` sits with its type.
   `PointExtrema` is a curve's own inversion and projection, as
   `parameters_of` and `project_point` are a surface's, not an operation
   (ADR 0033).
   The enums forward them through their variant lists (ADR 0017), so no
   caller matches on a variant to move or reverse a shape (ADR 0024). A method
   that more than one type answers is a trait method, documented once on
   the trait; a type overrides a provided method only for a faster closed
   form, and its inherent `impl` keeps only what is specific to it
   (`Circle::radius`, `BSplineCurve::insert_knot`). Callers bring the
   traits in with `use geomcore::prelude::*` (ADR 0015). A generic curve
   type implements the space and plane traits from one body, with the
   crate-private `impl_for_both_tiers!` of `base/tier.rs`.
3. **Operation traits are implemented in their operation's folder.**
   `HasBoundingBox` impls live in `bounding_box/`, `Discretize` impls in
   `operations/discretization/`, and `IntersectExact` impls in `operations/intersect/exact/`, all
   grouped by family. To see everything an operation does, open one folder.
   The one exception is the table of the crate-private `RangeOperations`
   in `geometry/curves/range.rs`: each curve's box, discretization and B-spline
   over a range, in one impl per type that calls the formulas of those
   folders. Trimmed curves and the whole bounded curves read their
   operations from it (ADR 0005).
4. **Exact intersection solvers are grouped by pair family.**
   `Line`–`Plane` is in `operations/intersect/exact/curve_surface.rs`, next to its
   result type `LinePlaneIntersection`. The general results that `Intersect`
   returns live in `operations/intersect/results/`.
5. **Internal math sits under the shapes.** `math/` (`elementary/`,
   `bspline/` and `numeric/`) works on poles, knots, numbers and the values
   of `base/`, not on the public shape types; a shape calls into it, never
   the other way round. No file under `math/` names `geometry`,
   `operations` or `bounding_box`, nor any name the crate root re-exports
   from them, and `tests/module_dependencies.rs` checks it. A value both
   sides need, such as `Side`, `Continuity` or `InlineList`, lives in
   `base/`. The only public items `math/` defines are the errors of the
   B-spline engine's checks (rule 1). `base/` sits under `math/` in turn:
   the pole algebra the engine runs on, which `TierPoint` requires, lives
   in `base/pole.rs`, and no file under `base/` names `math`, which
   `tests/module_dependencies.rs` also checks.
6. **Unit tests live in a sibling `tests.rs`.** `geometry/curves/line.rs` ends with
   `#[cfg(test)] mod tests;` and its tests are in `geometry/curves/line/tests.rs`.
   Tests that only
   use the public API go in `tests/` instead.
7. **A geometric question has one home.** Whether two directions are
   parallel or perpendicular, whether two axes are the same line, and the
   part of a vector across a direction are answered by
   `Direction::is_parallel` and `Direction::is_perpendicular` between
   directions, `Vector::is_parallel` and `Vector::is_perpendicular` between
   vectors of any length, `Vector::component_perpendicular_to`, and
   `Axis::is_parallel`, `Axis::contains` and `Axis::is_coaxial`, never written
   out at the call site (`tests/one_home_per_predicate.rs` checks it). The
   caller chooses only the angle: `tol.angular`
   for unbounded geometry, `tol.angular_for_extent(size)` where the
   geometry is bounded (ADR 0006). A vector the caller knows to be unit
   becomes a direction once, through `Direction::from_unit_vector`, which
   debug builds check (ADR 0054). Tier-generic code reaches the predicates
   through the tier's `VectorOperations`.

## Where to start

**A new curve or surface type.**
- Add `geometry/curves/<name>.rs` (or `geometry/surfaces/<name>.rs`) holding the type, its
  constructor error and its shape-trait impls. A curve that exists in both
  tiers is one `Generic<Name><P: TierPoint>` with the aliases `<Name>` and
  `<Name>2D`.
- Re-export it from the folder's module file (`curves.rs`, `surfaces.rs`) and from `lib.rs`, its one
  public path.
- Add it as a variant of the `Curve`/`ElementaryCurve` (or `Surface`) enum,
  and as one line in that enum's `*_variants!` list beside it: every
  forwarding method and `From` impl follows from the list (ADR 0017), and
  the compiler points at the matches that do real per-variant work.
- Make it an operand of `Intersect`: a `CurveOperand`/`SurfaceOperand`
  impl in `operations/intersect/routing.rs` saying which basis it is read as. The
  `Intersect` pair lists follow the variant list, so they need no entry
  (ADR 0023).
- Implement `RangeOperations` for it in `geometry/curves/range.rs`, calling a box
  formula from `bounding_box/`, a discretization from `operations/discretization/`
  and a B-spline conversion from `geometry/curves/conversion.rs`. If it is
  bounded, add it to `box_over_bounds!` in `bounding_box/curves.rs` and
  to `discretize_bounded!` in `operations/discretization.rs`; otherwise its
  `bounding_box` is `None`.
- Add a sample of it to `tests/common/samples.rs`; `cargo test` names
  the variant until there is one, and every property test then covers it
  (ADR 0027).
- Bind it in Python: see [The Python bindings](#the-python-bindings).

**A new exact intersection pair.**
- Write the solver as `impl IntersectExact<B> for A` in the family file
  under `operations/intersect/exact/`, with its result type, in both argument orders.
- Read its result into the family result with an `into_family` method
  in `operations/intersect/routing/exact_family.rs`, and list the pair in its
  family's entry of the closed-form pair table,
  `operations/intersect/routing/closed_form.rs`, so that `Intersect` routes the pair
  there instead of to the numeric search.
- Test its configurations in Rust, beside the solver. The bindings need
  no new test unless the pair returns a kind of result they have not
  converted before (see [The Python bindings](#the-python-bindings)).

**A change to marching.** A new seed source, or a change to when one
runs, goes in `operations/intersect/numeric/marching/seeding.rs`; a change to the
step loop, such as tracing through a point where the surfaces touch, in
`tracer.rs`; a new question asked of the two surfaces at a point on
both, in `pair.rs` (ADR 0025). A new surface kind joins the patch-pair
searches of marching and the curve searches by converting to patches in
`operations/intersect/numeric/surface_patches.rs` (ADR 0026).

**A new extrema query.** A new pair (curve–surface, surface–surface)
goes in its own file under `operations/extrema/`, next to its result
type, re-exported from `extrema.rs` and `lib.rs`. A numeric search on
B-spline pieces reuses the driver in
`operations/intersect/numeric/subdivision.rs`, as `curve_curve_pieces.rs`
does (ADR 0019, ADR 0033).

**A new operation.**
- Give it its own folder under `operations/`, like
  `operations/discretization/`: the trait in the module file beside the
  folder (`operations/discretization.rs`), the impls by family in the
  folder.
- Record the term in `CONTEXT.md`.

## The Python bindings

`geomcore-py/src/` uses the core crate's folder names. Each class sits in
the file its core type has, under `curves/` and `surfaces/`. The other
files are:

| File | Holds |
|---|---|
| `lib.rs` | the module wiring (`#[pymodule]`) and small shared argument checks |
| `values.rs` | `Point`, `Vector`, `Axis`, `Frame` and their 2D forms, `Tolerance`, `Transform` |
| `enums.rs` | the enumerations Python sees as `enum` classes: `Side`, `PointParametrization`, `IntersectionMethod`, `CompletenessStatus`, `IncompleteReason` |
| `bounding_box.rs`, `discretization.rs` | the classes of those operations |
| `approximation.rs` | `CurveApproximation`, `CurveApproximation2D`, `SurfaceApproximation`: what the operations within a deviation return |
| `fit.rs` | `PlaneFit`, `CircleFit`, `CircleFit2D`, `SphereFit`, `CylinderFit`: what the least-squares fits return |
| `arrays.rs` | the numpy bridge: the bulk queries that read and return float arrays |
| `pickling.rs` | pickling and the coordinate protocols the value classes share |
| `method_sets.rs` | the method sets the curve and surface classes share |
| `arguments.rs` | the class registry: the Python class of each core enum variant |
| `intersect.rs` | the routed `intersect` every class shares, and `curve_curve_extrema` |
| `results.rs` | the intersection, projection, extrema and curvature results |
| `errors.rs` | the Python exception hierarchy |
| `surfaces/uv_curve.rs` | the argument handling of `uv_curve_to_3d`, and the deprecated free functions of the UV curve mappings |

**Method sets.** A class writes its `#[pymethods]` block through
`shape_methods!`, naming the sets it takes (`curve`, `discretize`,
`surface`, ...) and listing only its own methods. Each set is the Python
face of a core trait (`curve` is `ParametricCurve`, `PointExtrema`,
`ReversibleCurve`, `Transformable`, `HasBoundingBox` and `Intersect`), with its docstrings and signatures
written once, so a class has a method exactly when its core type has it.
The surface sets split the same way: `parametric_surface` is what every
`ParametricSurface` with a box answers, `TrimmedSurface` included, and
`surface` what only a whole surface does (reversal, transform, iso curves,
`trimmed`, `intersect` with curves). `TrimmedSurface` writes its own
`transformed`, `trimmed` and `intersect`, whose docstrings say what
becomes of its ranges.
The table at the top of `method_sets.rs` says which set follows from
which trait.

The stubs mirror the sets (ADR 0029). Each set is a stub-only private
class in `geomcore/__init__.pyi` (`_CurveMethods`,
`_BoundedCurveMethods`, `_SurfaceMethods`, ...) holding the set's typed
signatures and its docstrings. A class stub inherits the sets its
`shape_methods!` line names and lists only its own methods. What differs
by class is a type argument: `_BoundedCurveMethods[Point2D, Discretization2D]`,
`_SurfaceMethods[Line, Circle]` (the classes of the `u` and `v` iso
curves).

**The class registry.** `shape_classes!` in `arguments.rs` lists, once per
core enum, the Python class of each variant (`Line => PyLine`). From that
list come the conversions both ways (`curve_classes::extract`,
`curve_classes::to_py`), the `TypeError` text naming the accepted
classes, and the classes' registration in their module.

**The interpreter lock.** `intersect` releases it around every core call,
so other Python threads run during a search (ADR 0021).

**What a Python test covers.** A test under `geomcore-py/tests/` covers
what the binding adds: the conversion of arguments and results, the
operands a method accepts, `repr`, pickling, the errors it raises, and a
new class or result field. Which configuration gives which curves,
points or overlaps is asserted in the core tests, beside the solver.
`intersect` has no path per pair: it calls the core router and converts
the family result (ADR 0018, ADR 0021), and `tests/test_parity.py` runs
every class against every other through it. A new closed-form
configuration therefore needs no Python test. The Python tests keep one
case per kind of result the binding converts (a contact point, an exact
conic, a fitted quartic, a quartic branch, a coincident pair, a marched
curve, an overlap, ...), so each conversion stays exercised end to end.

**Binding a new curve or surface type.**
- Add `geomcore-py/src/curves/<name>.rs` (or `surfaces/`) with the
  `#[pyclass]` holding the core type as field `0`, and its methods through
  `shape_methods!`: the sets its core traits give it, then its
  constructors and inherent methods. Declare the file and re-export the
  class in `curves/mod.rs` (or `surfaces/mod.rs`).
- Add one `Variant => PyClass` line to the registry in `arguments.rs`.
  `intersect`, every argument that takes a curve or surface, and every
  method that returns one then accept the class.
- Add its stub to `geomcore-py/geomcore/__init__.pyi`: a class
  inheriting the stub classes of the sets named in its `shape_methods!`
  line, with its own methods only. Add it to the `Curve`, `Curve2D` or
  `Surface` union there and in `__init__.py`, and to the stub's
  `__all__`. stubtest checks the stub's methods against the class.
- `tests/test_parity.py` fails until the new class has a fixture. It then
  checks that the class has every method its stub inherits, and runs it
  through `intersect` against every other class.

## Integration tests

Each file in `tests/` is its own test crate. Its first line says what it
checks.

| Area | Files |
|---|---|
| Evaluation and construction | `combined_evaluation.rs`, `construction_validation.rs`, `out_of_domain_evaluation.rs`, `precision.rs`, `surface_reversal.rs`, `transform_and_reversal.rs` |
| Golden fixtures | `vs_fixtures.rs`, `bezier_fixtures.rs`, `discretization_fixtures.rs`, `intersect_2d_fixtures.rs`, `swept_fixtures.rs` |
| B-splines and Bézier | `bspline_construction_and_shape.rs`, `bspline_tiers.rs`, `degree_elevation.rs`, `bezier_surface.rs` |
| Projection | `projection_point.rs` |
| Fitting | `least_squares_fits.rs` |
| UV curves | `uv_curve_to_3d.rs`, `uv_curve_composition.rs`, `uv_curve_to_3d_within.rs` |
| Intersection | `routed_intersections.rs`, `routed_surface_surface.rs`, `trimmed_surfaces.rs`, `bezier_intersections.rs`, `bezier_through_enums.rs`, `intersection_fingerprints.rs` |
| Bounding boxes | `bounding_box.rs`, `bounding_box_allocation.rs` |
| Discretization | `discretization.rs`, `discretization_allocation.rs` |
| Every variant | `shape_conformance.rs` |
| Module rules | `module_dependencies.rs` |

`tests/common/samples.rs` holds samples of every variant of every shape
enum; the tests that check a property on every shape take their shapes
from it. `src/geometry/enum_variants/tests.rs` expands each enum's variant list over
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
