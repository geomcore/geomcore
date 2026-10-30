# Golden fixtures

Each JSON file holds recorded inputs and the outputs geomcore must reproduce.
The integration tests in `tests/` load them through `tests/common/mod.rs`
(`load("name.json")`). Most files start with a `_note` key explaining their
conventions.

| File | Where the expected values come from | Replayed by |
|---|---|---|
| `construction.json` | Independent reference kernel | `vs_fixtures.rs` |
| `curves_analytic.json` | Independent reference kernel | `vs_fixtures.rs` |
| `curves_analytic_2d.json` | Independent reference kernel | `vs_fixtures.rs` |
| `curves_bspline.json` | Independent reference kernel; later cases also checked in exact rational arithmetic | `vs_fixtures.rs` |
| `parametrize.json` | Independent reference kernel | `vs_fixtures.rs` |
| `surfaces_analytic.json` | Independent reference kernel | `vs_fixtures.rs` |
| `surfaces_bspline.json` | Independent reference kernel | `vs_fixtures.rs` |
| `transforms.json` | Independent reference kernel | `vs_fixtures.rs` |
| `swept_surfaces.json` | Independent reference kernel | `swept_fixtures.rs` |
| `swept_plane_cuts.json` | Independent reference kernel | `swept_fixtures.rs` |
| `bezier.json` | Explicit Bernstein sums, computed outside geomcore | `bezier_fixtures.rs` |
| `discretization.json` | geomcore's own output, frozen to catch regressions | `discretization_fixtures.rs` |
| `intersect_2d.json` | Derived by hand; each case's `note` gives the derivation | `intersect_2d_fixtures.rs` |

## Updating a fixture

- **Reference-kernel files.** These are recorded by a generator script that
  is kept outside this repository. Don't edit their values by hand. If a
  case needs to change, ask a maintainer to re-record it.
- **`discretization.json`.** Regenerate with
  `UPDATE_DISCRETIZATION_FIXTURES=1 cargo test --test discretization_fixtures`,
  then review the diff: every change is a change in geomcore's output.
- **`bezier.json` and `intersect_2d.json`.** Edit by hand, and keep the
  derivation in the case's note.
