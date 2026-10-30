# Each public item has one path, at the crate root

Every public type, trait and function is reached by one path: the crate root, `geomcore::Line`, `geomcore::BSplineSurface`, `geomcore::CurveSurfaceIntersection`. The folders of the source (`curves/`, `surfaces/`, `intersect/`, `bounding_box/`, `discretization/`, and the value modules `point`, `vector`, `frame`, `transform`, `tolerance`, `projection`, `tier`) are private modules. The one public module is `prelude`, which re-exports the shape and operation traits for a glob import (ADR 0015); its re-exports are not inlined in the documentation, so each trait has one page, at the root.

Before this decision every item had two paths, `geomcore::Line` and `geomcore::curves::Line`. The flat one was used about six times as often in the repository, yet doc examples mixed both in one block, so a reader could not tell which was canonical. Worse, the folder paths made the source layout public: moving a file between folders, as the first layout pass did, was a breaking change. With the folders private, the layout is an implementation detail, and the crate root is the whole interface.

No folder is kept public as a namespace. The kernel has about 150 public names, the shape names are distinct across folders (only 2D types carry a suffix, ADR 0013), and a user never needs a folder to disambiguate. The intersection results are the largest group, and their names already say their pair (`LinePlaneIntersection`).

`ParametricCurve::as_bspline`, `ParametricCurve::as_elementary`, `ParametricCurve2D::as_bspline`, `ParametricSurface::as_bspline` and `ParametricSurface::as_elementary` stay on the public traits, hidden from the documentation, as routing hooks of the numeric search. The search is generic over any `ParametricCurve` and `ParametricSurface`: `intersect_curve_curve` and `intersect_curve_surface` accept a curve of the caller's own (the helix in `tests/routed_intersections.rs`), and `TrimmedCurve<C>` wraps any curve. Inside, the search must ask a `T: ParametricCurve` whether it is a B-spline, to take the Bézier-piece path that proves its hits complete. Without specialization, only a method of the trait can answer that for every `T`. Wrappers outside the crate must forward the hooks too: the Python bindings wrap a curve or a trimmed curve in one enum and pass it to the search.

## Considered options

- **Keep the folder modules public as well**: rejected. Two paths per item, and every file move is a breaking change.
- **Folders public, root flat names removed**: rejected. The flat names are the ones callers use, and the folders would turn the layout into the interface.
- **Move the hooks to a crate-private trait implemented by the crate's own types**: rejected. A caller's curve, and a `TrimmedCurve` over one, could not implement it, so the public generic search would need a blanket impl for every `ParametricCurve`, which conflicts with the crate's own impls without specialization. Wrapping the caller's curve in an adapter that answers `None` would lose the B-spline path for a caller passing `TrimmedCurve<&BSplineCurve>` and downgrade its completeness.
- **Downcast through `Any`**: rejected. It needs `'static` operands, and a trimmed curve borrowing its basis is not.

## Consequences

- Breaking: the `geomcore::curves::…`, `geomcore::surfaces::…`, `geomcore::intersect::…` and the other folder paths are gone; each item is at `geomcore::…` under the same name. Rust offers no cheap deprecation for a module path, so the change ships in one release.
- Moving a file between folders is never a public change again.
- The hooks remain callable, as hidden items with no stability promise. A future specialization, or a sealed curve trait, could remove them.
