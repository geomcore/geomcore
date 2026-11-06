# The Python package is flat, like the crate root

Every public Python name is reached by one path, at the package root: `geomcore.Circle`, `geomcore.BSplineSurface`, `geomcore.Point`, `geomcore.CurveCurveIntersection` (issue #203). This is ADR 0022 for the bindings. Before this decision the shapes lived in `geomcore.curves` and `geomcore.surfaces`. Values, results, enums and errors were at the root, and the native module was reachable as `geomcore.geomcore`, a second path to every name that IDE completion offered next to the first. A user had to know which of three modules held a name, and the UX study's first script imported `Plane` from `geomcore` and failed.

- The extension module is `geomcore._native`. `geomcore/__init__.py` re-exports it, so the native module's name is private and no longer doubles every path. Every class's `__module__` is `geomcore`, which is the path pickle records.
- `Curve`, `Curve2D` and `Surface` are public: each is the union of the classes of the Rust enum of the same name, and so their Python counterpart. In the stubs it is a type alias to annotate with. At runtime it is a `typing.Union` of the same classes, which `isinstance` takes from Python 3.10. On 3.9, `typing.get_args(Curve)` gives the tuple.
- `TrimmedCurve` is not a `Curve`, as in Rust, where it is a piece of one and not an enum variant. A method that takes a whole curve or a piece is typed `Curve | TrimmedCurve`. One that takes or returns a whole curve only, such as a trimmed curve's `basis` or a swept surface's profile, says `Curve`. The special case is spelled where it applies, rather than every whole-curve annotation needing a second, private name.
- `geomcore.curves` and `geomcore.surfaces` stay for one release as plain Python modules. Each name they held resolves through a module `__getattr__` with a `DeprecationWarning`, so code and pickles written against them keep working.
- The stubs are one file, `geomcore/__init__.pyi`. `curves.pyi` and `surfaces.pyi` only re-export, for type checkers reading the deprecated imports.

The root holds about 90 names. Like the crate root's 150 (ADR 0022), they do not collide, and a name says its kind (`…Intersection`, `…Error`, the `2D` suffix of ADR 0013).

## Considered options

- **Keep the split and document a module map**: rejected. Three modules to learn, against one in Rust, and the shapes were the only names a user had to look up.
- **Re-export the shapes at the root and keep the submodules for good**: rejected. Two paths per shape, the doubling ADR 0022 removed from Rust.
- **Runtime `Curve` as a tuple of classes**: rejected. It works with `isinstance` on every version, but stubtest and mypy see a tuple where the stub has a union, so the one name could not serve both uses.
- **`Curve` including `TrimmedCurve`, with a private `_WholeCurve` for the rest**: rejected. It made every whole-curve annotation (a trimmed curve's basis, a swept profile, the result of `uv_curve_to_3d`) use a private name in the public stubs, and the Python `Curve` differed from the Rust enum.
- **An abstract base class registering each class**: rejected. The stub classes would have to inherit it, and its methods would be typed twice, as the sets already are (ADR 0029).

## Consequences

- Breaking in the next release: `geomcore.curves` and `geomcore.surfaces` are removed, and `geomcore.geomcore` is gone now. It was never documented.
- `isinstance(x, geomcore.Curve)` raises on Python 3.9. `typing.get_args` is the portable form.
