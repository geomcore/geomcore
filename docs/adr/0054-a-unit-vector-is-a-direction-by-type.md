# A unit vector is a direction, by type

The directions of the kernel are unit vectors: an axis's direction, a frame's three axes, a circle's normal, a line's direction, a surface's normal, a curve's tangent and a surface's principal directions. They are values of their own type, `Direction` in space and `Direction2D` in the plane, which hold a vector of unit length by construction. `Vector::direction` builds one, normalizing once and returning `None` for a zero vector. The type's field is private, so no code outside `base/` can make a direction that is not unit.

The type carries a claim that the solvers otherwise make at each call site. Two directions are parallel when the sine of the angle between them, `|a × b|`, is at most `angular`, and perpendicular when its cosine, `|a · b|`, is. `Vector::is_parallel` must first divide by `|a| |b|`, a square root and a division on every call. The solvers ask these questions in their inner loops, and the instruction-count budget showed that paying for the normalization costs several percent. Before this type, the fast forms took plain vectors and relied on each caller to pass unit ones. A debug assertion checked the claim, at every call site separately. With the type, the claim is checked once, where the direction is built.

A direction reads as the unit vector it holds: it dereferences to `Vector`, so `d.x`, `d.dot(v)`, `d.cross(v)` and the other read-only operations are those of the vector and are not written twice. Its arithmetic gives vectors, `d * length`, `d1 + d2`, since a scaled or summed direction is not unit. Its cross product is a vector for the same reason: two directions' cross product is unit only when they are perpendicular, and the caller decides whether to make it a direction. Negation keeps a direction. A vector taken as an argument, by a constructor or by `dot`, `cross`, `angle`, `is_parallel` and `is_perpendicular`, may be a direction (`impl Into<Vector>`). `Vector::component_along` and `component_perpendicular_to` take only a direction, which is their precondition.

Inside the crate, a vector known to be unit by its construction (a speed divided by its length, a normal of a parametric surface) becomes a direction through the crate-private `from_unit_vector`, which checks the claim in debug builds, so the test suite checks every such site.

## Considered options

- **Plain vectors with crate-private unit forms of the predicates** (`is_parallel_unit` and the like, debug-checked): rejected. It keeps the claim at every call site, where a wrong one is found only by the tests that reach it, and public callers cannot use the fast forms at all.
- **A direction that duplicates the vector's operations**, as a separate type with its own `dot`, `cross` and components: rejected. It doubles the surface without adding a guarantee: reading a direction cannot make it non-unit.
- **A direction whose cross product renormalizes**: rejected. It pays a square root on every cross product, and has to fail when the directions are parallel.
- **A common trait of vectors and directions**: rejected. No code is generic over the two; the conversion `Into<Vector>` covers taking either.
- **A `Direction` class in Python**: not now. Python has no static types to carry the claim through, so a direction comes back as a unit `Vector`, and `component_along` there still takes a vector of any length.

## Consequences

The accessors that return a direction break code that bound their result to a `Vector`: it writes `d.vector()`. `ParametricSurface::normal`, the curves' `tangent` and the principal directions of `SurfaceCurvatures` return directions too (#401), so a solver asks a unit question of them without wrapping them. The solvers' own helpers that may return a gradient, or a NaN on an axis, stay `Vector`.
