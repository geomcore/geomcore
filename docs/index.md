# geomcore

geomcore is a Rust library with Python bindings for building, evaluating
and querying 2D and 3D curves and 3D surfaces. Use it to find
intersections, measure clearance, fit measured points or map curves
between 3D and surface coordinates. Its scope is geometry. CAD topology,
B-rep and Boolean operations are outside that scope.

```viewer
scene: assets/scenes/getting-started.json
title: A circle and an ellipse, with the four points where they meet
height: 420
```

## Find what you need

| I want to… | Start here |
|---|---|
| Install and run my first example | [Getting started](getting-started.md) |
| Draw geometry in the browser | [Visualize with the viewer](viewer.md) |
| Build an arc, sweep a profile or place a shape | [Build and place geometry](guides/construction.md) |
| Fit measured points or a surface grid | [Elementary shapes](geometry/elementary.md#fitting-to-points), [B-spline fitting](nurbs/editing-fitting.md#fitting) |
| Measure clearance between shapes | [Closest approach](guides/projection.md#closest-approach-between-shapes) |
| Measure length, find a bounding box or evaluate a normal | [Length and bounds](nurbs/editing-fitting.md#measure-and-bound), [Local properties](python.md#local-properties-and-one-sided-evaluation) |
| Project a point or test containment | [Projection & containment](guides/projection.md) |
| Intersect curves or surfaces | [Intersection](guides/intersection.md) |
| Turn a curve into points and chords | [Discretize a curve](guides/discretization.md) |
| Map a curve onto a surface | [Map a curve onto a surface](guides/parametrization.md) |
| Look up a shape, formula or domain | [Elementary curves & surfaces](geometry/elementary.md) |
| Work with B-spline curves and surfaces | [NURBS curves](nurbs/curves.md), [surfaces](nurbs/surfaces.md), [editing & fitting](nurbs/editing-fitting.md) |
| Use the Python package | [Python API notes](python.md) |
| Understand why the design is what it is | [Design decisions](adr/index.md), [Why ten special cases](explanation/why-ten-special-cases.md) |
| See where the kernel is going | [Roadmap](roadmap.md) |

## What is in the box

- **Evaluation** of lines, circles, ellipses, parabolas, hyperbolas,
  B-spline and Bézier curves, in 3D and in 2D; and of planes, cylinders,
  cones, spheres, tori, B-spline and Bézier surfaces, and surfaces
  extruded or revolved from any curve.
- **Transformations** — translation, rotation, mirroring and uniform
  scaling through one `Transform` type, applied to any curve or surface.
- **Queries** — containment, projection and extrema with a tolerance;
  one `intersect` call that routes each pair to its exact closed form or
  to a numeric search and says which ran, and whether every hit was
  found.
- **Curves on surfaces** — the exact UV curve where one exists, and one
  B-spline within a deviation everywhere else, in both directions.
- **B-spline fitting and editing** — interpolation, least-squares
  approximation, knot insertion and removal, degree elevation,
  splitting and exact conversions from the elementary shapes.

The Rust API reference is on [docs.rs](https://docs.rs/geomcore). Python
does not have a generated API reference yet; its type stubs and
docstrings describe the full surface, and the
[Python API notes](python.md) explain the conventions that differ from
Rust.

## Status

Pre-1.0: the API can still change between minor versions, with breaking
changes flagged in the [changelog](https://github.com/geomcore/geomcore/blob/main/CHANGELOG.md).
Evaluation, transformations and curve-on-surface parametrization are
validated against golden fixtures at CAD-kernel tolerances (1e-7). See
the [roadmap](roadmap.md) for what comes next.
