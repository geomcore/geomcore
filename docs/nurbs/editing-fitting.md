# Editing & fitting

Reshape B-splines without leaving the representation (editing), build
them from points (fitting), convert exact geometry into them, and
measure them. Curve API first, with per-direction surface equivalents
inline. The Python bindings mirror the Rust names 1:1, with errors as
`ValueError`.

## Knots and degree

Insertion, refinement, elevation and removal keep the shape exactly
(rational poles are edited in homogeneous form). Removal reports how
many copies actually came out — a knot whose removal would move the
curve past `max_deviation` is kept.

=== "Python"

    ```python
    refined = curve.insert_knot(0.5, 1)
    assert refined.eval_point(0.5) == curve.eval_point(0.5)

    dense = curve.refine_knots([0.25, 0.5, 0.75])

    elevated = curve.elevate_degree(1)
    assert elevated.degree() == 3

    simpler, removed = dense.remove_knot(0.5, 1, 1e-9)
    assert (removed, simpler.knots()) == (1, [0.0, 0.25, 0.75, 1.0])
    ```

=== "Rust"

    ```rust
    let refined = curve.insert_knot(0.5, 1).unwrap();
    assert_eq!(refined.eval_point(0.5), curve.eval_point(0.5));

    let dense = curve.refine_knots(&[0.25, 0.5, 0.75]).unwrap();

    let elevated = curve.elevate_degree(1).unwrap();
    assert_eq!(elevated.degree(), 3);

    let (simpler, removed) = dense.remove_knot(0.5, 1, 1e-9).unwrap();
    assert_eq!((removed, simpler.knots()), (1, &[0.0, 0.25, 0.75, 1.0][..]));
    ```

Surfaces offer the same per direction: `insert_u_knot` /
`insert_v_knot`, `refine_u_knots` / `refine_v_knots`,
`elevate_u_degree` / `elevate_v_degree`, `remove_u_knot(s)` /
`remove_v_knot(s)`. Each removal is kept only if a proven bound keeps
every point within `max_deviation`, with successive bounds adding up.

=== "Python"

    ```python
    grown = surf.insert_u_knot(0.5, 1).elevate_v_degree(1)
    ```

=== "Rust"

    ```rust
    let grown = surf.insert_u_knot(0.5, 1)?.elevate_v_degree(1)?;
    ```

## Split, extract, reverse

`split_at` cuts into two clamped pieces, `extract` keeps a sub-range,
`reversed` mirrors the parameter, `to_bezier_pieces` /
`to_bezier_patches` decompose into Bézier spans. Ranges may cross the
seam (up to one period) on periodic directions.

=== "Python"

    ```python
    left, right = curve.split_at(0.5)
    piece = curve.extract(0.0, 0.5)
    rev = curve.reversed()
    assert rev.eval_point(0.0) == curve.eval_point(1.0)
    assert len(curve.to_bezier_pieces()) == 1  # single span
    ```

=== "Rust"

    ```rust
    let (left, right) = curve.split_at(0.5).unwrap();
    let piece = curve.extract(0.0, 0.5).unwrap();
    let rev = curve.reversed();
    assert_eq!(rev.eval_point(0.0), curve.eval_point(1.0));
    assert_eq!(curve.to_bezier_pieces().len(), 1); // single span
    ```

The pieces are `BezierCurve` / `BezierCurve2D` (and, for surfaces,
`BezierSurface`): a Bézier type holds `degree + 1` poles, with no knot
vector, and keeps the span it came from as its bounds. `to_bspline()`
turns one back into a B-spline, and `from_bspline` accepts a B-spline
that is a single clamped span.

The pieces are what the viewer can *not* yet take directly — trimmed
curves have no scene mapping in the Python client
([viewer#6](https://github.com/geomcore/geomcore-viewer/issues/6)) —
so sample them to polylines for display.

## Fitting

Interpolation passes through the points; `approximate` fits with a
fixed pole count; `approximate_within` adds knots where the deviation
is largest until every point is within `max_deviation`.

=== "Python"

    ```python
    points = curve.eval_points([0.0, 0.25, 0.5, 0.75, 1.0])
    interp = BSplineCurve.interpolate(points, 3)
    approx = BSplineCurve.approximate(points, 2, 3)
    within = BSplineCurve.approximate_within(points, 2, 1e-6)
    ```

=== "Rust"

    ```rust
    use geomcore::PointParametrization;

    let points = curve.eval_points(&[0.0, 0.25, 0.5, 0.75, 1.0]);
    let interp = BSplineCurve::interpolate(&points, 3, PointParametrization::Centripetal).unwrap();
    let approx = BSplineCurve::approximate(&points, 2, 3, PointParametrization::Centripetal).unwrap();
    let within = BSplineCurve::approximate_within(&points, 2, 1e-6, PointParametrization::Centripetal).unwrap();
    ```

`BSplineSurface::approximate` fits a grid by separating rows and
columns (one factored system per direction); `approximate_within`
does the same within a deviation.

## Exact conversions

Circles and ellipses convert to exact periodic rational quadratics;
trimmed lines, arcs and conic pieces become clamped B-splines over the
same range; planes, cylinders, cones, spheres and tori to exact
patches.

=== "Python"

    ```python
    circle_bs = Circle(Point.origin(), Vector.unit_z(), 2.0).to_bspline()
    assert circle_bs.is_periodic()

    ellipse_bs = Ellipse(Point.origin(), Vector.unit_z(), Vector.unit_x(), 3.0, 1.0).to_bspline()

    seg_bs = TrimmedCurve.segment(Point(0, 0, 0), Point(1, 0, 0)).to_bspline()

    patch = Cylinder(Point.origin(), Vector.unit_z(), 2.0).to_bspline_patch((0.0, 1.0), (-1.0, 1.0))
    ```

=== "Rust"

    ```rust
    use geomcore::{Circle, Ellipse, Point, Vector};

    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let circle_bs = circle.to_bspline();
    assert!(circle_bs.is_periodic());

    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.0).unwrap();
    let ellipse_bs = ellipse.to_bspline();
    ```

## Measure and bound

Arc length integrates the speed to a requested error; `bounding_box`
gives the tight pole box (elementary shapes get closed-form boxes) —
available on every 3D curve and surface.

=== "Python"

    ```python
    assert round(curve.length(1e-9), 6) == 1.623225
    bb = curve.bounding_box()
    assert (bb.min(), bb.max()) == (Point(0, 0, 0), Point(1, 1, 0))
    ```

=== "Rust"

    ```rust
    use geomcore::HasBoundingBox;

    let len = curve.length(1e-9);
    let bb = curve.bounding_box().unwrap();
    assert_eq!((bb.min(), bb.max()), (Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0)));
    ```

Boxes have no viewer mapping yet
([viewer#8](https://github.com/geomcore/geomcore-viewer/issues/8)).
