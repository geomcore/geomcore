# NURBS curves

Create, evaluate and sample `BSplineCurve`. Editing (knots, degree,
splitting) and fitting live in [Editing & fitting](editing-fitting.md).

## Create

Three components: degree, poles, knots with multiplicities. A clamped,
single-span quadratic is a Bézier curve — the same one geomdl's
[basics](https://nurbs-python.readthedocs.io/en/5.x/basics.html)
start from.

=== "Python"

    ```python
    from geomcore import Point
    from geomcore.curves import BSplineCurve

    poles = [Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)]
    curve = BSplineCurve(2, poles, [0.0, 1.0], [3, 3], False)
    assert curve.degree() == 2
    ```

=== "Rust"

    ```rust
    use geomcore::prelude::*;
    use geomcore::{BSplineCurve, Point};

    let poles = vec![
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    ];
    let curve = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    assert_eq!(curve.degree(), 2);
    ```

Rational (NURBS) curves take one weight per pole via
`BSplineCurve.new_rational` / `new_rational`; periodic ones pass
`periodic=True` with a cyclic knot vector.

Textbooks write the knots flat, each repeated by its multiplicity:
`from_flat_knots` takes that form and `flat_knots()` gives it back.
`clamped_uniform(degree, poles)` spaces the knots evenly over `[0, 1]`,
clamped at both ends. `with_pole(index, point)` and
`with_weight(index, weight)` return the curve with one pole moved or
reweighted, so a pole can be dragged without rebuilding the curve.
A weight on a non-rational curve makes it rational, every other pole
weighing 1.

=== "Python"

    ```python
    flat = BSplineCurve.from_flat_knots(2, poles, [0, 0, 0, 1, 1, 1])
    assert flat == curve and curve.flat_knots() == [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    uniform = BSplineCurve.clamped_uniform(2, poles)    # the same curve again
    arc = curve.with_weight(1, 2 ** -0.5)               # a quarter circle
    dragged = curve.with_pole(-1, Point(0, 2, 0))
    ```

=== "Rust"

    ```rust
    let flat = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let same = BSplineCurve::from_flat_knots(2, curve.poles().to_vec(), flat, None, false).unwrap();
    assert_eq!(same, curve);
    let arc = curve.with_weight(1, std::f64::consts::FRAC_1_SQRT_2).unwrap();
    let dragged = curve.with_pole(2, Point::new(0.0, 2.0, 0.0)).unwrap();
    ```

## Evaluate

=== "Python"

    ```python
    p = curve.point_at(0.5)
    assert (p.x, p.y, p.z) == (0.75, 0.75, 0.0)
    pts = curve.points_at([0.0, 0.25, 0.5, 0.75, 1.0])
    ```

=== "Rust"

    ```rust
    let p = curve.point_at(0.5);
    assert_eq!(p, Point::new(0.75, 0.75, 0.0));
    let pts = curve.points_at(&[0.0, 0.25, 0.5, 0.75, 1.0]);
    ```

Derivatives come from `derivative` (order 1–2); `first_derivatives`/`second_derivatives`
return the point with its derivatives from a single basis pass.
`bounds()` gives the parameter range.

`continuity(u)` reads the representation: a knot of multiplicity `m`
gives `C^(degree - m)`, whatever the poles around it. Where a knot of
multiplicity `degree` lies between poles in line, the curve is still
smooth, and `geometric_continuity(u)` says so. It compares the two
sides' unit tangents (`G1`) and curvature vectors (`G2`), and never
reports less than `continuity`.

## Sampling density

geomdl controls plot density with `delta`/`sample_size`. Here sampling
with guarantees is `discretize`: points joined by segments, each
carrying a proven chord-deviation and angle-deviation bound — see
[Discretization](../guides/discretization.md).

=== "Python"

    ```python
    d = curve.discretize(0.01)
    assert d.points[0] == curve.point_at(0.0)
    ```

=== "Rust"

    ```rust
    use geomcore::DiscretizationCriteria;

    let d = curve.discretize(DiscretizationCriteria::new(0.01).unwrap());
    assert_eq!(d.points[0], curve.point_at(0.0));
    ```

## View it

From a script or notebook, send the real objects — B-splines travel
sampled (curve → 256-segment polyline) until native protocol types land
([viewer#5](https://github.com/geomcore/geomcore-viewer/issues/5)):

=== "Python"

    ```python
    import geomcore_viewer as gv

    view = gv.Viewer("Quadratic B-spline")
    view.add(curve, label="C")
    view.polyline([Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)],
                  label="control polygon", dashed=True)
    ```

=== "Rust"

    ```rust
    use geomcore_viewer::{Geom, Opts, Viewer};

    let mut view = Viewer::new("Quadratic B-spline")?;
    view.add(Geom::from(&curve), Opts::new().label("C"))?;
    ```
