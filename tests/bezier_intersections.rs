//! A Bézier curve or surface takes the same provably complete intersection
//! search as the B-spline it is.

use geomcore::{
    BSplineCurve, BezierCurve, BezierSurface, Point, Tolerance, curve_curve_extrema,
    intersect_curve_curve, intersect_curve_surface,
};

fn arch() -> BezierCurve {
    BezierCurve::new(
        2,
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
        ],
    )
    .unwrap()
}

fn chord() -> BezierCurve {
    BezierCurve::new(
        1,
        vec![Point::new(0.0, 0.5, 0.0), Point::new(2.0, 0.5, 0.0)],
    )
    .unwrap()
}

/// `z = u v` over the unit square.
fn saddle() -> BezierSurface {
    BezierSurface::new(
        1,
        1,
        vec![
            vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
            vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
        ],
    )
    .unwrap()
}

#[test]
fn test_bezier_curves_intersect_completely() {
    let tol = Tolerance::DEFAULT;
    let (arch, chord) = (arch(), chord());
    let as_splines = intersect_curve_curve(&arch.to_bspline(), &chord.to_bspline(), tol);
    let as_beziers = intersect_curve_curve(&arch, &chord, tol);
    assert!(as_splines.is_complete());
    assert!(as_beziers.is_complete());
    assert_eq!(as_beziers.hits().len(), as_splines.hits().len());
    let extrema = curve_curve_extrema(&arch, &chord, tol);
    assert!(extrema.is_complete());
}

#[test]
fn test_bezier_surfaces_intersect_completely() {
    let tol = Tolerance::DEFAULT;
    let vertical = BSplineCurve::new(
        1,
        vec![Point::new(0.5, 0.5, -1.0), Point::new(0.5, 0.5, 1.0)],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
    )
    .unwrap();
    let on_bezier = intersect_curve_surface(&vertical, &saddle(), tol);
    assert!(on_bezier.is_complete());
    assert_eq!(on_bezier.hits().len(), 1);
    let on_bezier_curve = intersect_curve_surface(&chord(), &saddle(), tol);
    assert!(on_bezier_curve.is_complete());
}
