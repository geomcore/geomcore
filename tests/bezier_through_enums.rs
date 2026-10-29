//! A Bézier curve or surface is a variant of the shape enums and an
//! operand of the routed `intersect`, both read as the B-spline it is.

use geomcore::curves::{Curve, Curve2D, TrimmedCurve, TrimmedCurve2D};
use geomcore::surfaces::Surface;
use geomcore::*;

const TOL: Tolerance = Tolerance::DEFAULT;

/// Cubic arch from (0,0,0) to (3,0,0), rising to y = 1.5.
fn arch() -> BezierCurve {
    BezierCurve::new(
        3,
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.0),
            Point::new(2.0, 2.0, 0.0),
            Point::new(3.0, 0.0, 0.0),
        ],
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

/// Quadratic arch from (0,0) to (2,0), rising to y = 1.
fn arch_2d() -> BezierCurve2D {
    BezierCurve2D::new(
        2,
        vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 2.0),
            Point2D::new(2.0, 0.0),
        ],
    )
    .unwrap()
}

#[test]
fn test_curve_enum_reads_the_bezier_curve() {
    let arch = arch();
    let curve = Curve::from(arch.clone());
    assert_eq!(curve, Curve::Bezier(arch.clone()));
    for u in [0.0, 0.3, 0.5, 1.0] {
        assert_eq!(curve.eval_point(u), arch.eval_point(u));
        assert_eq!(curve.eval_d2(u), arch.eval_d2(u));
    }
    assert_eq!(curve.bounds(), arch.bounds());
    assert_eq!(curve.as_elementary(), None);
    assert_eq!(curve.as_bspline(), Some(&arch.to_bspline()));
    assert_eq!(curve.to_bspline(), Ok(arch.to_bspline()));
    assert_eq!(curve.bounding_box(), arch.bounding_box());
    let query = Point::new(1.5, 3.0, 0.0);
    assert_eq!(
        curve.project_point(query, TOL),
        arch.project_point(query, TOL)
    );
}

#[test]
fn test_trimmed_curve_enum_discretizes_the_bezier_curve() {
    let criteria = DiscretizationCriteria::new(1e-3).unwrap();
    let through_enum = TrimmedCurve::new(Curve::from(arch()), 0.2, 0.9).unwrap();
    let direct = TrimmedCurve::new(arch(), 0.2, 0.9).unwrap();
    assert_eq!(
        through_enum.discretize(criteria),
        direct.discretize(criteria)
    );
    assert_eq!(through_enum.bounding_box(), direct.bounding_box());
}

#[test]
fn test_curve_2d_enum_reads_the_bezier_curve() {
    let arch = arch_2d();
    let curve = Curve2D::from(arch.clone());
    for u in [0.0, 0.3, 1.0] {
        assert_eq!(curve.eval_point(u), arch.eval_point(u));
    }
    assert_eq!(curve.as_bspline(), Some(&arch.to_bspline()));
    assert_eq!(curve.to_bspline(), Ok(arch.to_bspline()));
    assert_eq!(curve.bounding_box(), arch.bounding_box());
    let criteria = DiscretizationCriteria::new(1e-3).unwrap();
    assert_eq!(
        TrimmedCurve2D::new(curve, 0.0, 1.0)
            .unwrap()
            .discretize(criteria),
        TrimmedCurve2D::new(arch, 0.0, 1.0)
            .unwrap()
            .discretize(criteria)
    );
}

#[test]
fn test_surface_enum_reads_the_bezier_surface() {
    let saddle = saddle();
    let surface = Surface::from(&saddle);
    assert_eq!(surface, Surface::Bezier(saddle.clone()));
    assert_eq!(surface.eval_point(0.3, 0.4), saddle.eval_point(0.3, 0.4));
    assert_eq!(surface.eval_d1(0.3, 0.4), saddle.eval_d1(0.3, 0.4));
    assert_eq!(surface.as_elementary(), None);
    assert_eq!(surface.as_bspline(), Some(&saddle.to_bspline()));
    assert_eq!(surface.to_bspline(), Ok(saddle.to_bspline()));
    assert_eq!(surface.bounding_box(), saddle.bounding_box());
    assert_eq!(
        surface.u_iso(0.5).unwrap(),
        Curve::Bezier(saddle.u_iso(0.5).unwrap())
    );
    assert_eq!(
        surface.v_iso(0.25).unwrap(),
        Curve::Bezier(saddle.v_iso(0.25).unwrap())
    );
    assert_eq!(surface.reversed_u(), Surface::Bezier(saddle.reversed_u()));
    assert_eq!(surface.reversed_v(), Surface::Bezier(saddle.reversed_v()));
    assert_eq!(
        surface.u_reversed_parameter(0.3),
        saddle.u_reversed_parameter(0.3)
    );
}

#[test]
fn test_bezier_iso_curves_lie_on_the_surface() {
    let saddle = saddle();
    let (u_iso, v_iso) = (saddle.u_iso(0.5).unwrap(), saddle.v_iso(0.25).unwrap());
    for t in [0.0, 0.4, 1.0] {
        assert!(u_iso.eval_point(t).distance(saddle.eval_point(0.5, t)) < 1e-15);
        assert!(v_iso.eval_point(t).distance(saddle.eval_point(t, 0.25)) < 1e-15);
    }
    assert!(saddle.u_iso(2.0).is_err());
}

/// The Bézier curve meets every curve in both orders, with the hits the
/// B-spline it is meets.
#[test]
fn test_bezier_curve_intersects_curves_as_its_bspline() {
    let (arch, spline) = (arch(), arch().to_bspline());
    let line = Line::new(Point::new(0.0, 1.0, 0.0), Vector::X).unwrap();
    let circle = Circle::new(Point::new(1.5, 0.0, 0.0), Vector::Z, 1.0).unwrap();
    let other = BezierCurve::new(
        1,
        vec![Point::new(0.0, 1.0, 0.0), Point::new(3.0, 0.5, 0.0)],
    )
    .unwrap();
    let r = arch.intersect(&line, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert_eq!(r.hits().len(), 2);
    assert_eq!(r, spline.intersect(&line, TOL));
    assert_eq!(line.intersect(&arch, TOL), line.intersect(&spline, TOL));
    assert_eq!(arch.intersect(&circle, TOL), spline.intersect(&circle, TOL));
    assert_eq!(circle.intersect(&arch, TOL), circle.intersect(&spline, TOL));
    assert_eq!(arch.intersect(&other, TOL).hits().len(), 2);
    assert_eq!(
        arch.intersect(&other.to_bspline(), TOL),
        spline.intersect(&other.to_bspline(), TOL)
    );
    let through_enum = Curve::from(arch.clone());
    assert_eq!(
        through_enum.intersect(&Curve::from(line), TOL),
        spline.intersect(&line, TOL)
    );
    assert_eq!(arch.intersect(&through_enum, TOL).overlaps().len(), 1);
}

#[test]
fn test_trimmed_bezier_curve_keeps_only_its_hits() {
    let line = Line::new(Point::new(0.0, 1.0, 0.0), Vector::X).unwrap();
    let half = TrimmedCurve::new(arch(), 0.0, 0.5).unwrap();
    assert_eq!(half.intersect(&line, TOL).hits().len(), 1);
    assert_eq!(line.intersect(&half, TOL).hits().len(), 1);
    let plane = Plane::new(Point::new(0.0, 1.0, 0.0), Vector::Y).unwrap();
    assert_eq!(half.intersect(&plane, TOL).hits().len(), 1);
}

#[test]
fn test_bezier_curve_intersects_surfaces_as_its_bspline() {
    let (arch, spline) = (arch(), arch().to_bspline());
    let plane = Plane::new(Point::new(0.0, 1.0, 0.0), Vector::Y).unwrap();
    let sphere = Sphere::new(Point::new(1.5, 0.0, 0.0), 1.0).unwrap();
    let r = arch.intersect(&plane, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert_eq!(r.hits().len(), 2);
    assert_eq!(r, spline.intersect(&plane, TOL));
    assert_eq!(plane.intersect(&arch, TOL), plane.intersect(&spline, TOL));
    assert_eq!(arch.intersect(&sphere, TOL), spline.intersect(&sphere, TOL));
    assert_eq!(
        Curve::from(arch.clone()).intersect(&Surface::from(sphere), TOL),
        spline.intersect(&sphere, TOL)
    );
    let vertical = BezierCurve::new(
        1,
        vec![Point::new(0.5, 0.5, -1.0), Point::new(0.5, 0.5, 1.0)],
    )
    .unwrap();
    let on_saddle = vertical.intersect(&saddle(), TOL);
    assert!(on_saddle.is_complete());
    assert_eq!(on_saddle.hits().len(), 1);
    assert_eq!(saddle().intersect(&vertical, TOL).hits().len(), 1);
}

#[test]
fn test_bezier_surface_intersects_surfaces_as_its_bspline() {
    let (saddle, spline) = (saddle(), saddle().to_bspline());
    let plane = Plane::new(Point::new(0.0, 0.0, 0.25), Vector::Z).unwrap();
    let r = saddle.intersect(&plane, TOL);
    assert!(!r.curves().is_empty());
    assert_eq!(r, spline.intersect(&plane, TOL));
    assert_eq!(plane.intersect(&saddle, TOL), plane.intersect(&spline, TOL));
    assert_eq!(
        Surface::from(&saddle).intersect(&Surface::from(plane), TOL),
        spline.intersect(&plane, TOL)
    );
    let tilted = BezierSurface::new(
        1,
        1,
        vec![
            vec![Point::new(0.0, 0.0, 0.5), Point::new(0.0, 1.0, 0.5)],
            vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0)],
        ],
    )
    .unwrap();
    assert_eq!(
        saddle.intersect(&tilted, TOL),
        spline.intersect(&tilted.to_bspline(), TOL)
    );
}

#[test]
fn test_bezier_curve_2d_intersects_curves_as_its_bspline() {
    let (arch, spline) = (arch_2d(), arch_2d().to_bspline());
    let line = Line2D::new(Point2D::new(0.0, 0.5), Vector2D::new(1.0, 0.0)).unwrap();
    let r = arch.intersect(&line, TOL);
    assert_eq!(r.hits().len(), 2);
    assert_eq!(r, spline.intersect(&line, TOL));
    assert_eq!(line.intersect(&arch, TOL), line.intersect(&spline, TOL));
    assert_eq!(
        Curve2D::from(arch.clone()).intersect(&Curve2D::from(line), TOL),
        spline.intersect(&line, TOL)
    );
    let half = TrimmedCurve2D::new(arch, 0.0, 0.5).unwrap();
    assert_eq!(half.intersect(&line, TOL).hits().len(), 1);
}
