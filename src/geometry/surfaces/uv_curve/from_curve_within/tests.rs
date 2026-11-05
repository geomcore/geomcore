use std::f64::consts::{FRAC_PI_2, TAU};

use super::{Curve3DToUvApproximation, Curve3DToUvApproximationError, curve_3d_to_uv_within};
use crate::geometry::curves::{
    BSplineCurve, Circle, Ellipse, Line, Line2D, ParametricCurve, ParametricCurve2D, TrimmedCurve,
    TrimmedCurve2D,
};
use crate::geometry::surfaces::{BSplineSurface, Cylinder, ParametricSurface, Plane, Sphere};
use crate::{Point, Point2D, Tolerance, Vector, Vector2D, uv_curve_to_3d_within};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The largest same-parameter distance between the surface's image of the
/// UV curve and the 3D curve, at 2 000 parameters.
fn measured_deviation(
    surface: &impl ParametricSurface,
    curve: &impl ParametricCurve,
    approximation: &Curve3DToUvApproximation,
) -> f64 {
    let (first, last) = curve.bounds();
    let uv_curve = approximation.curve();
    (0..=2000)
        .map(|k| {
            let t = first + (last - first) * (k as f64 / 2000.0);
            let uv = uv_curve.point_at(t);
            surface.point_at(uv.x, uv.y).distance(curve.point_at(t))
        })
        .fold(0.0, f64::max)
}

/// Checks the result's shape and that it is within its deviation, and
/// returns it.
fn approximated(
    surface: &impl ParametricSurface,
    curve: &impl ParametricCurve,
    max_deviation: f64,
) -> Curve3DToUvApproximation {
    let approximation = curve_3d_to_uv_within(surface, curve, max_deviation, TOL).unwrap();
    let uv_curve = approximation.curve();
    assert_eq!(uv_curve.degree(), 5);
    assert_eq!(uv_curve.bounds(), curve.bounds());
    assert!(!uv_curve.is_periodic());
    assert!(approximation.deviation() <= max_deviation);
    assert!(measured_deviation(surface, curve, &approximation) <= max_deviation);
    approximation
}

#[test]
fn an_ellipse_across_a_cylinder_seam_has_one_continuous_uv_curve() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let ellipse = Ellipse::new(
        Point::ORIGIN,
        Vector::new(0.0, 1.0, 1.0),
        Vector::new(0.0, 1.0, -1.0),
        2.0 * 2f64.sqrt(),
        2.0,
    )
    .unwrap();
    let approximation = approximated(&cylinder, &ellipse, 1e-7);
    let uv_curve = approximation.curve();
    // `u` rises steadily from π/2 by one turn: no jump at the seam.
    let mut previous = uv_curve.point_at(0.0).x;
    assert!((previous - FRAC_PI_2).abs() < 1e-9);
    for k in 1..=400 {
        let u = uv_curve.point_at(TAU * (k as f64 / 400.0)).x;
        assert!(u > previous && u - previous < 0.1, "u jumps at step {k}");
        previous = u;
    }
    assert!((previous - FRAC_PI_2 - TAU).abs() < 1e-9);
}

#[test]
fn a_helix_maps_back_to_its_slanted_line() {
    // Three turns of a helix of pitch 4 on a radius-5 cylinder, as a
    // B-spline within 1e-10, then back to UV: the line it came from.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 5.0).unwrap();
    let slope = Vector2D::new(TAU, 4.0);
    let line = Line2D::new(Point2D::new(1.0, 0.0), slope).unwrap();
    let turns = TrimmedCurve2D::new(line, 0.0, 3.0 * slope.magnitude()).unwrap();
    let helix = uv_curve_to_3d_within(&cylinder, &turns, 1e-10, TOL)
        .unwrap()
        .into_curve();
    let approximation = approximated(&cylinder, &helix, 1e-7);
    let uv_curve = approximation.curve();
    for k in 0..=300 {
        let t = turns.bounds().1 * (k as f64 / 300.0);
        let (found, expected) = (uv_curve.point_at(t), turns.point_at(t));
        assert!(
            found.distance(expected) < 1e-7,
            "at {t}: {found:?} against {expected:?}"
        );
    }
}

#[test]
fn a_curve_on_a_spline_surface_maps_back_to_its_uv_curve() {
    // A curved patch, and the 3D curve of its diagonal.
    let poles = (0..4)
        .map(|i| {
            (0..4)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    Point::new(x, y, (x * 0.7).sin() + (y * 0.9).cos())
                })
                .collect()
        })
        .collect();
    let surface = BSplineSurface::new(
        3,
        3,
        poles,
        vec![0.0, 1.0],
        vec![4, 4],
        vec![0.0, 1.0],
        vec![4, 4],
        false,
        false,
    )
    .unwrap();
    let diagonal = Line2D::new(Point2D::new(0.1, 0.2), Vector2D::new(0.7, 0.6)).unwrap();
    let diagonal = TrimmedCurve2D::new(diagonal, 0.0, 1.0).unwrap();
    let curve = uv_curve_to_3d_within(&surface, &diagonal, 1e-10, TOL)
        .unwrap()
        .into_curve();
    let approximation = approximated(&surface, &curve, 1e-6);
    let uv_curve = approximation.curve();
    for k in 0..=100 {
        let t = k as f64 / 100.0;
        assert!(uv_curve.point_at(t).distance(diagonal.point_at(t)) < 1e-5);
    }
}

#[test]
fn a_polyline_keeps_its_corners() {
    // A degree-1 B-spline on a plane: C0 at its knot, kept by a knot of
    // multiplicity 5, and linear on each span, which degree 5 holds.
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let polyline = BSplineCurve::new(
        1,
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
        ],
        vec![0.0, 1.0, 2.0],
        vec![2, 1, 2],
        false,
    )
    .unwrap();
    let approximation = approximated(&plane, &polyline, 1e-9);
    let uv_curve = approximation.curve();
    let corner = uv_curve
        .knots()
        .iter()
        .position(|&knot| knot == 1.0)
        .unwrap();
    assert_eq!(uv_curve.multiplicities()[corner], 5);
    assert!(approximation.deviation() < 1e-12);
}

#[test]
fn the_degenerate_pairs_are_refused() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let coaxial = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    for deviation in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            curve_3d_to_uv_within(&cylinder, &coaxial, deviation, TOL),
            Err(Curve3DToUvApproximationError::InvalidDeviation)
        );
    }
    let axis = Line::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(
        curve_3d_to_uv_within(&cylinder, &axis, 1e-6, TOL),
        Err(Curve3DToUvApproximationError::Unbounded)
    );
    let wider = Circle::new(Point::ORIGIN, Vector::Z, 2.5).unwrap();
    assert_eq!(
        curve_3d_to_uv_within(&cylinder, &wider, 1e-6, TOL),
        Err(Curve3DToUvApproximationError::CurveNotOnSurface { distance: 0.5 })
    );
    // A quarter of a meridian of a sphere, ending at its south pole.
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let meridian = Circle::new(Point::ORIGIN, Vector::Y, 1.0).unwrap();
    let to_pole = TrimmedCurve::new(meridian, 0.0, FRAC_PI_2).unwrap();
    assert!(
        meridian
            .point_at(FRAC_PI_2)
            .distance(Point::new(0.0, 0.0, -1.0))
            < 1e-15
    );
    assert_eq!(
        curve_3d_to_uv_within(&sphere, &to_pole, 1e-6, TOL),
        Err(Curve3DToUvApproximationError::SingularParametrization)
    );
    // Short of the pole, it maps.
    let short = TrimmedCurve::new(meridian, 0.0, FRAC_PI_2 - 0.1).unwrap();
    approximated(&sphere, &short, 1e-7);
}
