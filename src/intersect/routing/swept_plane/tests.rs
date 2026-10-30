use super::*;
use crate::curves::{BSplineCurve, Ellipse, Parabola};
use crate::surfaces::{ParametricSurface, Surface};
use crate::{Axis, Vector};

const TOL: Tolerance = Tolerance::DEFAULT;

fn plane(point: (f64, f64, f64), normal: Vector) -> Plane {
    Plane::new(Point::new(point.0, point.1, point.2), normal).unwrap()
}

fn circle(center: (f64, f64, f64), normal: Vector, radius: f64) -> Curve {
    Circle::new(Point::new(center.0, center.1, center.2), normal, radius)
        .unwrap()
        .into()
}

/// A torus of major radius 3 and minor radius 1 about z.
fn torus() -> RevolutionSurface {
    let profile = circle((3.0, 0.0, 0.0), Vector::Y, 1.0);
    RevolutionSurface::new(profile, Axis::new(Point::ORIGIN, Vector::Z).unwrap()).unwrap()
}

/// A sphere of radius 2 about the origin, turned about z.
fn sphere() -> RevolutionSurface {
    let profile = circle((0.0, 0.0, 0.0), Vector::Y, 2.0);
    RevolutionSurface::new(profile, Axis::new(Point::ORIGIN, Vector::Z).unwrap()).unwrap()
}

/// Checks that samples of every curve lie on both surfaces.
fn assert_on_both(
    result: &SurfaceSurfaceIntersection,
    first: &impl ParametricSurface,
    second: &impl ParametricSurface,
) {
    for cut in result.curves() {
        let (start, end) = cut.curve.bounds();
        let (start, end) = (start.max(-5.0), end.min(5.0));
        for index in 0..=16 {
            let point = cut
                .curve
                .eval_point(start + (end - start) * index as f64 / 16.0);
            assert!(
                first.project_point(point, TOL).distance < 1e-9,
                "{point:?} off the first surface"
            );
            assert!(
                second.project_point(point, TOL).distance < 1e-9,
                "{point:?} off the second surface"
            );
        }
    }
}

#[test]
fn test_plane_holding_the_direction_cuts_generator_lines() {
    let surface =
        ExtrusionSurface::new(circle((0.0, 0.0, 0.0), Vector::Z, 2.0), Vector::Z).unwrap();
    let cutter = plane((1.0, 0.0, 0.0), Vector::X);
    let result = surface.intersect(&cutter, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert_eq!(result.curves().len(), 2);
    assert!(
        result
            .curves()
            .iter()
            .all(|cut| !cut.tangent && matches!(cut.curve, Curve::Line(_)))
    );
    assert_on_both(&result, &surface, &cutter);
}

#[test]
fn test_plane_touching_an_extrusion_cuts_one_tangent_line() {
    let surface =
        ExtrusionSurface::new(circle((0.0, 0.0, 0.0), Vector::Z, 2.0), Vector::Z).unwrap();
    let result = surface.intersect(&plane((2.0, 0.0, 0.0), Vector::X), TOL);
    assert_eq!(result.curves().len(), 1);
    assert!(result.curves()[0].tangent);
    let missing = surface.intersect(&plane((3.0, 0.0, 0.0), Vector::X), TOL);
    assert_eq!(missing.method(), IntersectionMethod::Exact);
    assert!(missing.is_empty());
}

#[test]
fn test_oblique_plane_cuts_the_carried_profile() {
    let profile: Curve = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.0)
        .unwrap()
        .into();
    let surface = ExtrusionSurface::new(profile, Vector::new(0.3, -0.2, 1.0)).unwrap();
    let cutter = plane((0.0, 0.0, 2.0), Vector::new(0.1, 0.4, 1.0));
    let result = surface.intersect(&cutter, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert_eq!(result.curves().len(), 1);
    assert!(matches!(result.curves()[0].curve, Curve::Ellipse(_)));
    assert_on_both(&result, &surface, &cutter);
}

#[test]
fn test_spline_profile_crossings_are_numeric_but_complete() {
    let poles = vec![
        Point::new(-2.0, 0.0, 0.0),
        Point::new(-1.0, 2.0, 0.0),
        Point::new(1.0, -2.0, 0.0),
        Point::new(2.0, 1.0, 0.0),
    ];
    let profile: Curve = BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false)
        .unwrap()
        .into();
    let surface = ExtrusionSurface::new(profile, Vector::Z).unwrap();
    let cutter = plane((0.0, 0.0, 0.0), Vector::Y);
    let result = surface.intersect(&cutter, TOL);
    assert_eq!(result.method(), IntersectionMethod::Numeric);
    assert!(result.is_complete());
    assert_eq!(result.curves().len(), 3);
    assert_on_both(&result, &surface, &cutter);
    let carried = surface.intersect(&plane((0.0, 0.0, 4.0), Vector::Z), TOL);
    assert_eq!(carried.method(), IntersectionMethod::Exact);
    assert!(matches!(carried.curves()[0].curve, Curve::BSpline(_)));
}

#[test]
fn test_planar_extrusion_of_a_conic_cuts_a_segment() {
    // A circle in the xz plane swept along z covers the strip |x| <= 2.
    let surface =
        ExtrusionSurface::new(circle((0.0, 0.0, 0.0), Vector::Y, 2.0), Vector::Z).unwrap();
    let cutter = plane((0.0, 0.0, 5.0), Vector::Z);
    let result = surface.intersect(&cutter, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    let segment = &result.curves()[0].curve;
    let ends = [segment.eval_point(0.0), segment.eval_point(1.0)];
    assert!(ends[0].distance(Point::new(-2.0, 0.0, 5.0)) < 1e-12);
    assert!(ends[1].distance(Point::new(2.0, 0.0, 5.0)) < 1e-12);
    // A parabola's strip is a half-plane: no whole curve, so it marches.
    let parabola: Curve = Parabola::new(Point::ORIGIN, Vector::Y, Vector::X, 1.0)
        .unwrap()
        .into();
    let open = ExtrusionSurface::new(parabola, Vector::Z).unwrap();
    assert_eq!(
        open.intersect(&cutter, TOL).method(),
        IntersectionMethod::Numeric
    );
}

#[test]
fn test_profile_in_the_plane_marches() {
    let profile: Curve = Line::new(Point::ORIGIN, Vector::Y).unwrap().into();
    let surface = ExtrusionSurface::new(profile, Vector::Z).unwrap();
    let result = surface.intersect(&plane((0.0, 0.0, 0.0), Vector::X), TOL);
    assert_eq!(result.method(), IntersectionMethod::Numeric);
}

#[test]
fn test_perpendicular_plane_cuts_circles_of_a_torus() {
    let surface = torus();
    let cutter = plane((0.0, 0.0, 0.5), Vector::Z);
    let result = surface.intersect(&cutter, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    let mut radii: Vec<f64> = result
        .curves()
        .iter()
        .map(|cut| match &cut.curve {
            Curve::Circle(circle) => circle.radius(),
            other => panic!("expected a circle, got {other:?}"),
        })
        .collect();
    radii.sort_by(f64::total_cmp);
    let half_width = 0.75_f64.sqrt();
    assert!((radii[0] - (3.0 - half_width)).abs() < 1e-12);
    assert!((radii[1] - (3.0 + half_width)).abs() < 1e-12);
    assert_on_both(&result, &surface, &cutter);

    let top = surface.intersect(&plane((0.0, 0.0, 1.0), Vector::Z), TOL);
    assert_eq!(top.curves().len(), 1);
    assert!(top.curves()[0].tangent);
}

#[test]
fn test_perpendicular_plane_cuts_a_sphere_once() {
    let surface = sphere();
    let result = surface.intersect(&plane((0.0, 0.0, 1.0), Vector::Z), TOL);
    assert_eq!(result.curves().len(), 1);
    assert!(!result.curves()[0].tangent);
    let pole = surface.intersect(&plane((0.0, 0.0, 2.0), Vector::Z), TOL);
    assert!(pole.curves().is_empty());
    assert_eq!(pole.contact_points().len(), 1);
    assert!(pole.contact_points()[0].distance(Point::new(0.0, 0.0, 2.0)) < 1e-12);
}

#[test]
fn test_crossing_on_the_axis_is_a_contact_point() {
    let profile: Curve = Line::new(Point::ORIGIN, Vector::new(1.0, 0.0, 1.0))
        .unwrap()
        .into();
    let cone =
        RevolutionSurface::new(profile, Axis::new(Point::ORIGIN, Vector::Z).unwrap()).unwrap();
    let result = cone.intersect(&plane((0.0, 0.0, 0.0), Vector::Z), TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.curves().is_empty());
    assert_eq!(result.contact_points(), &[Point::ORIGIN]);
}

#[test]
fn test_plane_through_the_axis_cuts_both_halves() {
    let surface = torus();
    let cutter = plane((0.0, 0.0, 0.0), Vector::new(1.0, 1.0, 0.0));
    let result = surface.intersect(&cutter, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert_eq!(result.curves().len(), 2);
    assert_on_both(&result, &surface, &cutter);

    let symmetric = sphere().intersect(&cutter, TOL);
    assert_eq!(symmetric.curves().len(), 1);
    assert_on_both(&symmetric, &sphere(), &cutter);
}

#[test]
fn test_other_positions_march() {
    let oblique = torus().intersect(&plane((0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 1.0)), TOL);
    assert_eq!(oblique.method(), IntersectionMethod::Numeric);
    let off_axis = torus().intersect(&plane((0.0, 0.5, 0.0), Vector::Y), TOL);
    assert_eq!(off_axis.method(), IntersectionMethod::Numeric);
}

#[test]
fn test_either_order_and_the_enums_agree() {
    let surface = torus();
    let cutter = plane((0.0, 0.0, 0.5), Vector::Z);
    let expected = surface.intersect(&cutter, TOL);
    assert_eq!(cutter.intersect(&surface, TOL), expected);
    let (wrapped, wrapped_plane) = (Surface::from(surface.clone()), Surface::from(cutter));
    assert_eq!(wrapped.intersect(&wrapped_plane, TOL), expected);
    assert_eq!(wrapped_plane.intersect(&wrapped, TOL), expected);
    let elementary = crate::surfaces::ElementarySurface::Plane(cutter);
    assert_eq!(elementary.intersect(&surface, TOL), expected);
}
