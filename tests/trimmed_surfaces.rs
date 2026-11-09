//! Intersections of trimmed surfaces: the intersection of the bases,
//! clipped to the parameter ranges of each trimmed side.

use geomcore::*;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_6, PI, TAU};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The points of `piece` at `count` parameters across its range.
fn samples(piece: &IntersectionCurve, count: usize) -> Vec<Point> {
    let (start, end) = piece.parameter_range;
    (0..count)
        .map(|i| {
            piece
                .curve
                .point_at(start + (end - start) * i as f64 / (count - 1) as f64)
        })
        .collect()
}

/// Every curve of `result` lies within `1e-6` of both patches, in their
/// ranges, across its whole parameter range.
fn assert_inside<A: ParametricSurface, B: ParametricSurface>(
    result: &SurfaceSurfaceIntersection,
    a: &A,
    b: &B,
) {
    for piece in result.curves() {
        let (start, end) = piece.parameter_range;
        assert!(start.is_finite() && end.is_finite() && start < end);
        for point in samples(piece, 41) {
            assert!(
                a.project_point(point, TOL).distance <= 1e-6,
                "{point:?} off the first patch"
            );
            assert!(
                b.project_point(point, TOL).distance <= 1e-6,
                "{point:?} off the second patch"
            );
        }
    }
}

/// Whether `point` lies on the boundary of the ranges of `patch`, a
/// periodic parameter compared modulo its period.
fn on_boundary<S: ParametricSurface>(patch: &TrimmedSurface<S>, point: Point) -> bool {
    let (u, v) = patch.basis().parameters_of(point, TOL);
    let near = |x: f64, end: f64, period: Option<f64>| match period {
        Some(period) => {
            let gap = (x - end).rem_euclid(period);
            gap.min(period - gap) <= 1e-7
        }
        None => (x - end).abs() <= 1e-7,
    };
    let basis = patch.basis();
    let (u0, u1) = patch.u_range();
    let (v0, v1) = patch.v_range();
    near(u, u0, basis.u_period())
        || near(u, u1, basis.u_period())
        || near(v, v0, basis.v_period())
        || near(v, v1, basis.v_period())
}

/// The ends of each open curve of `result` lie on the boundary of one of
/// the patches.
fn assert_ends_on_boundaries<A: ParametricSurface, B: ParametricSurface>(
    result: &SurfaceSurfaceIntersection,
    a: &TrimmedSurface<A>,
    b: &TrimmedSurface<B>,
) {
    for piece in result.curves() {
        let (start, end) = piece.parameter_range;
        let (first, last) = (piece.curve.point_at(start), piece.curve.point_at(end));
        if first.distance(last) <= 1e-9 {
            continue;
        }
        for point in [first, last] {
            assert!(
                on_boundary(a, point) || on_boundary(b, point),
                "{point:?} ends inside both patches"
            );
        }
    }
}

/// Every point of the untrimmed intersection inside both patches, away
/// from their boundaries, lies on a clipped curve: the clipping loses
/// nothing.
fn assert_covers<A: ParametricSurface, B: ParametricSurface>(
    clipped: &SurfaceSurfaceIntersection,
    untrimmed: &SurfaceSurfaceIntersection,
    a: &TrimmedSurface<A>,
    b: &TrimmedSurface<B>,
) {
    let polylines: Vec<Vec<Point>> = clipped
        .curves()
        .iter()
        .map(|piece| samples(piece, 2001))
        .collect();
    for piece in untrimmed.curves() {
        let (start, end) = piece.parameter_range;
        let (start, end) = (start.max(-50.0), end.min(50.0));
        for i in 0..400 {
            let point = piece
                .curve
                .point_at(start + (end - start) * i as f64 / 399.0);
            let inside = a.project_point(point, TOL).distance <= 1e-9
                && b.project_point(point, TOL).distance <= 1e-9
                && !on_boundary_within(a, point, 1e-3)
                && !on_boundary_within(b, point, 1e-3);
            if !inside {
                continue;
            }
            let nearest = polylines
                .iter()
                .flatten()
                .map(|sample| sample.distance(point))
                .fold(f64::INFINITY, f64::min);
            assert!(
                nearest <= 1e-2,
                "{point:?} of the untrimmed curve is missing ({nearest})"
            );
        }
    }
}

/// Whether `point` lies within `margin` in parameters of the boundary.
fn on_boundary_within<S: ParametricSurface>(
    patch: &TrimmedSurface<S>,
    point: Point,
    margin: f64,
) -> bool {
    let (u, v) = patch.basis().parameters_of(point, TOL);
    let near = |x: f64, end: f64, period: Option<f64>| {
        let gap = match period {
            Some(period) => {
                let gap = (x - end).rem_euclid(period);
                gap.min(period - gap)
            }
            None => (x - end).abs(),
        };
        gap <= margin
    };
    let basis = patch.basis();
    let ((u0, u1), (v0, v1)) = (patch.u_range(), patch.v_range());
    near(u, u0, basis.u_period())
        || near(u, u1, basis.u_period())
        || near(v, v0, basis.v_period())
        || near(v, v1, basis.v_period())
}

#[test]
fn test_two_trimmed_cylinders_meet_along_bounded_quartic_arcs() {
    let wide = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let narrow = Cylinder::new(Point::ORIGIN, Vector::X, 0.6).unwrap();
    let a = wide.trimmed((0.0, PI), (-3.0, 3.0)).unwrap();
    let b = narrow.trimmed((0.0, TAU), (-3.0, 0.5)).unwrap();
    let result = a.intersect(&b, TOL);
    assert!(result.is_complete());
    assert!(!result.curves().is_empty());
    assert_inside(&result, &a, &b);
    assert_ends_on_boundaries(&result, &a, &b);
    let untrimmed = wide.intersect(&narrow, TOL);
    assert!(result.curves().len() >= untrimmed.curves().len());
    assert_covers(&result, &untrimmed, &a, &b);
}

#[test]
fn test_a_plane_across_a_trimmed_cylinder_gives_an_ellipse_arc() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let half = cylinder.trimmed((0.0, PI), (-5.0, 5.0)).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0)).unwrap();
    let result = plane.intersect(&half, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert_eq!(result.curves().len(), 1);
    let piece = &result.curves()[0];
    assert!(matches!(piece.curve, Curve::Ellipse(_)));
    let (start, end) = piece.parameter_range;
    assert!(
        (end - start - PI).abs() <= 1e-9,
        "half the ellipse: {:?}",
        piece.parameter_range
    );
    assert_inside(&result, &plane, &half);
    for point in [piece.curve.point_at(start), piece.curve.point_at(end)] {
        assert!(on_boundary(&half, point));
    }
    // Either order gives the same pieces.
    let swapped = half.intersect(&plane, TOL);
    assert_eq!(swapped.curves().len(), 1);
    assert_eq!(swapped.curves()[0].curve, piece.curve);
    assert_eq!(swapped.curves()[0].parameter_range, piece.parameter_range);
}

#[test]
fn test_a_range_across_the_seam_keeps_one_piece() {
    // A band of the cylinder across its seam: the circle it cuts is one
    // arc running past the circle's period, not two pieces.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let band = cylinder.trimmed((1.5 * PI, 2.5 * PI), (-1.0, 1.0)).unwrap();
    let plane = Plane::new(Point::new(0.0, 0.0, 0.25), Vector::Z).unwrap();
    let result = band.intersect(&plane, TOL);
    assert_eq!(result.curves().len(), 1);
    let piece = &result.curves()[0];
    let (start, end) = piece.parameter_range;
    assert!(
        (end - start - PI).abs() <= 1e-9,
        "{:?}",
        piece.parameter_range
    );
    assert_inside(&result, &band, &plane);
    for point in [piece.curve.point_at(start), piece.curve.point_at(end)] {
        assert!(on_boundary(&band, point), "{point:?}");
    }
    // The middle of the arc is the seam's point, `u = 0` (= 2π).
    let middle = piece.curve.point_at(0.5 * (start + end));
    let seam = cylinder.point_at(0.0, 0.25);
    assert!(middle.distance(seam) <= 1e-9);
}

#[test]
fn test_a_line_through_a_trimmed_cylinder_becomes_a_segment() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let patch = cylinder.trimmed((0.0, PI), (0.0, 2.0)).unwrap();
    // The plane through the axis cuts two generators; one is in the patch.
    let plane = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    let result = plane.intersect(&patch, TOL);
    assert_eq!(result.curves().len(), 1);
    let piece = &result.curves()[0];
    let (start, end) = piece.parameter_range;
    let length = piece
        .curve
        .point_at(start)
        .distance(piece.curve.point_at(end));
    assert!((length - 2.0).abs() <= 1e-9, "{length}");
    assert_inside(&result, &plane, &patch);
}

#[test]
fn test_a_trimmed_cone_and_generator_cylinder_march_to_the_ranges() {
    // The intersection of the untrimmed pair runs to infinity, so the
    // untrimmed march is cut and incomplete; trimmed, no branch is cut (a
    // march samples, so its result is `Sampled`) and the curves end on the
    // ranges.
    let cone = Cone::new(Point::ORIGIN, Vector::Z, 0.3, 1.0).unwrap();
    let generator = Vector::new(0.3f64.sin(), 0.0, 0.3f64.cos());
    let cylinder = Cylinder::new(Point::new(1.0, 0.0, 0.0), generator, 0.5).unwrap();
    assert!(!cone.intersect(&cylinder, TOL).is_complete());
    let cone_patch = cone.trimmed((0.0, TAU), (0.5, 6.0)).unwrap();
    let cylinder_patch = cylinder.trimmed((0.0, TAU), (-2.0, 8.0)).unwrap();
    let result = cone_patch.intersect(&cylinder_patch, TOL);
    assert_eq!(result.method(), IntersectionMethod::Numeric);
    assert_eq!(*result.completeness(), Completeness::Sampled);
    assert!(!result.curves().is_empty());
    assert_inside(&result, &cone_patch, &cylinder_patch);
    assert_ends_on_boundaries(&result, &cone_patch, &cylinder_patch);
}

#[test]
fn test_two_trimmed_cones_whose_curve_runs_to_infinity_are_clipped_to_the_ranges() {
    // A wide cone across the upright one has generators parallel to its
    // own: the untrimmed quartic is exact, complete, and runs to infinity
    // in branches. Trimmed, the branches are clipped to the ranges, with
    // nothing inside both patches lost.
    let upright = Cone::new(
        Point::new(0.0, 0.0, 1.0 / FRAC_PI_6.tan()),
        Vector::Z,
        FRAC_PI_6,
        1.0,
    )
    .unwrap();
    let apex = Point::new(-3.0, 0.0, 1.0);
    let wide = Cone::new(apex + Vector::X, Vector::X, 1.2, 1.2f64.tan()).unwrap();
    let untrimmed = upright.intersect(&wide, TOL);
    assert_eq!(untrimmed.method(), IntersectionMethod::Exact);
    assert!(untrimmed.is_complete());
    assert!(
        untrimmed
            .curves()
            .iter()
            .any(|piece| matches!(piece.curve, Curve::QuarticBranch(_)))
    );
    let upright_patch = upright.trimmed((0.0, TAU), (-1.0, 8.0)).unwrap();
    let wide_patch = wide.trimmed((0.0, TAU), (-2.0, 6.0)).unwrap();
    let result = upright_patch.intersect(&wide_patch, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    assert!(!result.curves().is_empty());
    assert_inside(&result, &upright_patch, &wide_patch);
    assert_ends_on_boundaries(&result, &upright_patch, &wide_patch);
    assert_covers(&result, &untrimmed, &upright_patch, &wide_patch);
}

#[test]
fn test_contact_points_outside_the_ranges_are_dropped() {
    let sphere = Sphere::new(Point::new(0.0, 0.0, 1.0), 1.0).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(plane.intersect(&sphere, TOL).contact_points().len(), 1);
    let around = plane.trimmed((-1.0, 1.0), (-1.0, 1.0)).unwrap();
    assert_eq!(around.intersect(&sphere, TOL).contact_points().len(), 1);
    let aside = plane.trimmed((1.0, 2.0), (1.0, 2.0)).unwrap();
    let result = aside.intersect(&sphere, TOL);
    assert!(result.is_empty());
}

#[test]
fn test_patches_of_one_surface_are_coincident_where_they_overlap() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let a = plane.trimmed((0.0, 2.0), (0.0, 2.0)).unwrap();
    let overlapping = plane.trimmed((1.0, 3.0), (1.0, 3.0)).unwrap();
    assert!(a.intersect(&overlapping, TOL).is_coincident());
    let apart = plane.trimmed((5.0, 6.0), (5.0, 6.0)).unwrap();
    let result = a.intersect(&apart, TOL);
    assert!(!result.is_coincident());
    assert!(result.is_empty());
    // A patch lies on its untrimmed basis.
    assert!(a.intersect(&plane, TOL).is_coincident());
}

#[test]
fn test_patches_in_disjoint_boxes_do_not_meet() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 1.5).unwrap();
    let low = cylinder.trimmed((0.0, TAU), (-10.0, -5.0)).unwrap();
    let cap = sphere.trimmed((0.0, TAU), (1.0, FRAC_PI_2)).unwrap();
    let result = low.intersect(&cap, TOL);
    assert!(result.is_empty());
    assert!(result.is_complete());
}

#[test]
fn test_an_untrimmed_surface_keeps_its_result() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0)).unwrap();
    let result = plane.intersect(&cylinder, TOL);
    assert_eq!(result.curves().len(), 1);
    assert_eq!(
        result.curves()[0].parameter_range,
        result.curves()[0].curve.bounds()
    );
}
