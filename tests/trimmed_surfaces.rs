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
    // Only the loop at x < 0 reaches the narrow patch, and half of it,
    // y > 0, lies on the wide one: one arc, across the loop's start.
    assert_eq!(result.curves().len(), 1);
    let untrimmed = wide.intersect(&narrow, TOL);
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
fn test_a_trimmed_cone_and_generator_cylinder_cut_the_exact_branches_at_the_ranges() {
    // The intersection of the untrimmed pair runs to infinity, in exact
    // branches; trimmed, they are cut where they leave the ranges, and the
    // curves end there.
    let cone = Cone::new(Point::ORIGIN, Vector::Z, 0.3, 1.0).unwrap();
    let generator = Vector::new(0.3f64.sin(), 0.0, 0.3f64.cos());
    let cylinder = Cylinder::new(Point::new(1.0, 0.0, 0.0), generator, 0.5).unwrap();
    let whole = cone.intersect(&cylinder, TOL);
    assert!(whole.is_complete());
    assert!(
        whole
            .curves()
            .iter()
            .all(|piece| matches!(piece.curve, Curve::QuarticBranch(_)))
    );
    let cone_patch = cone.trimmed((0.0, TAU), (0.5, 6.0)).unwrap();
    let cylinder_patch = cylinder.trimmed((0.0, TAU), (-2.0, 8.0)).unwrap();
    let result = cone_patch.intersect(&cylinder_patch, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
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

#[test]
fn test_a_meridian_clipped_between_the_poles_has_its_curve_in_the_sphere_parameters() {
    // The plane cuts the whole sphere along a meridian, which runs through
    // the poles and has no curve in the sphere's parameters; the half of
    // it on the half sphere runs from pole to pole and has one (#313).
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let half = sphere.trimmed((0.0, PI), (-0.5 * PI, 0.5 * PI)).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    let pieces: Vec<_> = [
        half.intersect(&plane, Tolerance::DEFAULT)
            .curves()
            .iter()
            .map(|piece| (piece.clone(), piece.first_uv_curve.clone()))
            .collect::<Vec<_>>(),
        plane
            .intersect(&half, Tolerance::DEFAULT)
            .curves()
            .iter()
            .map(|piece| (piece.clone(), piece.second_uv_curve.clone()))
            .collect(),
    ]
    .into_iter()
    .flatten()
    .collect();
    assert_eq!(pieces.len(), 2);
    for (piece, uv_curve) in pieces {
        let uv_curve = uv_curve.expect("the half meridian has a curve in parameters");
        let (start, end) = piece.parameter_range;
        for parameter in [start, 0.5 * (start + end), end] {
            let uv = uv_curve.point_at(parameter);
            assert!(
                half.basis()
                    .point_at(uv.x, uv.y)
                    .distance(piece.curve.point_at(parameter))
                    <= 1e-9
            );
            assert!(uv.y.abs() <= 0.5 * PI + 1e-8);
        }
    }
}

#[test]
fn test_a_marched_curve_ends_on_the_face_boundary() {
    // A flat ring, a revolution of a line trimmed to radii (0.5, 2), cut by
    // the plane x = 0.25: the plane meets the basis along a straight line,
    // marched, and the two pieces on the ring end on its circles, at
    // y = ±sqrt(0.5² − 0.25²) inside and ±sqrt(2² − 0.25²) outside (#314).
    let profile = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let revolution = RevolutionSurface::new(Curve::Line(profile), axis).unwrap();
    let ring = revolution.trimmed((0.0, TAU), (0.5, 2.0)).unwrap();
    let plane = Plane::new(Point::new(0.25, 0.0, 0.0), Vector::X).unwrap();
    let result = ring.intersect(&plane, TOL);
    assert_eq!(result.curves().len(), 2);
    let mut ends: Vec<f64> = result
        .curves()
        .iter()
        .flat_map(|piece| {
            let (start, end) = piece.parameter_range;
            [piece.curve.point_at(start).y, piece.curve.point_at(end).y]
        })
        .collect();
    ends.sort_by(f64::total_cmp);
    let (inner, outer) = ((0.5f64.powi(2) - 0.0625).sqrt(), (4.0f64 - 0.0625).sqrt());
    for (end, expected) in ends.iter().zip([-outer, -inner, inner, outer]) {
        assert!(
            (end - expected).abs() <= TOL.confusion,
            "{ends:?} against ±{inner}, ±{outer}"
        );
    }
    assert_inside(&result, &ring, &plane);
}

#[test]
fn test_a_closed_fitted_loop_keeps_a_piece_across_its_start_whole() {
    // A narrow cylinder across a wide one meets it in two closed loops,
    // fitted from a closed form. The short wall keeps two arcs of each:
    // the arc across a loop's start is one piece, its range running past
    // the loop's period, as an arc of a circle across its seam (#315).
    let wide = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let narrow = Cylinder::new(Point::new(0.0, 0.0, 2.5), Vector::X, 1.0).unwrap();
    let side = wide.trimmed((0.0, TAU), (0.0, 5.0)).unwrap();
    let short_wall = narrow.trimmed((0.0, TAU), (-1.9, 1.9)).unwrap();
    let result = side.intersect(&short_wall, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert_eq!(result.curves().len(), 4, "two arcs of each loop");
    assert_inside(&result, &side, &short_wall);
    assert_ends_on_boundaries(&result, &side, &short_wall);
    let across: Vec<&IntersectionCurve> = result
        .curves()
        .iter()
        .filter(|piece| piece.parameter_range.1 > piece.curve.bounds().1)
        .collect();
    assert_eq!(across.len(), 2, "one arc across each loop's start");
    for piece in across {
        let period = piece.curve.period().expect("a closed loop is periodic");
        let (start, end) = piece.parameter_range;
        assert!(end - start < period);
        assert!(TrimmedCurve::new(piece.curve.clone(), start, end).is_ok());
        // Its curves in parameters hold over the whole piece, past the
        // loop's start too.
        let uv_curves = [&piece.first_uv_curve, &piece.second_uv_curve];
        for (surface, uv_curve) in [&wide as &dyn ParametricSurface, &narrow]
            .into_iter()
            .zip(uv_curves)
        {
            let uv_curve = uv_curve.as_ref().expect("a fitted loop has its curves");
            for i in 0..=20 {
                let t = start + (end - start) * f64::from(i) / 20.0;
                let uv = uv_curve.point_at(t);
                let off = surface
                    .point_at(uv.x, uv.y)
                    .distance(piece.curve.point_at(t));
                assert!(off <= 1e-6, "{off} at {t} in {:?}", piece.parameter_range);
            }
        }
    }
    let untrimmed = wide.intersect(&narrow, TOL);
    assert_covers(&result, &untrimmed, &side, &short_wall);
}

/// A bicubic B-spline surface over `[0, 3]²`, its 6 × 6 poles at `(i, j)`
/// lifted to `height(i, j)`.
fn bspline_surface(height: impl Fn(f64, f64) -> f64) -> BSplineSurface {
    let poles = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| {
                    let (x, y) = (f64::from(i), f64::from(j));
                    Point::new(x, y, height(x, y))
                })
                .collect()
        })
        .collect();
    let (knots, multiplicities) = (vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4]);
    BSplineSurface::new(
        3,
        3,
        poles,
        knots.clone(),
        multiplicities.clone(),
        knots,
        multiplicities,
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_two_trimmed_bspline_surfaces_are_marched_and_clipped_to_the_ranges() {
    // Two wavy sheets with no closed form cross in several curves; the
    // march runs on the bases and its curves are clipped to both patches.
    let wavy = bspline_surface(|x, y| (0.7 * x * y).sin());
    let bumpy = bspline_surface(|x, y| 0.5 * (1.3 * x).cos() * (1.1 * y).cos());
    let a = wavy.clone().trimmed((0.4, 2.6), (0.3, 2.2)).unwrap();
    let b = bumpy.clone().trimmed((0.6, 2.8), (0.5, 2.7)).unwrap();
    let result = a.intersect(&b, TOL);
    assert_eq!(result.method(), IntersectionMethod::Numeric);
    assert_eq!(*result.completeness(), Completeness::Sampled);
    assert!(!result.curves().is_empty());
    assert_inside(&result, &a, &b);
    assert_ends_on_boundaries(&result, &a, &b);
    let untrimmed = wavy.intersect(&bumpy, TOL);
    assert_covers(&result, &untrimmed, &a, &b);
}

#[test]
fn test_a_marched_curve_whose_fit_strays_is_cut_by_samples() {
    // A wall of degree 1 along x with a tent at x = 50 whose two spans
    // are 5e-9 wide in parameter but 0.07 long in space: rounding the
    // parameter there moves the point by more than the confusion
    // distance, so the march's fit in the wall's parameters strays by more
    // than it, and the curve cannot be mapped to them within a tenth of
    // it either. Trimmed in the middle of the tent's first span and cut by
    // the plane z = 0.3, the curve is cut where samples of it leave the
    // range, at (50.025, 0.025).
    let (at, half, cluster) = (50.0, 0.05, 1e-8);
    let poles = [0.0, at, at + half, at + 2.0 * half, 100.0]
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let y = if i == 2 { half } else { 0.0 };
            vec![Point::new(x, y, -1.0), Point::new(x, y, 1.0)]
        })
        .collect();
    let knots = vec![
        0.0,
        at,
        at + 0.5 * cluster,
        at + cluster,
        100.0 - 2.0 * half + cluster,
    ];
    let wall = BSplineSurface::new(
        1,
        1,
        poles,
        knots,
        vec![2, 1, 1, 1, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap();
    let band = wall
        .clone()
        .trimmed((at + 0.25 * cluster, 60.0), (0.1, 0.9))
        .unwrap();
    let plane = Plane::new(Point::new(0.0, 0.0, 0.3), Vector::Z).unwrap();
    let result = band.intersect(&plane, TOL);
    assert_eq!(result.method(), IntersectionMethod::Numeric);
    assert_inside(&result, &band, &plane);
    let cut = Point::new(50.025, 0.025, 0.3);
    let across: Vec<&IntersectionCurve> = result
        .curves()
        .iter()
        .filter(|piece| {
            let (start, end) = piece.parameter_range;
            [start, end]
                .iter()
                .any(|&t| piece.curve.point_at(t).distance(cut) <= 1e-8)
        })
        .collect();
    assert_eq!(across.len(), 1, "one piece ends on the cut");
    assert!(across[0].deviation > TOL.confusion);
    let whole_plane = plane.trimmed((-200.0, 200.0), (-200.0, 200.0)).unwrap();
    let untrimmed = wall.intersect(&plane, TOL);
    assert_covers(&result, &untrimmed, &band, &whole_plane);
}
