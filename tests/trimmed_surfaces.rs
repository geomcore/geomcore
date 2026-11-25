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
    // The fit farther off than the confusion is reported, through the clip
    // (#311).
    assert!(
        matches!(
            result.completeness(),
            Completeness::Incomplete {
                reason: IncompleteReason::LooseFit,
                ..
            }
        ),
        "{:?}",
        result.completeness()
    );
    let whole_plane = plane.trimmed((-200.0, 200.0), (-200.0, 200.0)).unwrap();
    let untrimmed = wall.intersect(&plane, TOL);
    assert_covers(&result, &untrimmed, &band, &whole_plane);
}

// ------------------------------------------------- a curve and a trimmed surface

/// A half cylinder of radius 1 across the seam, `u` from 3π/2 to 5π/2
/// (the side `x ≥ 0`), from height 0 to 5.
fn half_cylinder_across_the_seam() -> TrimmedSurface<Cylinder> {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    cylinder.trimmed((1.5 * PI, 2.5 * PI), (0.0, 5.0)).unwrap()
}

#[test]
fn test_a_line_meets_a_face_across_the_seam_where_the_face_is() {
    let face = half_cylinder_across_the_seam();
    // Along x at height 1: it meets the cylinder at x = -1 and x = 1, and
    // the face only at x = 1, where the basis's angle 0 is 2π on the face.
    let line = Line::new(Point::new(-5.0, 0.0, 1.0), Vector::X).unwrap();
    assert_eq!(line.intersect(face.basis(), TOL).hits().len(), 2);
    let result = line.intersect(&face, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    let [hit] = result.hits() else {
        panic!("one hit on the face: {result:?}")
    };
    assert!((hit.curve_parameter - 6.0).abs() < 1e-12);
    let (u, v) = hit.surface_parameters;
    assert!((u - TAU).abs() < 1e-12 && (v - 1.0).abs() < 1e-12);
    assert!(face.point_at(u, v).distance(hit.point) < 1e-12);
    assert_eq!(face.intersect(&line, TOL), result);
}

#[test]
fn test_hits_outside_the_ranges_are_dropped() {
    let face = half_cylinder_across_the_seam();
    let above = Line::new(Point::new(-5.0, 0.0, 7.0), Vector::X).unwrap();
    let result = above.intersect(&face, TOL);
    assert!(result.is_empty());
    assert_eq!(result.method(), IntersectionMethod::Exact);
    // A trimmed line reaching only the far side of the cylinder.
    let line = Line::new(Point::new(-5.0, 0.0, 1.0), Vector::X).unwrap();
    let short = TrimmedCurve::new(line, 0.0, 5.0).unwrap();
    assert!(short.intersect(&face, TOL).is_empty());
    let long = TrimmedCurve::new(line, 0.0, 10.0).unwrap();
    assert_eq!(long.intersect(&face, TOL).hits().len(), 1);
    assert_eq!(face.intersect(&long, TOL), long.intersect(&face, TOL));
}

#[test]
fn test_a_hit_within_confusion_past_an_edge_is_kept_on_the_edge() {
    let face = half_cylinder_across_the_seam();
    // At height 5 + 1e-9, past the top edge by less than the confusion
    // distance.
    let line = Line::new(Point::new(-5.0, 0.0, 5.0 + 1e-9), Vector::X).unwrap();
    let result = line.intersect(&face, TOL);
    let [hit] = result.hits() else {
        panic!("the hit on the edge is kept: {result:?}")
    };
    assert_eq!(hit.surface_parameters.1, 5.0);
}

#[test]
fn test_a_line_through_a_pole_meets_a_face_holding_it() {
    // A quarter of the upper half of a sphere, holding its pole, where
    // every u is the same point.
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let face = sphere.trimmed((0.0, FRAC_PI_2), (0.0, FRAC_PI_2)).unwrap();
    let axis = Line::new(Point::new(0.0, 0.0, -5.0), Vector::Z).unwrap();
    let result = axis.intersect(&face, TOL);
    let [hit] = result.hits() else {
        panic!("the pole only: {result:?}")
    };
    assert!(hit.point.distance(Point::new(0.0, 0.0, 1.0)) < 1e-12);
}

#[test]
fn test_a_generator_line_on_a_face_overlaps_it_over_its_height() {
    let face = half_cylinder_across_the_seam();
    // The generator at angle 0, which the face holds as 2π.
    let generator = Line::new(Point::new(1.0, 0.0, -10.0), Vector::Z).unwrap();
    let result = generator.intersect(&face, TOL);
    assert!(result.hits().is_empty());
    let [overlap] = result.overlaps() else {
        panic!("one overlap: {result:?}")
    };
    let (start, end) = overlap.curve_range;
    assert!((start - 10.0).abs() < 1e-9 && (end - 15.0).abs() < 1e-9);
    // The generator on the other side of the cylinder misses the face.
    let behind = Line::new(Point::new(-1.0, 0.0, -10.0), Vector::Z).unwrap();
    assert!(behind.intersect(&face, TOL).is_empty());
}

#[test]
fn test_a_circle_of_a_face_across_the_seam_overlaps_it_in_one_piece() {
    let face = half_cylinder_across_the_seam();
    let circle = Circle::new(Point::new(0.0, 0.0, 2.0), Vector::Z, 1.0).unwrap();
    let result = circle.intersect(&face, TOL);
    let [overlap] = result.overlaps() else {
        panic!("one piece across the circle's seam: {result:?}")
    };
    let (start, end) = overlap.curve_range;
    assert!((start - 1.5 * PI).abs() < 1e-9 && (end - 2.5 * PI).abs() < 1e-9);
    assert!(TrimmedCurve::new(circle, start, end).is_ok());
    for i in 0..=20 {
        let t = start + (end - start) * f64::from(i) / 20.0;
        assert!(face.contains(circle.point_at(t), TOL));
    }
}

#[test]
fn test_a_line_in_a_square_overlaps_it_across_the_square() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let square = plane
        .trimmed((0.0, 1.0), (0.0, 1.0))
        .unwrap()
        .into_surface();
    let line = Curve::Line(Line::new(Point::new(-3.0, 0.5, 0.0), Vector::X).unwrap());
    let result = line.intersect(&square, TOL);
    let [overlap] = result.overlaps() else {
        panic!("one overlap: {result:?}")
    };
    let (start, end) = overlap.curve_range;
    assert!((start - 3.0).abs() < 1e-9 && (end - 4.0).abs() < 1e-9);
    let aside = Line::new(Point::new(-3.0, 2.0, 0.0), Vector::X).unwrap();
    assert!(aside.intersect(&square, TOL).is_empty());
}

#[test]
fn test_a_spline_curve_meets_a_face_numerically_and_is_clipped() {
    // An arch over the plane, crossing it at x = 0.5 and x = 2.5; the
    // square holds the first crossing only.
    let arch = BSplineCurve::new(
        2,
        vec![
            Point::new(0.0, 0.5, -1.0),
            Point::new(1.5, 0.5, 3.0),
            Point::new(3.0, 0.5, -1.0),
        ],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
    )
    .unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(arch.intersect(&plane, TOL).hits().len(), 2);
    let square = plane.trimmed((0.0, 1.0), (0.0, 1.0)).unwrap();
    let result = arch.intersect(&square, TOL);
    let [hit] = result.hits() else {
        panic!("one hit in the square: {result:?}")
    };
    assert!(hit.point.x < 1.0);
    assert!(square.contains(hit.point, TOL));
}

// ------------------------------------------------- a transformed trimmed surface

/// One patch over each kind of basis, a range across a seam among them.
fn patches_of_every_basis() -> Vec<TrimmedSurface> {
    let axis = Axis::new(Point::new(0.5, 0.0, 0.0), Vector::Z).unwrap();
    let slanted = Line::new(Point::new(2.0, 0.0, 0.0), Vector::new(1.0, 0.0, 2.0)).unwrap();
    let arc_profile = Circle::new(Point::new(1.0, 2.0, 0.0), Vector::Z, 1.5).unwrap();
    let patch = |surface: Surface, u_range, v_range| {
        TrimmedSurface::new(surface, u_range, v_range).unwrap()
    };
    vec![
        patch(
            Plane::new(Point::new(1.0, 2.0, 3.0), Vector::new(1.0, 1.0, 1.0))
                .unwrap()
                .into(),
            (-1.0, 2.0),
            (0.5, 3.0),
        ),
        patch(
            Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 2.0)
                .unwrap()
                .into(),
            (1.5 * PI, 2.5 * PI),
            (1.0, 4.0),
        ),
        patch(
            Cone::new(Point::ORIGIN, Vector::X, 0.3, 1.0)
                .unwrap()
                .into(),
            (0.5, 2.0),
            (-1.0, 3.0),
        ),
        patch(
            Sphere::new(Point::new(0.0, 1.0, 0.0), 2.0).unwrap().into(),
            (5.0, 7.0),
            (-0.5, 1.0),
        ),
        patch(
            Torus::new(Point::ORIGIN, Vector::Y, 3.0, 1.0)
                .unwrap()
                .into(),
            (0.5, 2.0),
            (5.5, 7.5),
        ),
        patch(
            ExtrusionSurface::new(Curve::Line(slanted), Vector::new(0.0, 1.0, 1.0))
                .unwrap()
                .into(),
            (-1.0, 2.0),
            (0.5, 2.5),
        ),
        patch(
            ExtrusionSurface::new(Curve::Circle(arc_profile), Vector::Z)
                .unwrap()
                .into(),
            (5.0, 7.0),
            (-1.0, 1.0),
        ),
        patch(
            RevolutionSurface::new(Curve::Line(slanted), axis)
                .unwrap()
                .into(),
            (5.0, 7.0),
            (0.5, 2.0),
        ),
        patch(
            bspline_surface(|x, y| (0.7 * x * y).sin()).into(),
            (0.4, 2.6),
            (0.3, 2.2),
        ),
    ]
}

#[test]
fn test_a_transformed_patch_holds_the_image_of_each_of_its_points() {
    let axis = Axis::new(Point::new(0.5, -1.0, 0.0), Vector::new(1.0, 2.0, 2.0)).unwrap();
    let mirror = Transform::mirror_plane(
        Frame::new(
            Point::new(0.1, 0.2, 0.3),
            Vector::new(1.0, 1.0, 0.0),
            Vector::new(0.0, 0.0, 1.0),
        )
        .unwrap(),
    );
    let transforms = [
        ("rotation", Transform::rotation(axis, 0.7)),
        ("mirror", mirror),
        (
            "point mirror",
            Transform::mirror_point(Point::new(1.0, 0.0, -1.0)),
        ),
        (
            "scaling",
            Transform::scaling(Point::new(1.0, -2.0, 0.5), 2.5),
        ),
        (
            "shrinking mirror",
            Transform::scaling(Point::ORIGIN, 0.4).then(mirror),
        ),
    ];
    for patch in patches_of_every_basis() {
        let ((u0, u1), (v0, v1)) = (patch.u_range(), patch.v_range());
        for (name, transform) in transforms {
            let image = patch.transformed(transform).unwrap();
            // 200 points across the patch, its edges included.
            for i in 0..20 {
                for j in 0..10 {
                    let u = u0 + (u1 - u0) * f64::from(i) / 19.0;
                    let v = v0 + (v1 - v0) * f64::from(j) / 9.0;
                    let moved = transform.apply_point(patch.point_at(u, v));
                    assert!(
                        image.contains(moved, TOL),
                        "{name} of {:?}: the image of ({u}, {v}) is outside",
                        patch.basis()
                    );
                }
            }
        }
    }
}

#[test]
fn test_a_scaled_patch_keeps_its_shape_in_parameters() {
    let face = half_cylinder_across_the_seam();
    let scaled = face
        .transformed(Transform::scaling(Point::ORIGIN, 2.0))
        .unwrap();
    // The angle is kept, the height scales.
    assert_eq!(scaled.u_range(), face.u_range());
    assert_eq!(scaled.v_range(), (0.0, 10.0));
    assert_eq!(
        face.transformed(Transform::scaling(Point::ORIGIN, 0.0)),
        Err(TransformError::Collapsed)
    );
}

// ------------------------------------------------- a trimmed surface trimmed again

#[test]
fn test_a_face_trimmed_again_covers_the_overlap_of_the_ranges() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let square = plane.trimmed((0.0, 2.0), (0.0, 2.0)).unwrap();
    let corner = square.trimmed((1.0, 5.0), (-3.0, 0.5)).unwrap();
    assert_eq!(
        (corner.u_range(), corner.v_range()),
        ((1.0, 2.0), (0.0, 0.5))
    );
    assert_eq!(
        square.trimmed((3.0, 5.0), (0.0, 1.0)),
        Err(TrimmedSurfaceConstructionError::DisjointRanges)
    );
    assert_eq!(
        square.trimmed((1.0, 1.0), (0.0, 1.0)),
        Err(TrimmedSurfaceConstructionError::EmptyRange)
    );
}

#[test]
fn test_a_face_trimmed_again_reads_a_periodic_range_modulo_the_period() {
    let face = half_cylinder_across_the_seam();
    // Given across the seam the other way round, one period lower.
    let lower = face.trimmed((-0.25 * PI, 0.25 * PI), (0.0, 5.0)).unwrap();
    assert!((lower.u_range().0 - 1.75 * PI).abs() < 1e-12);
    assert!((lower.u_range().1 - 2.25 * PI).abs() < 1e-12);
    // A whole period keeps the face's range.
    let whole = face.trimmed((0.0, TAU), (0.0, 5.0)).unwrap();
    assert_eq!(whole.u_range(), face.u_range());
    assert_eq!(
        face.trimmed((0.75 * PI, 1.25 * PI), (0.0, 5.0)),
        Err(TrimmedSurfaceConstructionError::DisjointRanges)
    );
    // Most of a turn, met by a range across its gap on both sides.
    let cylinder = *face.basis();
    let wide = cylinder.trimmed((0.0, 1.9 * PI), (0.0, 5.0)).unwrap();
    assert_eq!(
        wide.trimmed((1.5 * PI, 2.2 * PI), (0.0, 5.0)),
        Err(TrimmedSurfaceConstructionError::RangeInTwoPieces)
    );
    assert_eq!(
        face.trimmed((0.0, 7.0), (0.0, 5.0)),
        Err(TrimmedSurfaceConstructionError::ExceedsPeriod)
    );
}

/// Whether `x` lies in `range`, within `1e-9`, shifted by whole periods
/// when `period` is given.
fn in_range(x: f64, (start, end): (f64, f64), period: Option<f64>) -> bool {
    let x = match period {
        Some(period) => x - period * ((x - start + 1e-9) / period).floor(),
        None => x,
    };
    x >= start - 1e-9 && x <= end + 1e-9
}

/// Whether `uv` lies in the ranges of `patch`, a periodic parameter
/// shifted by whole periods when `modulo` is set.
fn in_ranges<S: ParametricSurface>(patch: &TrimmedSurface<S>, uv: Point2D, modulo: bool) -> bool {
    let basis = patch.basis();
    let (u_period, v_period) = if modulo {
        (basis.u_period(), basis.v_period())
    } else {
        (None, None)
    };
    in_range(uv.x, patch.u_range(), u_period) && in_range(uv.y, patch.v_range(), v_period)
}

/// The one curve of `result`, an edge shared by `a` and `b`: exact,
/// tangent, its curve in each patch's parameters on the curve and inside
/// that patch's ranges, its middle as they are written, no period away.
fn shared_edge<A: ParametricSurface, B: ParametricSurface>(
    result: &SurfaceSurfaceIntersection,
    a: &TrimmedSurface<A>,
    b: &TrimmedSurface<B>,
) -> IntersectionCurve {
    assert!(!result.is_coincident(), "{result:?}");
    assert!(result.contact_points().is_empty(), "{result:?}");
    let [edge] = result.curves() else {
        panic!("one shared edge: {result:?}")
    };
    assert!(edge.tangent);
    assert!(edge.deviation <= 1e-7);
    let (start, end) = edge.parameter_range;
    let first = edge.first_uv_curve.as_ref().unwrap();
    let second = edge.second_uv_curve.as_ref().unwrap();
    for i in 0..=20 {
        let t = start + (end - start) * f64::from(i) / 20.0;
        let point = edge.curve.point_at(t);
        let (on_a, on_b) = (first.point_at(t), second.point_at(t));
        assert!(a.point_at(on_a.x, on_a.y).distance(point) < 1e-12);
        assert!(b.point_at(on_b.x, on_b.y).distance(point) <= edge.deviation + 1e-12);
        let modulo = i != 10;
        assert!(in_ranges(a, on_a, modulo) && in_ranges(b, on_b, modulo));
    }
    edge.clone()
}

#[test]
fn test_stacked_faces_of_a_cylinder_meet_along_their_circle() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let low = cylinder.trimmed((0.0, TAU), (0.0, 1.0)).unwrap();
    let high = cylinder.trimmed((0.0, TAU), (1.0, 2.0)).unwrap();
    let result = low.intersect(&high, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    let edge = shared_edge(&result, &low, &high);
    let Curve::Circle(circle) = &edge.curve else {
        panic!("a circle: {edge:?}")
    };
    assert!((circle.center().z - 1.0).abs() < 1e-12 && (circle.radius() - 1.0).abs() < 1e-12);
    assert_eq!(edge.parameter_range, (0.0, TAU));
    // In the other order, the curves in parameters swap.
    let swapped = shared_edge(&high.intersect(&low, TOL), &high, &low);
    assert_eq!(swapped.first_uv_curve, edge.second_uv_curve);
    assert_eq!(swapped.second_uv_curve, edge.first_uv_curve);
}

#[test]
fn test_adjacent_squares_of_a_plane_meet_along_their_side() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let left = plane.trimmed((0.0, 1.0), (0.0, 1.0)).unwrap();
    let right = plane.trimmed((1.0, 2.0), (0.0, 1.0)).unwrap();
    let edge = shared_edge(&left.intersect(&right, TOL), &left, &right);
    let (start, end) = edge.parameter_range;
    let ends = [edge.curve.point_at(start), edge.curve.point_at(end)];
    assert!(ends[0].distance(Point::new(1.0, 0.0, 0.0)) < 1e-12);
    assert!(ends[1].distance(Point::new(1.0, 1.0, 0.0)) < 1e-12);
    // A side shared in part: only the common piece.
    let lower_right = plane.trimmed((1.0, 2.0), (-0.5, 0.5)).unwrap();
    let edge = shared_edge(&left.intersect(&lower_right, TOL), &left, &lower_right);
    assert_eq!(edge.parameter_range, (0.0, 0.5));
}

#[test]
fn test_an_edge_across_the_seam_is_one_curve_in_each_faces_parameters() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    // The same half, written below the seam on one face, above it on the
    // other.
    let low = cylinder
        .trimmed((-FRAC_PI_2, FRAC_PI_2), (0.0, 1.0))
        .unwrap();
    let high = cylinder.trimmed((1.5 * PI, 2.5 * PI), (1.0, 2.0)).unwrap();
    let edge = shared_edge(&low.intersect(&high, TOL), &low, &high);
    assert_eq!(edge.parameter_range, (-FRAC_PI_2, FRAC_PI_2));
    // Two halves side by side share their two generators, one at the
    // seam of the second.
    let east = cylinder.trimmed((1.5 * PI, 2.5 * PI), (0.0, 1.0)).unwrap();
    let west = cylinder.trimmed((FRAC_PI_2, 1.5 * PI), (0.0, 1.0)).unwrap();
    let result = east.intersect(&west, TOL);
    assert!(!result.is_coincident() && result.contact_points().is_empty());
    let mut feet: Vec<Point> = result
        .curves()
        .iter()
        .map(|edge| {
            assert!(edge.tangent);
            let Curve::Line(line) = &edge.curve else {
                panic!("a generator: {edge:?}")
            };
            assert_eq!(edge.parameter_range, (0.0, 1.0));
            line.point_at(0.0)
        })
        .collect();
    feet.sort_by(|p, q| p.y.total_cmp(&q.y));
    assert_eq!(feet.len(), 2);
    assert!(feet[0].distance(Point::new(0.0, -1.0, 0.0)) < 1e-12);
    assert!(feet[1].distance(Point::new(0.0, 1.0, 0.0)) < 1e-12);
}

#[test]
fn test_faces_meeting_at_a_corner_touch_at_a_point() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let a = plane.trimmed((0.0, 1.0), (0.0, 1.0)).unwrap();
    let b = plane.trimmed((1.0, 2.0), (1.0, 2.0)).unwrap();
    let result = a.intersect(&b, TOL);
    assert!(!result.is_coincident() && result.curves().is_empty());
    assert_eq!(result.contact_points(), [Point::new(1.0, 1.0, 0.0)]);
    // Two quarters of a cylinder, diagonal across height 1, touch at the
    // two corners where their angles meet.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let low = cylinder.trimmed((0.0, PI), (0.0, 1.0)).unwrap();
    let high = cylinder.trimmed((PI, TAU), (1.0, 2.0)).unwrap();
    let result = low.intersect(&high, TOL);
    assert!(!result.is_coincident() && result.curves().is_empty());
    let mut points = result.contact_points().to_vec();
    points.sort_by(|p, q| p.x.total_cmp(&q.x));
    assert_eq!(points.len(), 2);
    assert!(points[0].distance(Point::new(-1.0, 0.0, 1.0)) < 1e-12);
    assert!(points[1].distance(Point::new(1.0, 0.0, 1.0)) < 1e-12);
}

#[test]
fn test_faces_a_hair_apart_or_overlapping_by_a_hair_share_their_edge() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let left = plane.trimmed((0.0, 1.0), (0.0, 1.0)).unwrap();
    for gap in [-5e-8, 5e-8] {
        let right = plane.trimmed((1.0 + gap, 2.0), (0.0, 1.0)).unwrap();
        let edge = shared_edge(&left.intersect(&right, TOL), &left, &right);
        // Reported on the first face's side: a gap is its deviation from
        // the second, which a hair of overlap does not have.
        assert!((edge.curve.point_at(0.0).x - 1.0).abs() < 1e-15);
        assert!((edge.deviation - gap.max(0.0)).abs() < 1e-15, "{edge:?}");
    }
    // Farther apart than the confusion distance: nothing.
    let right = plane.trimmed((1.0 + 1e-6, 2.0), (0.0, 1.0)).unwrap();
    assert!(left.intersect(&right, TOL).is_empty());
}

#[test]
fn test_faces_overlapping_with_area_stay_coincident() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let left = plane.trimmed((0.0, 1.0), (0.0, 1.0)).unwrap();
    // A strip a thousandth wide, which no grid of samples would see.
    let strip = plane.trimmed((0.999, 2.0), (0.4, 0.401)).unwrap();
    let result = left.intersect(&strip, TOL);
    assert!(result.is_coincident() && result.curves().is_empty());
    // Overlapping, and sharing a further edge across the seam.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let most = cylinder.trimmed((0.0, 1.5 * PI), (0.0, 1.0)).unwrap();
    let rest = cylinder.trimmed((PI, TAU), (0.0, 1.0)).unwrap();
    let result = most.intersect(&rest, TOL);
    assert!(result.is_coincident() && result.curves().is_empty());
    // A face across the seam inside a whole turn.
    let around = cylinder.trimmed((0.0, TAU), (-1.0, 1.0)).unwrap();
    let across = cylinder.trimmed((1.5 * PI, 2.5 * PI), (0.0, 0.5)).unwrap();
    assert!(around.intersect(&across, TOL).is_coincident());
}

#[test]
fn test_faces_of_one_cylinder_built_on_other_frames_share_their_edge() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let low = cylinder.trimmed((0.0, TAU), (0.0, 1.0)).unwrap();
    // The axis reversed, the frame turned and moved along it: the height
    // on it is 2 - z, and its angle runs the other way.
    let frame = Frame::new(
        Point::new(0.0, 0.0, 2.0),
        -Vector::Z,
        Vector::new(1.0, 1.0, 0.0),
    )
    .unwrap();
    let other = Cylinder::from_frame(frame, 1.0).unwrap();
    let high = other.trimmed((0.0, TAU), (0.0, 1.0)).unwrap();
    let edge = shared_edge(&low.intersect(&high, TOL), &low, &high);
    assert_eq!(edge.parameter_range, (0.0, TAU));
    // Half of each: they share the eighth of the circle from angle 0 to
    // π/4, where the other's angle 0 lies.
    let low_half = cylinder.trimmed((0.0, PI), (0.0, 1.0)).unwrap();
    let high_half = other.trimmed((0.0, PI), (0.0, 1.0)).unwrap();
    let edge = shared_edge(&low_half.intersect(&high_half, TOL), &low_half, &high_half);
    let (start, end) = edge.parameter_range;
    assert!((end - start - 0.25 * PI).abs() < 1e-9, "{edge:?}");
    // Overlapping with area, on the other frame: coincident.
    let overlapping = other.trimmed((0.0, TAU), (0.0, 1.5)).unwrap();
    assert!(low.intersect(&overlapping, TOL).is_coincident());
}

#[test]
fn test_squares_of_one_plane_on_a_frame_turned_a_quarter_share_their_side() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let left = plane.trimmed((0.0, 1.0), (0.0, 1.0)).unwrap();
    // u along y and v along -x: the square x in (1, 2), y in (0, 1).
    let frame = Frame::new(Point::ORIGIN, Vector::Z, Vector::Y).unwrap();
    let turned = Plane::from_frame(frame);
    let right = turned.trimmed((0.0, 1.0), (-2.0, -1.0)).unwrap();
    let edge = shared_edge(&left.intersect(&right, TOL), &left, &right);
    assert_eq!(edge.parameter_range, (0.0, 1.0));
}

#[test]
fn test_faces_of_one_spline_surface_share_their_edge() {
    let sheet = bspline_surface(|x, y| (0.7 * x * y).sin());
    let near = sheet.clone().trimmed((0.0, 1.5), (0.5, 2.5)).unwrap();
    let far = sheet.trimmed((1.5, 3.0), (0.0, 2.0)).unwrap();
    let result = near.intersect(&far, TOL);
    let edge = shared_edge(&result, &near, &far);
    assert_eq!(edge.parameter_range, (0.5, 2.0));
}

/// A deterministic sequence of numbers spread evenly over `[-1, 1)`.
struct Noise(u64);

impl Noise {
    fn next(&mut self) -> f64 {
        // xorshift64
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    }

    fn vector(&mut self) -> Vector {
        Vector::new(self.next(), self.next(), self.next())
    }
}

/// How a plane cuts a sphere: across its axis at a slant, through its
/// centre holding the axis, or through one pole.
#[derive(Clone, Copy, Debug)]
enum Cut {
    Oblique,
    Meridian,
    ThroughOnePole,
}

/// A sphere in a tilted frame from `noise` and a plane cutting it as
/// `cut` says, with the distance from the centre to the plane.
fn sphere_and_cutting_plane(noise: &mut Noise, cut: Cut) -> (Sphere, Plane, f64) {
    let frame = loop {
        if let Ok(frame) = Frame::new(
            Point::ORIGIN + 3.0 * noise.vector(),
            noise.vector(),
            noise.vector(),
        ) {
            break frame;
        }
    };
    let radius = 2.0 + noise.next();
    let sphere = Sphere::from_frame(frame, radius).unwrap();
    let (x, y, z) = (
        Vector::from(frame.x_direction()),
        Vector::from(frame.y_direction()),
        Vector::from(frame.z_direction()),
    );
    let angle = PI * noise.next();
    let across_axis = x * angle.cos() + y * angle.sin();
    let (point, normal) = match cut {
        Cut::Oblique => (
            frame.origin() + 0.6 * radius * noise.next() * z,
            z + 0.5 * (1.5 + noise.next()) * across_axis,
        ),
        Cut::Meridian => (frame.origin(), across_axis),
        Cut::ThroughOnePole => (
            frame.origin() + radius * z,
            across_axis + 0.5 * noise.next() * z,
        ),
    };
    let plane = Plane::new(point, normal).unwrap();
    let normal = normal / normal.magnitude();
    let distance = (frame.origin() - point).dot(normal).abs();
    (sphere, plane, distance)
}

#[test]
fn test_a_plane_patch_cuts_a_whole_sphere_patch_along_one_closed_circle() {
    // Whatever the plane's direction, through a pole or both, the section
    // is one whole circle, cut at the sphere's seam and joined across it
    // (#509).
    let mut noise = Noise(0x5EED_0509);
    for cut in [Cut::Oblique, Cut::Meridian, Cut::ThroughOnePole] {
        for _ in 0..12 {
            let (sphere, plane, distance) = sphere_and_cutting_plane(&mut noise, cut);
            let radius = sphere.radius();
            let whole = sphere.trimmed((0.0, TAU), (-FRAC_PI_2, FRAC_PI_2)).unwrap();
            let patch = plane
                .trimmed((-3.0 * radius, 3.0 * radius), (-3.0 * radius, 3.0 * radius))
                .unwrap();
            let result = patch.intersect(&whole, TOL);
            let context = format!("{cut:?} {sphere:?} {plane:?}");
            assert_eq!(result.method(), IntersectionMethod::Exact, "{context}");
            assert!(result.is_complete(), "{context}");
            assert!(result.contact_points().is_empty(), "{context}");
            assert_eq!(result.curves().len(), 1, "{context}");
            let piece = &result.curves()[0];
            assert!(matches!(piece.curve, Curve::Circle(_)), "{context}");
            let (start, end) = piece.parameter_range;
            assert!(
                (end - start - TAU).abs() <= 1e-9,
                "{context}: {start} {end}"
            );
            let points = samples(piece, 4001);
            for &point in &points {
                let off_sphere = (point.distance(sphere.center()) - radius).abs();
                let off_plane = plane.project_point(point, TOL).distance;
                assert!(
                    off_sphere <= 1e-9 && off_plane <= 1e-9,
                    "{context}: {point:?}"
                );
            }
            let length: f64 = points
                .windows(2)
                .map(|pair| pair[0].distance(pair[1]))
                .sum();
            let expected = TAU * (radius * radius - distance * distance).sqrt();
            assert!(
                (length - expected).abs() <= 1e-6 * expected,
                "{context}: length {length}, expected {expected}"
            );
            assert_covers(&result, &plane.intersect(&sphere, TOL), &patch, &whole);
            // Either order gives the same piece.
            let swapped = whole.intersect(&patch, TOL);
            assert_eq!(swapped.curves().len(), 1, "{context}");
            assert_eq!(swapped.curves()[0].parameter_range, piece.parameter_range);
        }
    }
}

#[test]
fn test_a_plane_cuts_a_sphere_face_into_arcs_ending_on_its_boundary() {
    // A face of the sphere short of its seam and of both poles: the circle
    // is cut where it leaves the face, through the planes of its sides.
    let mut noise = Noise(0x0509_FACE);
    for cut in [Cut::Oblique, Cut::Meridian, Cut::ThroughOnePole] {
        let mut pieces = 0;
        for _ in 0..12 {
            let (sphere, plane, _) = sphere_and_cutting_plane(&mut noise, cut);
            let face = sphere.trimmed((0.5, 4.5), (-1.2, 1.3)).unwrap();
            let result = plane.intersect(&face, TOL);
            assert!(result.is_complete(), "{cut:?} {sphere:?} {plane:?}");
            assert_inside(&result, &plane, &face);
            for piece in result.curves() {
                let (start, end) = piece.parameter_range;
                for point in [piece.curve.point_at(start), piece.curve.point_at(end)] {
                    assert!(
                        on_boundary(&face, point),
                        "{cut:?} {point:?} inside the face"
                    );
                }
            }
            let untrimmed = plane.intersect(&sphere, TOL);
            let everywhere = plane.trimmed((-1e3, 1e3), (-1e3, 1e3)).unwrap();
            assert_covers(&result, &untrimmed, &everywhere, &face);
            pieces += result.curves().len();
        }
        assert!(pieces >= 12, "{cut:?}: {pieces} arcs");
    }
}
