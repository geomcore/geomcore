//! The routed intersection: exact when a closed form exists, numeric
//! otherwise, and the result says which method ran.

use geomcore::{Curve, Surface, TrimmedCurve, TrimmedCurve2D, *};
use std::f64::consts::{PI, TAU};

const TOL: Tolerance = Tolerance::DEFAULT;

fn assert_hits_lie_on<C: ParametricCurve, S: ParametricSurface>(
    r: &CurveSurfaceIntersection,
    c: &C,
    s: &S,
) {
    for h in r.hits() {
        assert!(
            c.point_at(h.curve_parameter).distance(h.point) <= 1e-7,
            "curve parameter does not give the point"
        );
        let (u, v) = h.surface_parameters;
        assert!(
            s.point_at(u, v).distance(h.point) <= 1e-7,
            "surface parameters do not give the point"
        );
    }
}

#[test]
fn test_line_sphere_routes_exact_and_leaves_surface_parameters_to_the_caller() {
    let line = Line::new(Point::new(0.0, 0.0, -5.0), Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let r = line.intersect(&sphere, TOL);
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert!(r.is_complete());
    assert_eq!(r.hits().len(), 2); // through both poles
    assert_hits_lie_on(&r, &line, &sphere);
    for h in r.hits() {
        // The line parameter comes out of the solve; the sphere's are
        // recovered from the point, at a pole as well.
        assert_eq!(
            h.surface_parameters,
            sphere.parameters_of(h.point, Tolerance::DEFAULT)
        );
    }
}

#[test]
fn test_reverse_order_reads_the_same() {
    let line = Line::new(Point::new(-5.0, 0.5, 0.0), Vector::X).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    assert_eq!(
        line.intersect(&cylinder, TOL),
        cylinder.intersect(&line, TOL)
    );
}

#[test]
fn test_routed_matches_exact_hits() {
    let circle = Circle::new(Point::new(0.0, 0.0, 0.5), Vector::X, 2.0).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let routed = circle.intersect(&sphere, TOL);
    let CircleSurfaceIntersection::Hits(exact) = circle.intersect_exact(&sphere, TOL) else {
        panic!("expected hits")
    };
    let routed_parameters: Vec<f64> = routed.hits().iter().map(|h| h.curve_parameter).collect();
    let exact_parameters: Vec<f64> = exact.iter().map(|h| h.parameter).collect();
    assert_eq!(routed_parameters, exact_parameters);
}

#[test]
fn test_line_in_plane_is_an_unbounded_overlap() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let r = line.intersect(&plane, TOL);
    assert!(r.hits().is_empty());
    assert_eq!(
        r.overlaps(),
        &[CurveSurfaceOverlap {
            curve_range: (f64::NEG_INFINITY, f64::INFINITY)
        }]
    );
    assert!(!r.is_empty());
}

#[test]
fn test_circle_on_sphere_is_a_full_period_overlap() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let r = circle.intersect(&sphere, TOL);
    assert_eq!(
        r.overlaps(),
        &[CurveSurfaceOverlap {
            curve_range: (0.0, TAU)
        }]
    );
}

#[test]
fn test_circle_torus_routes_exact_and_complete() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let circle = Circle::new(Point::new(3.0, 0.0, 0.0), Vector::Z, 1.5).unwrap();
    let r = circle.intersect(&torus, TOL);
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert!(r.is_complete());
    assert_eq!(r.hits().len(), 4);
    assert_hits_lie_on(&r, &circle, &torus);
    // The same pair the other way round.
    let swapped = torus.intersect(&circle, TOL);
    assert_eq!(swapped.method(), IntersectionMethod::Exact);
    assert_eq!(swapped.hits().len(), 4);
    // A tube circle is one overlap over the whole circle.
    let tube = torus.u_iso(1.0).unwrap();
    let on = tube.intersect(&torus, TOL);
    assert_eq!(on.method(), IntersectionMethod::Exact);
    assert_eq!(
        on.overlaps(),
        &[CurveSurfaceOverlap {
            curve_range: (0.0, TAU)
        }]
    );
}

#[test]
fn test_ellipse_torus_routes_numeric() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let ellipse = Ellipse::from_frame(Frame::WORLD, 5.0, 1.0).unwrap();
    let r = ellipse.intersect(&torus, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert!(!r.hits().is_empty());
    assert_hits_lie_on(&r, &ellipse, &torus);
}

#[test]
fn test_bspline_through_enum_routes_numeric() {
    let circle = Circle::new(Point::ORIGIN, Vector::X, 2.0).unwrap();
    let curve = Curve::from(circle.to_bspline());
    let plane = Surface::from(Plane::new(Point::ORIGIN, Vector::Z).unwrap());
    let r = curve.intersect(&plane, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert_eq!(r.hits().len(), 2);
    assert_hits_lie_on(&r, &curve, &plane);
}

#[test]
fn test_trimmed_arc_across_the_seam_keeps_only_its_hits() {
    // A unit circle in z = 0 meets the plane x = cos(0.3) at ±0.3 (≡ 2π - 0.3).
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let plane = Plane::new(Point::new(0.3f64.cos(), 0.0, 0.0), Vector::X).unwrap();
    let arc = TrimmedCurve::new(circle, -1.0, 1.0).unwrap();
    let r = arc.intersect(&plane, TOL);
    let mut parameters: Vec<f64> = r.hits().iter().map(|h| h.curve_parameter).collect();
    parameters.sort_by(f64::total_cmp);
    assert_eq!(parameters.len(), 2);
    assert!((parameters[0] + 0.3).abs() < 1e-9 && (parameters[1] - 0.3).abs() < 1e-9);
    let far = TrimmedCurve::new(circle, 1.0, PI).unwrap();
    assert!(far.intersect(&plane, TOL).is_empty());
}

#[test]
fn test_tangent_hit_has_multiplicity_two() {
    let line = Line::new(Point::new(-5.0, 1.0, 0.0), Vector::X).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let r = line.intersect(&cylinder, TOL);
    assert_eq!(r.hits().len(), 1);
    assert!(r.hits()[0].multiplicity >= 2);
}

fn assert_curve_hits_lie_on<A: ParametricCurve, B: ParametricCurve>(
    r: &CurveCurveIntersection,
    a: &A,
    b: &B,
) {
    for h in r.hits() {
        assert!(a.point_at(h.first_parameter).distance(h.point) <= 1e-7);
        assert!(b.point_at(h.second_parameter).distance(h.point) <= 1e-7);
    }
}

/// Every overlap range is monotone and lies in its curve's range.
fn assert_overlaps_within(r: &CurveCurveIntersection, first: (f64, f64), second: (f64, f64)) {
    let inside = |(low, high): (f64, f64), (start, end): (f64, f64)| {
        low <= high && low >= start - 1e-12 && high <= end + 1e-12
    };
    for o in r.overlaps() {
        assert!(
            inside(o.first_range, first),
            "{:?} not in {:?}",
            o.first_range,
            first
        );
        assert!(
            inside(o.second_range, second),
            "{:?} not in {:?}",
            o.second_range,
            second
        );
    }
}

fn close(a: (f64, f64), b: (f64, f64)) -> bool {
    (a.0 - b.0).abs() < 1e-9 && (a.1 - b.1).abs() < 1e-9
}

#[test]
fn test_circle_line_reads_from_the_receiver() {
    let line = Line::new(Point::new(-5.0, 0.0, 0.0), Vector::X).unwrap();
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let r = circle.intersect(&line, TOL);
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert_curve_hits_lie_on(&r, &circle, &line);
    let forward = line.intersect(&circle, TOL);
    assert_eq!(r.hits().len(), forward.hits().len());
    for y in forward.hits() {
        assert!(
            r.hits()
                .iter()
                .any(|x| x.first_parameter == y.second_parameter
                    && x.second_parameter == y.first_parameter)
        );
    }
}

#[test]
fn test_coincident_lines_overlap_everywhere() {
    let a = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let r = a.intersect(&a, TOL);
    assert_eq!(r.overlaps().len(), 1);
    assert_eq!(
        r.overlaps()[0].first_range,
        (f64::NEG_INFINITY, f64::INFINITY)
    );
}

#[test]
fn test_overlapping_segments_clip_to_both_ranges() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let a = TrimmedCurve::new(line, 0.0, 2.0).unwrap();
    let b = TrimmedCurve::new(line, 1.0, 3.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert_eq!(r.overlaps().len(), 1);
    let o = r.overlaps()[0];
    assert!((o.first_range.0 - 1.0).abs() < 1e-12 && (o.first_range.1 - 2.0).abs() < 1e-12);
    assert!((o.second_range.0 - 1.0).abs() < 1e-12 && (o.second_range.1 - 2.0).abs() < 1e-12);
    assert!(o.same_direction);
}

#[test]
fn test_bspline_circle_routes_numeric() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let other = Circle::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    let bspline = other.to_bspline();
    let r = circle.intersect(&bspline, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert_eq!(r.hits().len(), 2);
    assert_curve_hits_lie_on(&r, &circle, &bspline);
}

#[test]
fn test_skew_lines_are_empty_and_complete() {
    let a = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let b = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.is_empty() && r.is_complete());
}

#[test]
fn test_coincident_arcs_overlap_on_both_sides_of_the_seam() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let a = TrimmedCurve::new(circle, 0.0, 5.0).unwrap();
    let b = TrimmedCurve::new(circle, 4.0, 7.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert_eq!(r.overlaps().len(), 2);
    assert_overlaps_within(&r, (0.0, 5.0), (4.0, 7.0));
    let o = r.overlaps();
    assert!(close(o[0].first_range, (0.0, 7.0 - TAU)));
    assert!(close(o[0].second_range, (TAU, 7.0)));
    assert!(close(o[1].first_range, (4.0, 5.0)));
    assert!(close(o[1].second_range, (4.0, 5.0)));
    assert!(o.iter().all(|piece| piece.same_direction));
    // The other order reads the same overlaps from b's side.
    let reverse = b.intersect(&a, TOL);
    assert_eq!(reverse.overlaps().len(), 2);
    assert_overlaps_within(&reverse, (4.0, 7.0), (0.0, 5.0));
}

#[test]
fn test_coincident_lines_with_opposite_directions() {
    let a = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let b = Line::new(Point::new(5.0, 0.0, 0.0), -Vector::X).unwrap();
    // a covers x in [0, 4], b covers x in [2, 5].
    let first = TrimmedCurve::new(a, 0.0, 4.0).unwrap();
    let second = TrimmedCurve::new(b, 0.0, 3.0).unwrap();
    let r = first.intersect(&second, TOL);
    assert_eq!(r.overlaps().len(), 1);
    assert!(close(r.overlaps()[0].first_range, (2.0, 4.0)));
    assert!(close(r.overlaps()[0].second_range, (1.0, 3.0)));
    assert!(!r.overlaps()[0].same_direction);
    // Unbounded and opposite: both ranges are the whole line.
    let whole = a.intersect(&b, TOL);
    assert_eq!(whole.overlaps().len(), 1);
    let o = whole.overlaps()[0];
    assert_eq!(o.first_range, (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(o.second_range, (f64::NEG_INFINITY, f64::INFINITY));
    assert!(!o.same_direction);
}

#[test]
fn test_coincident_full_circles_with_different_frames() {
    let a = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let rotated = Frame::new(Point::ORIGIN, Vector::Z, Vector::new(1.0, 1.0, 0.0)).unwrap();
    let flipped = Frame::new(Point::ORIGIN, -Vector::Z, Vector::Y).unwrap();
    for (frame, same_direction) in [(rotated, true), (flipped, false)] {
        let b = Circle::from_frame(frame, 2.0).unwrap();
        let r = a.intersect(&b, TOL);
        assert_eq!(r.method(), IntersectionMethod::Exact);
        assert_eq!(r.overlaps().len(), 1, "{frame:?}");
        let o = r.overlaps()[0];
        assert_eq!(o.same_direction, same_direction, "{frame:?}");
        assert!(close(o.first_range, a.bounds()));
        assert!((o.second_range.1 - o.second_range.0 - TAU).abs() < 1e-9);
        assert_overlaps_within(&r, a.bounds(), b.bounds());
    }
}

#[test]
fn test_trimmed_curves_keep_only_their_hits() {
    let line = Line::new(Point::new(-5.0, 0.0, 0.0), Vector::X).unwrap();
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    // The arc across the seam keeps the hit at angle 0, not the one at PI.
    let arc = TrimmedCurve::new(circle, -1.0, 1.0).unwrap();
    let r = line.intersect(&arc, TOL);
    assert_eq!(r.hits().len(), 1);
    assert!((r.hits()[0].first_parameter - 7.0).abs() < 1e-12);
    // The trimmed arc stores its angle, which agrees with the circle's own.
    let angle = r.hits()[0].second_parameter;
    assert!(angle.abs() < 1e-12);
    assert!((angle - circle.parameter_of(r.hits()[0].point, Tolerance::DEFAULT)).abs() < 1e-12);
    assert_curve_hits_lie_on(&r, &line, &arc);
    // A segment that stops before the circle meets it once, at PI.
    let segment = TrimmedCurve::new(line, 0.0, 5.0).unwrap();
    let r = circle.intersect(&segment, TOL);
    assert_eq!(r.hits().len(), 1);
    let hit = r.hits()[0];
    assert!((hit.first_parameter - PI).abs() < 1e-12);
    assert!((hit.second_parameter - 3.0).abs() < 1e-12);
    assert_curve_hits_lie_on(&r, &circle, &segment);
}

#[test]
fn test_routed_curve_curve_matches_exact_hits() {
    let line = Line::new(Point::new(-5.0, 0.5, 0.0), Vector::X).unwrap();
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.0).unwrap();
    let routed = Curve::from(line).intersect(&Curve::from(ellipse), TOL);
    let LineConicIntersection::Points(first, second) = line.intersect_exact(&ellipse, TOL) else {
        panic!("expected two hits");
    };
    let routed: Vec<_> = routed
        .hits()
        .iter()
        .map(|h| {
            (
                h.first_parameter,
                h.point,
                h.second_parameter,
                h.multiplicity,
            )
        })
        .collect();
    let exact = [first, second].map(|h| (h.line_parameter, h.point, h.conic_parameter, 1));
    assert_eq!(routed, exact);
}

#[test]
fn test_arc_across_the_seam_of_a_full_circle_splits_there() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let arc = TrimmedCurve::new(circle, -1.0, 1.0).unwrap();
    let r = arc.intersect(&circle, TOL);
    assert_eq!(r.overlaps().len(), 2);
    assert_overlaps_within(&r, (-1.0, 1.0), circle.bounds());
    let o = r.overlaps();
    assert!(close(o[0].first_range, (-1.0, 0.0)));
    assert!(close(o[0].second_range, (TAU - 1.0, TAU)));
    assert!(close(o[1].first_range, (0.0, 1.0)));
    assert!(close(o[1].second_range, (0.0, 1.0)));
    // Read from the full circle: one overlap per side of the arc's seam.
    let reverse = circle.intersect(&arc, TOL);
    assert_eq!(reverse.overlaps().len(), 2);
    assert_overlaps_within(&reverse, circle.bounds(), (-1.0, 1.0));
}

fn assert_hits_2d_lie_on<A: ParametricCurve2D, B: ParametricCurve2D>(
    r: &CurveCurveIntersection2D,
    first: &A,
    second: &B,
) {
    for h in r.hits() {
        assert!(first.point_at(h.first_parameter).distance(h.point) <= 1e-7);
        assert!(second.point_at(h.second_parameter).distance(h.point) <= 1e-7);
    }
}

#[test]
fn test_2d_line_circle_hits() {
    let line = Line2D::new(Point2D::new(-5.0, 0.0), Vector2D::new(1.0, 0.0)).unwrap();
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let r = circle.intersect(&line, TOL);
    assert_eq!(r.hits().len(), 2);
    assert_hits_2d_lie_on(&r, &circle, &line);
    // The line first reads the same hits from its own side.
    let from_line = line.intersect(&circle, TOL);
    assert_eq!(from_line.hits().len(), 2);
    assert_hits_2d_lie_on(&from_line, &line, &circle);
    assert!((from_line.hits()[0].first_parameter - 3.0).abs() < 1e-9);
}

#[test]
fn test_2d_trimmed_segment_misses() {
    let line = Line2D::new(Point2D::new(-5.0, 0.0), Vector2D::new(1.0, 0.0)).unwrap();
    let segment = TrimmedCurve2D::new(line, 0.0, 1.0).unwrap();
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    assert!(segment.intersect(&circle, TOL).is_empty());
    assert!(circle.intersect(&segment, TOL).is_empty());
}

#[test]
fn test_2d_tangent_circles() {
    let a = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    let b = Circle2D::new(Point2D::new(2.0, 0.0), 1.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert_eq!(r.hits().len(), 1);
    assert_eq!(r.hits()[0].multiplicity, 2);
    assert_hits_2d_lie_on(&r, &a, &b);
}

#[test]
fn test_2d_tangent_line_has_multiplicity_two() {
    let line = Line2D::new(Point2D::new(-5.0, 2.0), Vector2D::X).unwrap();
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let r = line.intersect(&circle, TOL);
    assert_eq!(r.hits().len(), 1);
    assert_eq!(r.hits()[0].multiplicity, 2);
}

#[test]
fn test_2d_lines_cross_and_parallels_are_empty() {
    let a = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let b = Line2D::new(Point2D::new(1.0, -1.0), Vector2D::Y).unwrap();
    let r = a.intersect(&b, TOL);
    assert_eq!(r.hits().len(), 1);
    assert!((r.hits()[0].first_parameter - 1.0).abs() < 1e-12);
    assert!((r.hits()[0].second_parameter - 1.0).abs() < 1e-12);
    let parallel = Line2D::new(Point2D::new(0.0, 1.0), Vector2D::X).unwrap();
    let r = a.intersect(&parallel, TOL);
    assert!(r.is_empty() && r.is_complete());
}

#[test]
fn test_2d_routed_hits_equal_exact_hits() {
    let line = Line2D::new(Point2D::new(-5.0, 0.5), Vector2D::X).unwrap();
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let routed = Curve2D::from(line).intersect(&Curve2D::from(circle), TOL);
    let LineCircle2DIntersection::Points((t1, p1), (t2, p2)) = line.intersect_exact(&circle, TOL)
    else {
        panic!("expected two hits");
    };
    assert_eq!(routed.hits()[0].first_parameter, t1);
    assert_eq!(routed.hits()[0].point, p1);
    assert_eq!(routed.hits()[1].first_parameter, t2);
    assert_eq!(routed.hits()[1].point, p2);
}

#[test]
fn test_2d_trimmed_arc_keeps_only_its_hits() {
    let line = Line2D::new(Point2D::new(-5.0, 0.0), Vector2D::X).unwrap();
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    // The arc across the seam keeps the hit at angle 0, not the one at PI.
    let arc = TrimmedCurve2D::new(circle, -1.0, 1.0).unwrap();
    let r = arc.intersect(&line, TOL);
    assert_eq!(r.hits().len(), 1);
    assert!(r.hits()[0].first_parameter.abs() < 1e-12);
    assert!((r.hits()[0].second_parameter - 7.0).abs() < 1e-12);
}

#[test]
fn test_2d_coincident_segments_with_opposite_directions() {
    let a = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let b = Line2D::new(Point2D::new(5.0, 0.0), -Vector2D::X).unwrap();
    let first = TrimmedCurve2D::new(a, 0.0, 4.0).unwrap();
    let second = TrimmedCurve2D::new(b, 0.0, 3.0).unwrap();
    let r = first.intersect(&second, TOL);
    assert_eq!(r.overlaps().len(), 1);
    assert!(close(r.overlaps()[0].first_range, (2.0, 4.0)));
    assert!(close(r.overlaps()[0].second_range, (1.0, 3.0)));
    assert!(!r.overlaps()[0].same_direction);
    let whole = a.intersect(&a, TOL);
    assert_eq!(whole.overlaps().len(), 1);
    assert_eq!(
        whole.overlaps()[0].first_range,
        (f64::NEG_INFINITY, f64::INFINITY)
    );
    assert!(whole.overlaps()[0].same_direction);
}

#[test]
fn test_2d_coincident_arcs_overlap_on_both_sides_of_the_seam() {
    let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    let a = TrimmedCurve2D::new(circle, 0.0, 5.0).unwrap();
    let b = TrimmedCurve2D::new(circle, 4.0, 7.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert_eq!(r.overlaps().len(), 2);
    let o = r.overlaps();
    assert!(close(o[0].first_range, (0.0, 7.0 - TAU)));
    assert!(close(o[0].second_range, (TAU, 7.0)));
    assert!(close(o[1].first_range, (4.0, 5.0)));
    assert!(close(o[1].second_range, (4.0, 5.0)));
    assert!(o.iter().all(|piece| piece.same_direction));
}

#[test]
fn test_2d_coincident_circles_in_opposite_directions() {
    let a = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let frame = Frame2D::new(Point2D::ORIGIN, Vector2D::X, Vector2D::new(0.0, -1.0)).unwrap();
    let b = Circle2D::from_frame(frame, 2.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert_eq!(r.overlaps().len(), 1);
    let o = r.overlaps()[0];
    assert!(!o.same_direction);
    assert!(close(o.first_range, a.bounds()));
    assert!((o.second_range.1 - o.second_range.0 - TAU).abs() < 1e-9);
}

#[test]
fn test_consecutive_segments_meet_at_their_shared_end() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let a = TrimmedCurve::new(line, 0.0, 1.0).unwrap();
    let b = TrimmedCurve::new(line, 1.0, 2.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 1);
    let hit = r.hits()[0];
    assert!((hit.first_parameter - 1.0).abs() < 1e-12);
    assert!((hit.second_parameter - 1.0).abs() < 1e-12);
    assert_eq!(hit.multiplicity, 1);
    assert!(hit.point.distance(Point::new(1.0, 0.0, 0.0)) < 1e-12);
}

#[test]
fn test_consecutive_arcs_meet_at_their_shared_end() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let a = TrimmedCurve::new(circle, 0.0, PI / 2.0).unwrap();
    let b = TrimmedCurve::new(circle, PI / 2.0, PI).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 1);
    let hit = r.hits()[0];
    assert!((hit.first_parameter - PI / 2.0).abs() < 1e-12);
    assert!((hit.second_parameter - PI / 2.0).abs() < 1e-12);
    assert!(hit.point.distance(Point::new(0.0, 1.0, 0.0)) < 1e-12);
}

#[test]
fn test_consecutive_arcs_meet_at_the_seam() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let a = TrimmedCurve::new(circle, 1.5 * PI, TAU).unwrap();
    let b = TrimmedCurve::new(circle, 0.0, PI / 2.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 1);
    let hit = r.hits()[0];
    assert!((hit.first_parameter - TAU).abs() < 1e-12);
    assert!(hit.second_parameter.abs() < 1e-12);
    assert!(hit.point.distance(Point::new(1.0, 0.0, 0.0)) < 1e-12);
    // Read from the other side, the seam is the start of the receiver.
    let reverse = b.intersect(&a, TOL);
    assert_eq!(reverse.hits().len(), 1);
    assert!(reverse.hits()[0].first_parameter.abs() < 1e-12);
    assert!((reverse.hits()[0].second_parameter - TAU).abs() < 1e-12);
}

#[test]
fn test_arcs_overlapping_on_a_span_report_no_extra_hit() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let a = TrimmedCurve::new(circle, 0.0, PI).unwrap();
    let b = TrimmedCurve::new(circle, PI / 2.0, 1.5 * PI).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.hits().is_empty());
    assert_eq!(r.overlaps().len(), 1);
    assert!(close(r.overlaps()[0].first_range, (PI / 2.0, PI)));
    // An arc against its full circle overlaps it without a hit at the seam.
    let arc = TrimmedCurve::new(circle, 0.0, 1.0).unwrap();
    let r = arc.intersect(&circle, TOL);
    assert!(r.hits().is_empty());
    assert_eq!(r.overlaps().len(), 1);
}

#[test]
fn test_2d_consecutive_segments_meet_at_their_shared_end() {
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let a = TrimmedCurve2D::new(line, 0.0, 1.0).unwrap();
    let b = TrimmedCurve2D::new(line, 1.0, 2.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 1);
    let hit = r.hits()[0];
    assert!((hit.first_parameter - 1.0).abs() < 1e-12);
    assert!((hit.second_parameter - 1.0).abs() < 1e-12);
    assert_eq!(hit.multiplicity, 1);
    assert!(hit.point.distance(Point2D::new(1.0, 0.0)) < 1e-12);
}

#[test]
fn test_2d_consecutive_arcs_meet_at_their_shared_end() {
    let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    let a = TrimmedCurve2D::new(circle, 0.0, PI / 2.0).unwrap();
    let b = TrimmedCurve2D::new(circle, PI / 2.0, PI).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 1);
    let hit = r.hits()[0];
    assert!((hit.first_parameter - PI / 2.0).abs() < 1e-12);
    assert!((hit.second_parameter - PI / 2.0).abs() < 1e-12);
    assert!(hit.point.distance(Point2D::new(0.0, 1.0)) < 1e-12);
}

#[test]
fn test_2d_arcs_overlapping_on_a_span_report_no_extra_hit() {
    let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    let a = TrimmedCurve2D::new(circle, 0.0, PI).unwrap();
    let b = TrimmedCurve2D::new(circle, PI / 2.0, 1.5 * PI).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.hits().is_empty());
    assert_eq!(r.overlaps().len(), 1);
}

#[test]
fn test_hit_just_below_the_seam_has_its_angle_at_the_seam() {
    // The line crosses the circle a hair below angle 0: the circle's angle
    // is recovered from the point, near the seam.
    let near_seam = |angle: f64| angle.min(TAU - angle) < 1e-9;
    let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    let line = Line2D::new(Point2D::new(-5.0, -1e-11), Vector2D::X).unwrap();
    let r = circle.intersect(&line, TOL);
    assert_eq!(r.hits().len(), 2);
    let seam_hit = r.hits().iter().find(|h| h.point.x > 0.0).unwrap();
    assert!(near_seam(seam_hit.first_parameter));
    let circle_3d = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let line_3d = Line::new(Point::new(-5.0, -1e-11, 0.0), Vector::X).unwrap();
    let r_3d = circle_3d.intersect(&line_3d, TOL);
    let seam_hit = r_3d.hits().iter().find(|h| h.point.x > 0.0).unwrap();
    assert!(near_seam(seam_hit.first_parameter));
    // Two circles meeting a hair below the seam of the first.
    let other = Circle2D::new(Point2D::new(2.0, -2e-11), 1.0).unwrap();
    let r = circle.intersect(&other, TOL);
    assert_eq!(r.hits().len(), 1);
    assert!(near_seam(r.hits()[0].first_parameter));
}

#[test]
fn test_concrete_pairs_route_like_their_enums_bit_for_bit() {
    let line = Line::new(Point::new(-5.0, 0.3, 0.2), Vector::new(1.0, 0.1, 0.4)).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::new(0.1, 0.0, 1.0), 1.5).unwrap();
    let typed = line.intersect(&cylinder, TOL);
    assert_eq!(typed.hits().len(), 2);
    assert_eq!(
        typed,
        Curve::from(line).intersect(&Surface::from(cylinder), TOL)
    );
    let segment = TrimmedCurve::new(line, 0.0, 5.0).unwrap();
    assert_eq!(
        segment.intersect(&cylinder, TOL),
        TrimmedCurve::new(Curve::from(line), 0.0, 5.0)
            .unwrap()
            .intersect(&Surface::from(cylinder), TOL)
    );

    let first = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let second = Circle::new(Point::new(1.5, 0.0, 0.0), Vector::Z, 1.0).unwrap();
    let typed = first.intersect(&second, TOL);
    assert_eq!(typed.hits().len(), 2);
    assert_eq!(
        typed,
        Curve::from(first).intersect(&Curve::from(second), TOL)
    );

    let first = Circle2D::new(Point2D::new(0.0, 0.0), 2.0).unwrap();
    let second = Circle2D::new(Point2D::new(1.5, 0.5), 1.0).unwrap();
    let typed = first.intersect(&second, TOL);
    assert_eq!(typed.hits().len(), 2);
    assert_eq!(
        typed,
        Curve2D::from(first).intersect(&Curve2D::from(second), TOL)
    );
}

/// A curve type defined outside the crate: a helix of radius 2 about the
/// z axis, rising 1 per turn, over three turns.
struct Helix;

impl ParametricCurve for Helix {
    fn point_at(&self, u: f64) -> Point {
        Point::new(2.0 * u.cos(), 2.0 * u.sin(), u / TAU)
    }
    fn derivative(&self, u: f64, order: u32) -> Vector {
        match order % 4 {
            1 => Vector::new(-2.0 * u.sin(), 2.0 * u.cos(), 1.0 / TAU),
            2 => Vector::new(-2.0 * u.cos(), -2.0 * u.sin(), 0.0),
            3 => Vector::new(2.0 * u.sin(), -2.0 * u.cos(), 0.0),
            _ => Vector::new(2.0 * u.cos(), 2.0 * u.sin(), 0.0),
        }
    }
    fn bounds(&self) -> (f64, f64) {
        (0.0, 3.0 * TAU)
    }
    fn period(&self) -> Option<f64> {
        None
    }
}

#[test]
fn test_caller_defined_curve_intersects_through_the_public_numeric_functions() {
    // A slanted plane through the axis at height 1.5 crosses each turn.
    let plane = Plane::new(Point::new(0.0, 0.0, 1.5), Vector::new(0.05, 0.0, 1.0)).unwrap();
    let r = intersect_curve_surface(&Helix, &plane, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 1);
    assert_hits_lie_on(&r, &Helix, &plane);

    // A vertical line on the helix's cylinder meets it once per turn.
    let line = Line::new(Point::new(0.0, 2.0, 0.0), Vector::Z).unwrap();
    let r = intersect_curve_curve(&Helix, &line, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    let parameters: Vec<f64> = r.hits().iter().map(|h| h.first_parameter).collect();
    assert_eq!(parameters.len(), 3);
    for (turn, parameter) in parameters.iter().enumerate() {
        assert!((parameter - (PI / 2.0 + turn as f64 * TAU)).abs() < 1e-9);
    }
}

#[test]
fn test_segments_apart_by_less_than_the_confusion_touch() {
    // A unit-speed line: a parameter gap of 5e-9 is 5e-9 in space.
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let a = TrimmedCurve::new(line, 0.0, 1.0).unwrap();
    let b = TrimmedCurve::new(line, 1.0 + 5e-9, 2.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 1);
    assert_eq!(r.hits()[0].multiplicity, 1);
    assert!(r.hits()[0].point.distance(Point::new(1.0, 0.0, 0.0)) < 1e-8);
}

#[test]
fn test_segments_apart_by_more_than_the_confusion_do_not_meet() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let a = TrimmedCurve::new(line, 0.0, 1.0).unwrap();
    let b = TrimmedCurve::new(line, 1.0 + 1e-3, 2.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.is_empty());
    assert!(r.is_complete());
}

#[test]
fn test_2d_consecutive_arcs_meet_at_the_seam() {
    let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    let a = TrimmedCurve2D::new(circle, 1.5 * PI, TAU).unwrap();
    let b = TrimmedCurve2D::new(circle, 0.0, PI / 2.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 1);
    let hit = r.hits()[0];
    assert!((hit.first_parameter - TAU).abs() < 1e-12);
    assert!(hit.second_parameter.abs() < 1e-12);
    assert!(hit.point.distance(Point2D::new(1.0, 0.0)) < 1e-12);
}

#[test]
fn test_touch_between_coincident_curves_running_opposite_ways() {
    // The reversed circle has the point at angle a at parameter -a.
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let reversed = Circle::from_frame(
        Frame::new(Point::ORIGIN, -Vector::Z, Vector::X).unwrap(),
        1.0,
    )
    .unwrap();
    let a = TrimmedCurve::new(circle, 0.0, PI / 2.0).unwrap();
    // Angles [PI / 2, PI] on the reversed circle.
    let b = TrimmedCurve::new(reversed, PI, 1.5 * PI).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 1);
    let hit = r.hits()[0];
    assert!((hit.first_parameter - PI / 2.0).abs() < 1e-12);
    assert!((hit.second_parameter - 1.5 * PI).abs() < 1e-12);
    assert!(hit.point.distance(Point::new(0.0, 1.0, 0.0)) < 1e-12);

    // A segment of the line running the other way: x in [1, 2].
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let backwards = Line::new(Point::ORIGIN, -Vector::X).unwrap();
    let a = TrimmedCurve::new(line, 0.0, 1.0).unwrap();
    let b = TrimmedCurve::new(backwards, -2.0, -1.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 1);
    assert!((r.hits()[0].first_parameter - 1.0).abs() < 1e-12);
    assert!((r.hits()[0].second_parameter + 1.0).abs() < 1e-12);
}

#[test]
fn test_arcs_making_up_the_whole_circle_touch_twice() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let a = TrimmedCurve::new(circle, 0.0, PI).unwrap();
    let b = TrimmedCurve::new(circle, PI, TAU).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.overlaps().is_empty());
    let pairs: Vec<(f64, f64)> = r
        .hits()
        .iter()
        .map(|h| (h.first_parameter, h.second_parameter))
        .collect();
    assert_eq!(pairs.len(), 2);
    assert!(close(pairs[0], (0.0, TAU)));
    assert!(close(pairs[1], (PI, PI)));

    let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    let a = TrimmedCurve2D::new(circle, 0.0, PI).unwrap();
    let b = TrimmedCurve2D::new(circle, PI, TAU).unwrap();
    let r = a.intersect(&b, TOL);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 2);
}

#[test]
fn test_routed_circle_circle_recovers_both_angles() {
    let first = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let second = Circle::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    let r = first.intersect(&second, TOL);
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert_eq!(r.hits().len(), 2);
    for h in r.hits() {
        assert!(first.point_at(h.first_parameter).distance(h.point) < 1e-12);
        assert!(second.point_at(h.second_parameter).distance(h.point) < 1e-12);
    }
}

#[test]
fn test_trimmed_operands_store_their_parameter_within_range() {
    let first = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let second = Circle::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    // An arc across the seam of the first circle holds both hits.
    let arc = TrimmedCurve::new(first, -PI / 2.0, PI / 2.0).unwrap();
    let r = arc.intersect(&second, TOL);
    assert_eq!(r.hits().len(), 2);
    for h in r.hits() {
        let angle = h.first_parameter;
        assert!((-PI / 2.0..=PI / 2.0).contains(&angle));
        // The stored angle is the circle's own, up to a whole turn.
        let own = first.parameter_of(h.point, Tolerance::DEFAULT);
        assert!(((angle - own).rem_euclid(TAU) + 1e-12) % TAU < 1e-9);
        // The untrimmed side has its own angle.
        assert!(second.point_at(h.second_parameter).distance(h.point) < 1e-12);
    }
    // The hits are sorted by their first parameter.
    assert!(r.hits()[0].first_parameter < r.hits()[1].first_parameter);
}

#[test]
fn test_numeric_path_stores_every_parameter() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let other = Circle::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    let r = circle.intersect(&other.to_bspline(), TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert_eq!(r.hits().len(), 2);
    assert_curve_hits_lie_on(&r, &circle, &other.to_bspline());
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let ellipse = Ellipse::from_frame(Frame::WORLD, 5.0, 1.0).unwrap();
    let r = ellipse.intersect(&torus, TOL);
    assert!(!r.hits().is_empty());
    assert_hits_lie_on(&r, &ellipse, &torus);
}

/// A hit within `tol.confusion` in space of the seam, just outside an arc
/// that starts or ends there, on a unit circle in the xy plane.
struct SeamCase {
    /// The arc's range.
    range: (f64, f64),
    /// The angle of the hit: a hair past the arc's seam end.
    angle: f64,
    /// The direction angle of the line through the hit; its second hit on
    /// the circle falls outside the arc.
    line_angle: f64,
    /// The angle of the centre of a unit circle through the hit; its
    /// second hit falls outside the arc.
    crossing_angle: f64,
    /// The arc parameter the hit keeps: the seam's representative.
    expected: f64,
}

const SEAM_CASES: [SeamCase; 2] = [
    SeamCase {
        range: (0.0, PI),
        angle: -5e-9,
        line_angle: PI / 3.0,
        crossing_angle: -PI / 3.0,
        expected: 0.0,
    },
    SeamCase {
        range: (PI, TAU),
        angle: 5e-9,
        line_angle: 2.0 * PI / 3.0,
        crossing_angle: PI / 3.0,
        expected: TAU,
    },
];

fn assert_seam_parameter(parameter: f64, case: &SeamCase) {
    assert!(
        (parameter - case.expected).abs() < 1e-12,
        "parameter {parameter}, expected {}",
        case.expected
    );
}

#[test]
fn test_hit_a_hair_past_the_seam_end_of_an_arc_is_kept_in_3d() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    for case in &SEAM_CASES {
        let arc = TrimmedCurve::new(circle, case.range.0, case.range.1).unwrap();
        let hit = circle.point_at(case.angle);
        let (sin, cos) = case.line_angle.sin_cos();
        let line = Line::new(hit, Vector::new(cos, sin, 0.0)).unwrap();
        let (sin, cos) = (case.crossing_angle + case.angle).sin_cos();
        let crossing = Circle::new(Point::new(cos, sin, 0.0), Vector::Z, 1.0).unwrap();
        let (sin, cos) = case.line_angle.sin_cos();
        let plane = Plane::new(hit, Vector::new(-sin, cos, 0.0)).unwrap();

        let r = arc.intersect(&line, TOL);
        assert_eq!(r.hits().len(), 1, "arc with line, range {:?}", case.range);
        assert_seam_parameter(r.hits()[0].first_parameter, case);
        let r = line.intersect(&arc, TOL);
        assert_eq!(r.hits().len(), 1, "line with arc, range {:?}", case.range);
        assert_seam_parameter(r.hits()[0].second_parameter, case);
        let r = arc.intersect(&crossing, TOL);
        assert_eq!(r.hits().len(), 1, "arc with circle, range {:?}", case.range);
        assert_seam_parameter(r.hits()[0].first_parameter, case);
        let r = crossing.intersect(&arc, TOL);
        assert_eq!(r.hits().len(), 1, "circle with arc, range {:?}", case.range);
        assert_seam_parameter(r.hits()[0].second_parameter, case);
        let r = arc.intersect(&plane, TOL);
        assert_eq!(r.hits().len(), 1, "arc with plane, range {:?}", case.range);
        assert_seam_parameter(r.hits()[0].curve_parameter, case);
    }
}

#[test]
fn test_hit_a_hair_past_the_seam_end_of_an_arc_is_kept_in_2d() {
    let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    for case in &SEAM_CASES {
        let arc = TrimmedCurve2D::new(circle, case.range.0, case.range.1).unwrap();
        let hit = circle.point_at(case.angle);
        let (sin, cos) = case.line_angle.sin_cos();
        let line = Line2D::new(hit, Vector2D::new(cos, sin)).unwrap();
        let (sin, cos) = (case.crossing_angle + case.angle).sin_cos();
        let crossing = Circle2D::new(Point2D::new(cos, sin), 1.0).unwrap();

        let r = arc.intersect(&line, TOL);
        assert_eq!(r.hits().len(), 1, "arc with line, range {:?}", case.range);
        assert_seam_parameter(r.hits()[0].first_parameter, case);
        let r = line.intersect(&arc, TOL);
        assert_eq!(r.hits().len(), 1, "line with arc, range {:?}", case.range);
        assert_seam_parameter(r.hits()[0].second_parameter, case);
        let r = arc.intersect(&crossing, TOL);
        assert_eq!(r.hits().len(), 1, "arc with circle, range {:?}", case.range);
        assert_seam_parameter(r.hits()[0].first_parameter, case);
        let r = crossing.intersect(&arc, TOL);
        assert_eq!(r.hits().len(), 1, "circle with arc, range {:?}", case.range);
        assert_seam_parameter(r.hits()[0].second_parameter, case);
    }
}

#[test]
fn test_hit_farther_than_the_confusion_past_the_seam_is_dropped() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let arc = TrimmedCurve::new(circle, 0.0, PI).unwrap();
    let hit = circle.point_at(-1e-6);
    let (sin, cos) = (PI / 3.0).sin_cos();
    let line = Line::new(hit, Vector::new(cos, sin, 0.0)).unwrap();
    assert!(arc.intersect(&line, TOL).hits().is_empty());
}

/// Two arcs of one circle whose ends meet at the seam, `gap` apart in
/// angle: `[gap, PI/2]` and `[3PI/2, 2PI]`, then the second shifted one
/// turn up against `[PI, 2PI]`.
fn seam_touch_ranges(gap: f64) -> [((f64, f64), (f64, f64)); 2] {
    [
        ((gap, PI / 2.0), (3.0 * PI / 2.0, TAU)),
        ((TAU + gap, TAU + 1.0), (PI, TAU)),
    ]
}

fn assert_touch_parameters(
    first: f64,
    second: f64,
    eval: (Point, Point),
    point: Point,
    ranges: ((f64, f64), (f64, f64)),
) {
    let within = |p: f64, (lo, hi): (f64, f64)| p >= lo - 1e-9 && p <= hi + 1e-9;
    assert!(
        within(first, ranges.0),
        "first {first} outside {:?}",
        ranges.0
    );
    assert!(
        within(second, ranges.1),
        "second {second} outside {:?}",
        ranges.1
    );
    assert!(
        eval.0.distance(point) <= 1e-7,
        "first parameter is not the point"
    );
    assert!(
        eval.1.distance(point) <= 1e-7,
        "second parameter is not the point"
    );
}

#[test]
fn test_arcs_touching_across_the_seam_store_the_near_ends() {
    for radius in [1.0, 0.5, 0.01] {
        for gap_in_space in [0.0, 1e-10, 5e-9, 5e-8, 9e-8] {
            let gap = gap_in_space / radius;
            for (range_a, range_b) in seam_touch_ranges(gap) {
                let circle = Circle::new(Point::ORIGIN, Vector::Z, radius).unwrap();
                let a = TrimmedCurve::new(circle, range_a.0, range_a.1).unwrap();
                let b = TrimmedCurve::new(circle, range_b.0, range_b.1).unwrap();
                for (first, second, ranges) in
                    [(&a, &b, (range_a, range_b)), (&b, &a, (range_b, range_a))]
                {
                    let r = first.intersect(second, TOL);
                    assert!(r.overlaps().is_empty());
                    assert_eq!(r.hits().len(), 1, "radius {radius}, gap {gap}");
                    let h = r.hits()[0];
                    assert_touch_parameters(
                        h.first_parameter,
                        h.second_parameter,
                        (
                            circle.point_at(h.first_parameter),
                            circle.point_at(h.second_parameter),
                        ),
                        h.point,
                        ranges,
                    );
                }
                let circle = Circle2D::new(Point2D::ORIGIN, radius).unwrap();
                let a = TrimmedCurve2D::new(circle, range_a.0, range_a.1).unwrap();
                let b = TrimmedCurve2D::new(circle, range_b.0, range_b.1).unwrap();
                let lift = |p: Point2D| Point::new(p.x, p.y, 0.0);
                for (first, second, ranges) in
                    [(&a, &b, (range_a, range_b)), (&b, &a, (range_b, range_a))]
                {
                    let r = first.intersect(second, TOL);
                    assert!(r.overlaps().is_empty());
                    assert_eq!(r.hits().len(), 1, "2D radius {radius}, gap {gap}");
                    let h = r.hits()[0];
                    assert_touch_parameters(
                        h.first_parameter,
                        h.second_parameter,
                        (
                            lift(circle.point_at(h.first_parameter)),
                            lift(circle.point_at(h.second_parameter)),
                        ),
                        lift(h.point),
                        ranges,
                    );
                }
            }
        }
    }
}

#[test]
fn test_reviewer_seam_touch_stores_two_pi_on_the_second_arc() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let a = TrimmedCurve::new(circle, 5e-8, PI / 2.0).unwrap();
    let b = TrimmedCurve::new(circle, 3.0 * PI / 2.0, TAU).unwrap();
    let r = a.intersect(&b, TOL);
    assert_eq!(r.hits().len(), 1);
    assert!((r.hits()[0].second_parameter - TAU).abs() < 1e-12);
}

/// A caller-defined surface: a plane the crate does not know.
struct Floor(Plane);

impl ParametricSurface for Floor {
    fn point_at(&self, u: f64, v: f64) -> Point {
        self.0.point_at(u, v)
    }
    fn derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        self.0.derivative(u, v, du, dv)
    }
    fn u_bounds(&self) -> (f64, f64) {
        self.0.u_bounds()
    }
    fn v_bounds(&self) -> (f64, f64) {
        self.0.v_bounds()
    }
    fn u_period(&self) -> Option<f64> {
        None
    }
    fn v_period(&self) -> Option<f64> {
        None
    }
    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        self.0.project_point(point, tol)
    }
    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        self.0.contains(point, tol)
    }
}

#[test]
fn test_caller_defined_surface_intersects_through_the_public_numeric_function() {
    let floor = Floor(Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap());
    let circle = Circle::new(Point::ORIGIN, Vector::X, 2.0).unwrap();
    let r = intersect_curve_surface(&circle, &floor, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert!(r.overlaps().is_empty());
    assert_eq!(r.hits().len(), 2);
    assert_hits_lie_on(&r, &circle, &floor);
    for h in r.hits() {
        assert!((h.point.z - 1.0).abs() < 1e-9);
    }
}

/// Hits of the arc `[0, PI]` of a circle of `radius` in the xy plane with
/// a line, a plane and (2D) a line, all through the point at angle
/// `PI + excess`, a hair past the arc's end; their other hit on the circle
/// falls outside the arc. Returns the stored arc parameters, `None` for a
/// dropped hit.
fn hits_past_the_end_of_an_arc(radius: f64, excess: f64) -> Vec<Option<f64>> {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, radius).unwrap();
    let arc = TrimmedCurve::new(circle, 0.0, PI).unwrap();
    let hit = circle.point_at(PI + excess);
    let (sin, cos) = (-PI / 3.0).sin_cos();
    let line = Line::new(hit, Vector::new(cos, sin, 0.0)).unwrap();
    let plane = Plane::new(hit, Vector::new(-sin, cos, 0.0)).unwrap();
    let circle_2d = Circle2D::new(Point2D::ORIGIN, radius).unwrap();
    let arc_2d = TrimmedCurve2D::new(circle_2d, 0.0, PI).unwrap();
    let line_2d = Line2D::new(Point2D::new(hit.x, hit.y), Vector2D::new(cos, sin)).unwrap();
    let only = |parameters: Vec<f64>| {
        assert!(parameters.len() <= 1, "radius {radius}, excess {excess}");
        parameters.first().copied()
    };
    vec![
        only(
            arc.intersect(&line, TOL)
                .hits()
                .iter()
                .map(|h| h.first_parameter)
                .collect(),
        ),
        only(
            line.intersect(&arc, TOL)
                .hits()
                .iter()
                .map(|h| h.second_parameter)
                .collect(),
        ),
        only(
            arc.intersect(&plane, TOL)
                .hits()
                .iter()
                .map(|h| h.curve_parameter)
                .collect(),
        ),
        only(
            arc_2d
                .intersect(&line_2d, TOL)
                .hits()
                .iter()
                .map(|h| h.first_parameter)
                .collect(),
        ),
        only(
            line_2d
                .intersect(&arc_2d, TOL)
                .hits()
                .iter()
                .map(|h| h.second_parameter)
                .collect(),
        ),
    ]
}

#[test]
fn test_hit_a_hair_past_a_plain_end_of_an_arc_is_kept_at_the_end() {
    for parameter in hits_past_the_end_of_an_arc(1.0, 5e-9) {
        assert_eq!(parameter, Some(PI));
    }
    for parameter in hits_past_the_end_of_an_arc(1.0, 1e-6) {
        assert_eq!(parameter, None);
    }
}

#[test]
fn test_excess_past_a_range_end_is_measured_in_space() {
    // 5e-8 in space is kept, 2e-7 is not, whatever the radius.
    for radius in [10.0, 0.01] {
        for parameter in hits_past_the_end_of_an_arc(radius, 5e-8 / radius) {
            assert_eq!(parameter, Some(PI), "radius {radius}");
        }
        for parameter in hits_past_the_end_of_an_arc(radius, 2e-7 / radius) {
            assert_eq!(parameter, None, "radius {radius}");
        }
    }
    // A segment: 5e-8 past its end in space (fifty times `tol.parametric`)
    // is kept at the end; 2e-7 is not.
    let segment =
        TrimmedCurve::new(Line::new(Point::ORIGIN, Vector::X).unwrap(), 0.0, 1.0).unwrap();
    for (excess, expected) in [(5e-8, Some(1.0)), (2e-7, None)] {
        let crossing = Line::new(Point::new(1.0 + excess, -1.0, 0.0), Vector::Y).unwrap();
        let r = segment.intersect(&crossing, TOL);
        assert_eq!(r.hits().first().map(|h| h.first_parameter), expected);
        let r = crossing.intersect(&segment, TOL);
        assert_eq!(r.hits().first().map(|h| h.second_parameter), expected);
        let wall = Plane::new(Point::new(1.0 + excess, 0.0, 0.0), Vector::X).unwrap();
        let r = segment.intersect(&wall, TOL);
        assert_eq!(r.hits().first().map(|h| h.curve_parameter), expected);
    }
}

// ---- 2D B-splines: the numeric search, embedded at z = 0 ----

/// Every hit's parameters give its point on both curves.
fn assert_hits_lie_on_2d<A: ParametricCurve2D, B: ParametricCurve2D>(
    r: &CurveCurveIntersection2D,
    a: &A,
    b: &B,
) {
    for h in r.hits() {
        let (s, t) = (h.first_parameter, h.second_parameter);
        assert!(a.point_at(s).distance(h.point) <= 1e-7);
        assert!(b.point_at(t).distance(h.point) <= 1e-7);
    }
}

/// Asserts that `r`'s hits are at `expected`, in any order.
#[track_caller]
fn assert_hit_points_2d(r: &CurveCurveIntersection2D, expected: &[Point2D]) {
    let actual: Vec<Point2D> = r.hits().iter().map(|h| h.point).collect();
    assert_eq!(actual.len(), expected.len(), "{actual:?}");
    for q in expected {
        assert!(
            actual.iter().any(|p| p.distance(*q) <= 1e-9),
            "{q:?} not among {actual:?}"
        );
    }
}

#[test]
fn test_bspline_2d_against_a_line_is_numeric_and_complete() {
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap().to_bspline();
    let line = Line2D::new(Point2D::new(-5.0, 0.0), Vector2D::X).unwrap();
    for r in [circle.intersect(&line, TOL), line.intersect(&circle, TOL)] {
        assert_eq!(r.method(), IntersectionMethod::Numeric);
        assert!(r.is_complete());
        assert_hit_points_2d(&r, &[Point2D::new(-2.0, 0.0), Point2D::new(2.0, 0.0)]);
    }
    assert_hits_lie_on_2d(&circle.intersect(&line, TOL), &circle, &line);
}

#[test]
fn test_two_bspline_2d_circles_meet_twice() {
    let a = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap().to_bspline();
    let b = Circle2D::new(Point2D::new(1.0, 0.0), 1.0)
        .unwrap()
        .to_bspline();
    let r = a.intersect(&b, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert!(r.is_complete());
    let h = 0.75f64.sqrt();
    assert_hit_points_2d(&r, &[Point2D::new(0.5, -h), Point2D::new(0.5, h)]);
    assert_hits_lie_on_2d(&r, &a, &b);
}

#[test]
fn test_bspline_2d_variant_against_a_circle() {
    let spline: Curve2D = Circle2D::new(Point2D::ORIGIN, 1.0)
        .unwrap()
        .to_bspline()
        .into();
    let circle = Circle2D::new(Point2D::new(1.0, 0.0), 1.0).unwrap();
    let r = spline.intersect(&circle, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert!(r.is_complete());
    assert_eq!(r.hits().len(), 2);
    assert_hits_lie_on_2d(&r, &spline, &circle);
    // Two elementary variants keep the closed form.
    let exact: Curve2D = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap().into();
    assert_eq!(
        exact.intersect(&circle, TOL).method(),
        IntersectionMethod::Exact
    );
}

#[test]
fn test_trimmed_bspline_2d_keeps_the_hits_in_its_range() {
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap().to_bspline();
    let line = Line2D::new(Point2D::new(-5.0, 0.0), Vector2D::X).unwrap();
    // The half of the curve through (2, 0) only.
    let start = circle.project_point(Point2D::new(0.0, -2.0), TOL).parameter;
    let end = circle.project_point(Point2D::new(0.0, 2.0), TOL).parameter;
    let (first, last) = if start < end {
        (start, end)
    } else {
        (start, end + circle.bounds().1 - circle.bounds().0)
    };
    let half = TrimmedCurve2D::new(circle.clone(), first, last).unwrap();
    let r = half.intersect(&line, TOL);
    assert_eq!(r.method(), IntersectionMethod::Numeric);
    assert_hit_points_2d(&r, &[Point2D::new(2.0, 0.0)]);
    assert_hits_lie_on_2d(&r, &half, &line);
}

// ---- 2D conics: every pair of lines and conics is exact ----

fn exact_hits_2d(result: Conic2DIntersection) -> Vec<CurveCurveHit2D> {
    match result {
        Conic2DIntersection::Hits(hits) => hits.to_vec(),
        other => panic!("expected hits, got {other:?}"),
    }
}

#[test]
fn test_2d_conic_pairs_are_exact_through_every_type() {
    let line = Line2D::new(Point2D::new(0.0, 1.0), Vector2D::X).unwrap();
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 2.0).unwrap();
    let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0, 0.5).unwrap();
    let all: [Curve2D; 5] = [
        line.into(),
        circle.into(),
        ellipse.into(),
        parabola.into(),
        hyperbola.into(),
    ];
    for a in &all {
        for b in &all {
            let r = a.intersect(b, TOL);
            assert_eq!(r.method(), IntersectionMethod::Exact);
            assert!(r.is_complete());
            // The reverse order meets at the same points.
            let back = b.intersect(a, TOL);
            assert_eq!(back.method(), IntersectionMethod::Exact);
            assert_eq!(back.hits().len(), r.hits().len());
            for hit in r.hits() {
                assert!(
                    back.hits()
                        .iter()
                        .any(|other| other.point.distance(hit.point) <= 1e-9
                            && other.multiplicity == hit.multiplicity)
                );
            }
        }
    }
}

#[test]
fn test_2d_typed_conics_route_exact_and_match_the_typed_solver() {
    let line = Line2D::new(Point2D::new(0.0, 1.0), Vector2D::X).unwrap();
    let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 2.0).unwrap();
    let routed = line.intersect(&ellipse, TOL);
    let typed = exact_hits_2d(line.intersect_exact(&ellipse, TOL));
    assert_eq!(routed.method(), IntersectionMethod::Exact);
    assert_eq!(routed.hits(), typed.as_slice());
    assert_hits_2d_lie_on(&routed, &line, &ellipse);
    // The reverse order swaps which parameter slot holds the line vs conic.
    let reverse = exact_hits_2d(ellipse.intersect_exact(&line, TOL));
    assert_eq!(typed.len(), reverse.len());
    for a in &typed {
        let b = reverse
            .iter()
            .find(|other| other.point.distance(a.point) < 1e-12)
            .expect("reverse order reports the same point");
        assert_eq!(a.first_parameter, b.second_parameter);
        assert_eq!(a.second_parameter, b.first_parameter);
        assert_eq!(a.multiplicity, b.multiplicity);
    }
}

#[test]
fn test_2d_conic_tangency_multiplicities() {
    let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 2.0).unwrap();
    let tangent = Line2D::new(Point2D::new(0.0, 2.0), Vector2D::X).unwrap();
    let hits = exact_hits_2d(tangent.intersect_exact(&ellipse, TOL));
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].multiplicity, 2);
    // The circle of the minor radius touches the ellipse at its two minor vertices.
    let inner = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let hits = exact_hits_2d(inner.intersect_exact(&ellipse, TOL));
    assert_eq!(hits.len(), 2);
    assert!(hits.iter().all(|hit| hit.multiplicity == 2));
    let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0, 1.0).unwrap();
    let apex_line = Line2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    assert_eq!(
        exact_hits_2d(parabola.intersect_exact(&apex_line, TOL))[0].multiplicity,
        2
    );
    let vertex_line = Line2D::new(Point2D::new(1.0, 0.0), Vector2D::Y).unwrap();
    assert_eq!(
        exact_hits_2d(hyperbola.intersect_exact(&vertex_line, TOL))[0].multiplicity,
        2
    );
}

#[test]
fn test_2d_conic_coincidence_and_disjoint_typed_results() {
    let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 2.0).unwrap();
    assert_eq!(
        ellipse.intersect_exact(&ellipse, TOL),
        Conic2DIntersection::Coincident
    );
    let far = Ellipse2D::new(Point2D::new(10.0, 0.0), Vector2D::X, 3.0, 2.0).unwrap();
    assert_eq!(
        ellipse.intersect_exact(&far, TOL),
        Conic2DIntersection::Empty
    );
    let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0, 1.0).unwrap();
    assert_eq!(
        hyperbola.intersect_exact(&hyperbola, TOL),
        Conic2DIntersection::Coincident
    );
}

#[test]
fn test_2d_trimmed_conics_clip_in_both_orders() {
    let line = Line2D::new(Point2D::new(1.5, 0.0), Vector2D::Y).unwrap();
    let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 2.0).unwrap();
    let upper = TrimmedCurve2D::new(ellipse, 0.0, PI).unwrap();
    let wanted = [Point2D::new(1.5, 3.0f64.sqrt())];
    for r in [line.intersect(&upper, TOL), upper.intersect(&line, TOL)] {
        assert_eq!(r.method(), IntersectionMethod::Exact);
        assert_hit_points_2d(&r, &wanted);
    }
    let segment = TrimmedCurve2D::new(line, 0.0, 1.0).unwrap();
    // The segment stays below the hit.
    assert!(segment.intersect(&upper, TOL).is_empty());
    assert!(upper.intersect(&segment, TOL).is_empty());
    // A trimmed conic against a trimmed conic.
    let parabola = Parabola2D::new(Point2D::new(4.0, 0.0), -Vector2D::X, 1.0).unwrap();
    let arc = TrimmedCurve2D::new(parabola, 0.0, 5.0).unwrap();
    let r = upper.intersect(&arc, TOL);
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert_hits_2d_lie_on(&r, &upper, &arc);
    assert_eq!(r.hits().len(), arc.intersect(&upper, TOL).hits().len());
}

#[test]
fn test_2d_coincident_conic_pieces_overlap() {
    let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    let a = TrimmedCurve2D::new(parabola, -2.0, 2.0).unwrap();
    let b = TrimmedCurve2D::new(parabola, 0.0, 5.0).unwrap();
    let r = a.intersect(&b, TOL);
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert_eq!(r.overlaps().len(), 1);
    assert!(close(r.overlaps()[0].first_range, (0.0, 2.0)));
    assert!(close(r.overlaps()[0].second_range, (0.0, 2.0)));
    // A whole parabola against a piece of itself overlaps over the piece.
    let whole = parabola.intersect(&b, TOL);
    assert_eq!(whole.overlaps().len(), 1);
    assert!(close(whole.overlaps()[0].first_range, (0.0, 5.0)));
    let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 2.0).unwrap();
    let full = ellipse.intersect(&ellipse, TOL);
    assert_eq!(full.overlaps().len(), 1);
    assert!(close(full.overlaps()[0].first_range, (0.0, TAU)));
}

#[test]
fn test_2d_bspline_against_a_conic_is_numeric() {
    let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 2.0).unwrap();
    let spline = Ellipse2D::new(Point2D::ORIGIN, Vector2D::Y, 3.0, 2.0)
        .unwrap()
        .to_bspline();
    for r in [
        spline.intersect(&ellipse, TOL),
        ellipse.intersect(&spline, TOL),
    ] {
        assert_eq!(r.method(), IntersectionMethod::Numeric);
        assert!(r.is_complete());
        // x^2/9 + y^2/4 = 1 and x^2/4 + y^2/9 = 1: x^2 = y^2 = 36/13.
        let c = (36.0f64 / 13.0).sqrt();
        assert_hit_points_2d(
            &r,
            &[
                Point2D::new(c, c),
                Point2D::new(-c, c),
                Point2D::new(c, -c),
                Point2D::new(-c, -c),
            ],
        );
    }
    assert_hits_lie_on_2d(&spline.intersect(&ellipse, TOL), &spline, &ellipse);
}

/// The x coordinates of the hits, sorted.
fn hit_xs(r: &CurveSurfaceIntersection) -> Vec<f64> {
    let mut xs: Vec<f64> = r.hits().iter().map(|h| h.point.x).collect();
    xs.sort_by(f64::total_cmp);
    xs
}

/// A vase: an interpolated profile turned about the z axis.
fn vase() -> RevolutionSurface {
    let points = [
        Point::new(3.0, 0.0, 0.0),
        Point::new(4.0, 0.0, 3.0),
        Point::new(2.0, 0.0, 6.0),
        Point::new(2.5, 0.0, 9.0),
        Point::new(3.2, 0.0, 10.0),
    ];
    let profile = BSplineCurve::interpolate(&points, 3, PointParametrization::Centripetal).unwrap();
    RevolutionSurface::new(profile.into(), Axis::new(Point::ORIGIN, Vector::Z).unwrap()).unwrap()
}

#[test]
fn test_unbounded_line_meets_a_revolution_surface_on_both_sides_of_its_origin() {
    // #174: the line was seeded over a fixed window of its parameter
    // around its origin, so the crossing more than 10 away was missed.
    let vase = vase();
    let reference = TrimmedCurve::segment(Point::new(-30.0, 0.5, 2.0), Point::new(30.0, 0.5, 2.0))
        .unwrap()
        .intersect(&vase, TOL);
    assert_eq!(reference.hits().len(), 2);
    let expected = hit_xs(&reference);
    for x0 in [-30.0, -10.0, 0.0, 10.0, 30.0] {
        let line = Line::new(Point::new(x0, 0.5, 2.0), Vector::X).unwrap();
        for r in [line.intersect(&vase, TOL), vase.intersect(&line, TOL)] {
            let xs = hit_xs(&r);
            assert_eq!(xs.len(), 2, "origin x = {x0}: {xs:?}");
            for (x, e) in xs.iter().zip(&expected) {
                assert!(
                    (x - e).abs() <= 1e-7,
                    "origin x = {x0}: {xs:?}, expected {expected:?}"
                );
            }
            assert_hits_lie_on(&r, &line, &vase);
        }
    }
}

#[test]
fn test_hyperbola_meets_a_revolution_surface_where_its_crossing_is_alone() {
    // #310: the hyperbola runs close to the vase around three crossings,
    // and the closest seeds all crowded there, so the lone fourth
    // crossing was missed.
    let vase = vase();
    let hyperbola = Hyperbola::new(
        Point::new(-1.34857433611751, -1.251850001069137, 7.548267129822567),
        Vector::new(0.975886192657318, 0.1577412274699706, -0.15087691717778046),
        Vector::new(
            -0.1893834652173253,
            0.9555606775133391,
            -0.22591523784934603,
        ),
        2.1756456246775375,
        2.8592710888926964,
    )
    .unwrap();
    let expected = [-1.497277, -0.421211, 0.151598, 0.731247];
    let parameters = |r: &CurveSurfaceIntersection| {
        let mut ts: Vec<f64> = r.hits().iter().map(|h| h.curve_parameter).collect();
        ts.sort_by(f64::total_cmp);
        ts
    };
    let trimmed = TrimmedCurve::new(hyperbola, -3.0, 3.0).unwrap();
    for r in [
        hyperbola.intersect(&vase, TOL),
        vase.intersect(&hyperbola, TOL),
        trimmed.intersect(&vase, TOL),
    ] {
        let ts = parameters(&r);
        assert_eq!(ts.len(), 4, "{ts:?}");
        for (t, e) in ts.iter().zip(expected) {
            assert!((t - e).abs() <= 1e-6, "{ts:?}, expected {expected:?}");
        }
    }
    assert_hits_lie_on(&hyperbola.intersect(&vase, TOL), &hyperbola, &vase);
}

/// The hits of `r` lie within `range` of the curve's parameter.
fn assert_hits_within(r: &CurveSurfaceIntersection, (first, last): (f64, f64), case: &str) {
    for h in r.hits() {
        assert!(
            first <= h.curve_parameter && h.curve_parameter <= last,
            "{case}: hit at {} outside [{first}, {last}]",
            h.curve_parameter
        );
    }
}

#[test]
fn test_trimmed_curve_meets_a_swept_surface_only_within_its_range() {
    // #308: the numeric search let Newton run off the piece, and returned
    // hits of the basis curve beyond its ends.
    let pin = TrimmedCurve::segment(Point::new(0.0, 0.3, 2.0), Point::new(1.0, 0.3, 2.0)).unwrap();
    let ring = Circle::new(Point::ORIGIN, Vector::Z, 4.0).unwrap();
    let wall = ExtrusionSurface::new(ring.into(), Vector::Z).unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let tube = Circle::new(Point::new(3.0, 0.0, 2.0), Vector::Y, 1.0).unwrap();
    let ring_surface = RevolutionSurface::new(tube.into(), axis).unwrap();
    for r in [
        pin.intersect(&wall, TOL),
        wall.intersect(&pin, TOL),
        pin.intersect(&ring_surface, TOL),
        ring_surface.intersect(&pin, TOL),
    ] {
        assert!(r.hits().is_empty(), "{:?}", hit_xs(&r));
    }
    // A piece reaching the wall keeps its hit, at its end.
    let reaching =
        TrimmedCurve::segment(Point::new(0.0, 0.3, 2.0), Point::new(4.0, 0.3, 2.0)).unwrap();
    let x = (16.0f64 - 0.09).sqrt();
    let r = reaching.intersect(&wall, TOL);
    assert_eq!(hit_xs(&r).len(), 1);
    assert!((hit_xs(&r)[0] - x).abs() <= 1e-7, "{:?}", hit_xs(&r));
    // Segments and arcs across a vase and an extruded B-spline wall.
    let vase = vase();
    let profile = BSplineCurve::interpolate(
        &[
            Point::new(-3.0, 0.0, 0.0),
            Point::new(-1.0, 2.0, 0.0),
            Point::new(1.0, -1.0, 0.0),
            Point::new(3.0, 1.0, 0.0),
        ],
        3,
        PointParametrization::Centripetal,
    )
    .unwrap();
    let swept = ExtrusionSurface::new(profile.into(), Vector::Z).unwrap();
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let point = |next: &mut dyn FnMut() -> f64| {
        Point::new(8.0 * next() - 4.0, 8.0 * next() - 4.0, 10.0 * next())
    };
    for k in 0..200 {
        let (a, b) = (point(&mut next), point(&mut next));
        let Ok(segment) = TrimmedCurve::segment(a, b) else {
            continue;
        };
        let range = (segment.start_parameter(), segment.end_parameter());
        let circle = Circle::new(a, b - a, 0.5 + 3.0 * next()).unwrap();
        let start = 6.0 * next();
        let arc = TrimmedCurve::new(circle, start, start + 0.3 + 2.0 * next()).unwrap();
        let arc_range = (arc.start_parameter(), arc.end_parameter());
        for (name, surface) in [
            ("vase", Surface::from(vase.clone())),
            ("swept", swept.clone().into()),
        ] {
            let r = segment.intersect(&surface, TOL);
            assert_hits_within(&r, range, &format!("segment {k} on the {name}"));
            assert_hits_lie_on(&r, &segment, &surface);
            let r = arc.intersect(&surface, TOL);
            assert_hits_within(&r, arc_range, &format!("arc {k} on the {name}"));
            assert_hits_lie_on(&r, &arc, &surface);
        }
    }
}

#[test]
fn test_unbounded_line_meets_a_revolved_circle_like_the_equal_torus() {
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let tube = Circle::new(Point::new(3.0, 0.0, 0.0), Vector::Y, 1.0).unwrap();
    let ring = RevolutionSurface::new(tube.into(), axis).unwrap();
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let line = Line::new(Point::new(-10.0, 0.3, 0.2), Vector::X).unwrap();
    let expected = hit_xs(&line.intersect(&torus, TOL));
    assert_eq!(expected.len(), 4);
    let xs = hit_xs(&line.intersect(&ring, TOL));
    assert_eq!(xs.len(), 4, "{xs:?}");
    for (x, e) in xs.iter().zip(&expected) {
        assert!((x - e).abs() <= 1e-7, "{xs:?}, expected {expected:?}");
    }
}

#[test]
fn test_unbounded_line_meets_an_extrusion_far_from_its_origin() {
    // An extrusion is unbounded along its direction, so its box is
    // unknown, yet a line across it meets it only near its profile.
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 4.5, 3.5).unwrap();
    let wall = ExtrusionSurface::new(ellipse.into(), Vector::Z).unwrap();
    let crossing = |p: Point, d: Vector| {
        let line = Line::new(p, d).unwrap();
        let r = line.intersect(&wall, TOL);
        assert_hits_lie_on(&r, &line, &wall);
        hit_xs(&r)
    };
    // y = 0.5 on the ellipse x²/4.5² + y²/3.5² = 1.
    let x = 4.5 * (1.0 - (0.5f64 / 3.5).powi(2)).sqrt();
    for (p, d) in [
        (Point::new(-30.0, 0.5, 1.0), Vector::X),
        (Point::new(30.0, 0.5, 50.0), Vector::X),
        (Point::new(-30.0, 0.5, -20.0), Vector::new(1.0, 0.0, 1.0)),
    ] {
        let xs = crossing(p, d);
        assert_eq!(xs.len(), 2, "{p:?} along {d:?}: {xs:?}");
        assert!(
            (xs[0] + x).abs() <= 1e-7 && (xs[1] - x).abs() <= 1e-7,
            "{xs:?}"
        );
    }
    // A line along the direction, beside the wall, never meets it.
    assert!(crossing(Point::new(9.0, 0.0, 0.0), Vector::Z).is_empty());
    // A line that passes beside the prism holding the wall never meets it.
    assert!(crossing(Point::new(-30.0, 9.0, 0.0), Vector::X).is_empty());
}

#[test]
fn test_surfaces_on_one_control_net_meet_along_each_whole_edge_at_any_size() {
    // Every row and column of the net is straight, so the plain and the
    // rational surface share their four straight edges, and cross inside
    // along two more curves.
    let poles: Vec<Vec<Point>> = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| Point::new(i as f64, j as f64, ((i + 1) * (j + 2)) as f64 * 0.1))
                .collect()
        })
        .collect();
    let weights = vec![
        vec![1.0, 2.0, 1.0],
        vec![0.5, 3.0, 1.0],
        vec![1.0, 0.2, 4.0],
    ];
    let plain = BezierSurface::new(2, 2, poles.clone()).unwrap();
    let rational = BezierSurface::new_rational(2, 2, poles.clone(), weights).unwrap();
    let corner = |i: usize, j: usize| poles[i][j];
    let edges = [
        (corner(0, 0), corner(0, 2)),
        (corner(2, 0), corner(2, 2)),
        (corner(0, 0), corner(2, 0)),
        (corner(0, 2), corner(2, 2)),
    ];
    for (scale, shift) in [(1.0, 0.0), (2f64.powi(-20), 0.0), (1.0, 1e7)] {
        let offset = Vector::new(shift, -2.0 * shift, 1.5 * shift);
        // Scaled first, then moved: a point `p` goes to `p * scale + offset`.
        let similarity =
            Transform::scaling(Point::ORIGIN, scale).then(Transform::translation(offset));
        let moved = |p: Point| Point::ORIGIN + (p - Point::ORIGIN) * scale + offset;
        let tol = Tolerance::new(TOL.confusion * scale, TOL.angular, TOL.parametric);
        let a = plain.clone().transformed(similarity).unwrap();
        let b = rational.clone().transformed(similarity).unwrap();
        let result = Surface::from(a).intersect(&Surface::from(b), tol);
        assert_eq!(result.curves().len(), 6, "at scale {scale}, shift {shift}");
        for (start, end) in edges {
            let (start, end) = (moved(start), moved(end));
            let whole = result.curves().iter().any(|piece| {
                let (low, high) = piece.parameter_range;
                let (p, q) = (piece.curve.point_at(low), piece.curve.point_at(high));
                let close = |x: Point, y: Point| x.distance(y) <= 1e-6 * scale;
                (close(p, start) && close(q, end)) || (close(p, end) && close(q, start))
            });
            assert!(
                whole,
                "edge {start:?} to {end:?} at scale {scale}, shift {shift}"
            );
        }
    }
}
