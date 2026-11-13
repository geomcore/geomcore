//! Building a B-spline the textbook way and reading its shape: flat knot
//! vectors, clamped uniform knots, single-pole and single-weight edits,
//! and geometric continuity at a knot.

use geomcore::prelude::*;
use geomcore::{
    BSplineConstructionError, BSplineCurve, BSplineCurve2D, Continuity, Point, Point2D,
    PoleEditError, Tolerance,
};

const TOL: Tolerance = Tolerance::DEFAULT;

fn wavy_poles() -> Vec<Point> {
    vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.5),
        Point::new(3.0, 2.0, 1.0),
        Point::new(4.0, 0.0, 0.3),
        Point::new(5.0, 1.0, 2.0),
        Point::new(6.0, -1.0, 0.0),
    ]
}

#[test]
fn flat_knots_round_trip_every_kind_of_curve() {
    let clamped = BSplineCurve::new(
        3,
        wavy_poles(),
        vec![0.0, 0.4, 0.7, 1.0],
        vec![4, 1, 1, 4],
        false,
    )
    .unwrap();
    assert_eq!(
        clamped.flat_knots(),
        &[0.0, 0.0, 0.0, 0.0, 0.4, 0.7, 1.0, 1.0, 1.0, 1.0]
    );
    let rational = BSplineCurve::new_rational(
        2,
        wavy_poles(),
        vec![1.0, 0.5, 2.0, 1.0, 0.7, 1.0],
        vec![0.0, 1.0, 2.0, 3.0, 4.0],
        vec![3, 1, 1, 1, 3],
        false,
    )
    .unwrap();
    let periodic = BSplineCurve::new(
        2,
        wavy_poles(),
        vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        vec![1, 1, 1, 1, 1, 1, 1],
        true,
    )
    .unwrap();
    for curve in [clamped, rational, periodic] {
        let rebuilt = BSplineCurve::from_flat_knots(
            curve.degree(),
            curve.poles().to_vec(),
            curve.flat_knots().to_vec(),
            curve.weights().map(<[f64]>::to_vec),
            curve.is_periodic(),
        )
        .unwrap();
        assert_eq!(rebuilt, curve);
    }
}

#[test]
fn an_unclamped_flat_knot_vector_builds_the_uniform_textbook_curve() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 2.0),
        Point2D::new(4.0, 0.0),
    ];
    let flat = (0..8).map(f64::from).collect();
    let curve = BSplineCurve2D::from_flat_knots(3, poles, flat, None, false).unwrap();
    assert_eq!(curve.bounds(), (3.0, 4.0));
    // The uniform cubic starts at (P0 + 4 P1 + P2) / 6.
    let start = curve.point_at(3.0);
    assert!(start.distance(Point2D::new(7.0 / 6.0, 10.0 / 6.0)) < 1e-15);
}

#[test]
fn a_decreasing_or_short_flat_knot_vector_is_refused() {
    let decreasing = vec![0.0, 0.0, 0.0, 0.6, 0.4, 1.0, 1.0, 1.0, 1.0, 1.0];
    assert_eq!(
        BSplineCurve::from_flat_knots(3, wavy_poles(), decreasing, None, false),
        Err(BSplineConstructionError::KnotsNotIncreasing)
    );
    let short = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    assert_eq!(
        BSplineCurve::from_flat_knots(3, wavy_poles(), short, None, false),
        Err(BSplineConstructionError::PoleCountMismatch)
    );
}

#[test]
fn clamped_uniform_spaces_the_knots_evenly_over_the_unit_interval() {
    let curve = BSplineCurve::clamped_uniform(3, wavy_poles()).unwrap();
    assert_eq!(curve.knots(), &[0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]);
    assert_eq!(curve.multiplicities(), &[4, 1, 1, 4]);
    assert_eq!(curve.point_at(0.0), wavy_poles()[0]);
    assert_eq!(curve.point_at(1.0), wavy_poles()[5]);
    // As many poles as the degree plus one: the Bézier curve.
    let bezier = BSplineCurve::clamped_uniform(5, wavy_poles()).unwrap();
    assert_eq!(bezier.knots(), &[0.0, 1.0]);
    assert_eq!(
        BSplineCurve::clamped_uniform(6, wavy_poles()),
        Err(BSplineConstructionError::PoleCountMismatch)
    );
    assert_eq!(
        BSplineCurve::clamped_uniform(0, wavy_poles()),
        Err(BSplineConstructionError::InvalidDegree)
    );
}

#[test]
fn moving_one_pole_changes_only_its_spans() {
    let curve = BSplineCurve::clamped_uniform(3, wavy_poles()).unwrap();
    let moved = curve.with_pole(5, Point::new(6.0, 3.0, 0.0)).unwrap();
    assert_eq!(moved.poles()[5], Point::new(6.0, 3.0, 0.0));
    assert_eq!(moved.knots(), curve.knots());
    // Pole 5 has no say before the last span starts at 2/3.
    assert_eq!(moved.point_at(0.5), curve.point_at(0.5));
    assert_ne!(moved.point_at(0.9), curve.point_at(0.9));
    assert_eq!(curve.poles(), wavy_poles().as_slice());
    assert_eq!(
        curve.with_pole(6, Point::ORIGIN).err(),
        Some(PoleEditError::IndexOutOfRange)
    );
}

#[test]
fn setting_a_weight_makes_the_curve_rational_and_pulls_it_toward_the_pole() {
    let curve = BSplineCurve::clamped_uniform(3, wavy_poles()).unwrap();
    let heavy = curve.with_weight(2, 4.0).unwrap();
    assert_eq!(heavy.weights(), Some(&[1.0, 1.0, 4.0, 1.0, 1.0, 1.0][..]));
    let pole = wavy_poles()[2];
    let parameter = curve.project_point(pole, TOL).parameter;
    assert!(heavy.point_at(parameter).distance(pole) < curve.point_at(parameter).distance(pole));
    let back = heavy.with_weight(2, 1.0).unwrap();
    assert!(back.is_rational());
    assert!(back.point_at(0.3).distance(curve.point_at(0.3)) < 1e-15);
    for weight in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            curve.with_weight(2, weight).err(),
            Some(PoleEditError::InvalidWeight)
        );
    }
    assert_eq!(
        curve.with_weight(9, 1.0).err(),
        Some(PoleEditError::IndexOutOfRange)
    );
}

/// The quadratic through `poles` with a double knot at 1: C0 there by
/// the knot.
fn double_knot_quadratic(poles: [(f64, f64); 5]) -> BSplineCurve2D {
    let poles = poles.iter().map(|&(x, y)| Point2D::new(x, y)).collect();
    BSplineCurve2D::new(2, poles, vec![0.0, 1.0, 2.0], vec![3, 2, 3], false).unwrap()
}

#[test]
fn geometric_continuity_reads_the_shape_where_the_knot_says_c0() {
    // The tangent turns: G0.
    let kinked =
        double_knot_quadratic([(0.0, 0.0), (1.0, 1.0), (2.0, 1.0), (3.0, 2.0), (4.0, 0.0)]);
    assert_eq!(kinked.continuity(1.0, TOL), Continuity::Order(0));
    assert_eq!(kinked.geometric_continuity(1.0, TOL), Continuity::Order(0));
    // Poles in line around the knot but unevenly spaced: G1.
    let tangent =
        double_knot_quadratic([(0.0, 0.0), (1.0, 1.0), (2.0, 1.0), (4.0, 1.0), (5.0, 0.0)]);
    assert_eq!(tangent.continuity(1.0, TOL), Continuity::Order(0));
    assert_eq!(tangent.geometric_continuity(1.0, TOL), Continuity::Order(1));
    // Away from the knot, the shape is as smooth as the representation.
    assert_eq!(tangent.geometric_continuity(0.5, TOL), Continuity::Smooth);
}

#[test]
fn a_knot_inserted_up_to_the_degree_leaves_the_shape_curvature_continuous() {
    let smooth = BSplineCurve::clamped_uniform(3, wavy_poles()).unwrap();
    let split = smooth.insert_knot(0.5, 3).unwrap();
    assert_eq!(split.continuity(0.5, TOL), Continuity::Order(0));
    assert_eq!(split.geometric_continuity(0.5, TOL), Continuity::Order(2));
    // Every knot the curve is already C2 at reports its parametric order.
    let fine = smooth.insert_knot(0.5, 1).unwrap();
    assert_eq!(fine.geometric_continuity(0.5, TOL), Continuity::Order(2));
}

#[test]
fn a_periodic_seam_is_compared_like_any_other_knot() {
    let square = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
    ];
    let poles = square.iter().map(|&(x, y)| Point2D::new(x, y)).collect();
    let knots = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let curve = BSplineCurve2D::new(2, poles, knots, vec![2, 2, 2, 2, 2], true).unwrap();
    // A double knot of a quadratic at every edge midpoint: C0 by the
    // representation, while the curve passes the midpoint between two
    // mirror-image pieces, curvature-continuous, the seam included.
    for knot in [0.0, 1.0, 4.0] {
        assert_eq!(curve.continuity(knot, TOL), Continuity::Order(0));
        assert_eq!(curve.geometric_continuity(knot, TOL), Continuity::Order(2));
    }
    // Sliding the next pole along the seam tangent keeps the tangent and
    // changes the curvature of the first piece alone.
    let slid = curve.with_pole(1, Point2D::new(1.0, 2.0)).unwrap();
    assert_eq!(slid.geometric_continuity(4.0, TOL), Continuity::Order(1));
}

#[test]
fn a_side_with_no_speed_keeps_the_parametric_order() {
    // A repeated pole at the knot stops the curve there from both sides.
    let stopping =
        double_knot_quadratic([(0.0, 0.0), (1.0, 1.0), (1.0, 1.0), (1.0, 1.0), (2.0, 0.0)]);
    assert_eq!(
        stopping.geometric_continuity(1.0, TOL),
        Continuity::Order(0)
    );
}

/// The point a basis gives: the sum of `values[j]` times pole
/// `first + j`, wrapping past the last pole.
fn point_from_basis(curve: &BSplineCurve, u: f64) -> Point {
    let (first, values) = curve.basis_functions(u);
    let poles = curve.poles();
    let mut sum = [0.0; 3];
    for (j, value) in values.iter().enumerate() {
        let pole = poles[(first + j) % poles.len()];
        sum[0] += value * pole.x;
        sum[1] += value * pole.y;
        sum[2] += value * pole.z;
    }
    Point::new(sum[0], sum[1], sum[2])
}

#[test]
fn the_basis_functions_rebuild_the_point_and_the_span_holds_it() {
    let clamped = BSplineCurve::new(
        3,
        wavy_poles(),
        vec![0.0, 0.4, 0.7, 1.0],
        vec![4, 1, 1, 4],
        false,
    )
    .unwrap();
    let periodic = BSplineCurve::new(
        3,
        wavy_poles(),
        vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        vec![1; 7],
        true,
    )
    .unwrap();
    for curve in [&clamped, &periodic] {
        let (first, last) = curve.bounds();
        let flat = curve.flat_knots();
        for k in 0..=40 {
            let u = first + (last - first) * f64::from(k) / 40.0;
            let (_, values) = curve.basis_functions(u);
            assert_eq!(values.len(), curve.degree() + 1);
            assert!((values.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            assert!(point_from_basis(curve, u).distance(curve.point_at(u)) < 1e-12);
            let s = curve.span_index(u);
            assert!(flat[s] <= u && flat[s] < flat[s + 1], "u = {u}, s = {s}");
            assert!(u < flat[s + 1] || u == last, "u = {u}, s = {s}");
        }
    }
    // Past the bounds a non-periodic curve clamps, as evaluation does,
    // and a periodic one wraps.
    assert_eq!(clamped.span_index(-1.0), clamped.span_index(0.0));
    assert_eq!(clamped.basis_functions(2.0), clamped.basis_functions(1.0));
    assert_eq!(periodic.span_index(7.5), periodic.span_index(1.5));
    assert!(point_from_basis(&periodic, 7.5).distance(periodic.point_at(7.5)) < 1e-12);
    // The 2D tier gives the same basis.
    let flat_poles: Vec<Point2D> = wavy_poles()
        .iter()
        .map(|p| Point2D::new(p.x, p.y))
        .collect();
    let clamped_2d = BSplineCurve2D::new(
        3,
        flat_poles,
        vec![0.0, 0.4, 0.7, 1.0],
        vec![4, 1, 1, 4],
        false,
    )
    .unwrap();
    assert_eq!(
        clamped_2d.basis_functions(0.55),
        clamped.basis_functions(0.55)
    );
    assert_eq!(clamped_2d.span_index(0.55), clamped.span_index(0.55));
}
