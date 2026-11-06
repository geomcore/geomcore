use std::f64::consts::FRAC_PI_2;

use crate::geometry::curves::{
    BSplineCurve, BSplineCurve2D, Circle, Curve, ParametricCurve, ParametricCurve2D,
    PointParametrization, TrimmedCurve,
};
use crate::geometry::surfaces::{BSplineSurface, ParametricSurface, RevolutionSurface, Torus};
use crate::{
    Axis, BoundingBox, BoundingBox2D, HasBoundingBox, HasBoundingBox2D, Point, Point2D, Tolerance,
    Vector,
};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The rounding of evaluating a shape of unit size.
const ROUNDING: f64 = 1e-13;

/// The box of `count + 1` evenly spaced points of `point_at` over `[first, last]`.
fn sampled_box(
    point_at: impl Fn(f64) -> Point,
    (first, last): (f64, f64),
    count: usize,
) -> BoundingBox {
    BoundingBox::from_points(
        (0..=count).map(|k| point_at(first + (last - first) * (k as f64 / count as f64))),
    )
    .unwrap()
}

/// Checks that `tight` holds `sampled`, up to the rounding of evaluating
/// the samples, and that each side is within the confusion distance and
/// `slack` of it, the sampling's own shortfall.
fn assert_tight(tight: BoundingBox, sampled: BoundingBox, slack: f64) {
    let (low, high) = (tight.min_corner(), tight.max_corner());
    let (sampled_low, sampled_high) = (sampled.min_corner(), sampled.max_corner());
    for (side, sampled_side) in [
        (low.x, sampled_low.x),
        (low.y, sampled_low.y),
        (low.z, sampled_low.z),
    ] {
        assert!(side <= sampled_side + ROUNDING, "{tight:?} {sampled:?}");
        assert!(
            sampled_side - side <= TOL.confusion + slack,
            "{tight:?} {sampled:?}"
        );
    }
    for (side, sampled_side) in [
        (high.x, sampled_high.x),
        (high.y, sampled_high.y),
        (high.z, sampled_high.z),
    ] {
        assert!(side >= sampled_side - ROUNDING, "{tight:?} {sampled:?}");
        assert!(
            side - sampled_side <= TOL.confusion + slack,
            "{tight:?} {sampled:?}"
        );
    }
}

#[test]
fn an_interpolated_curve_is_boxed_to_its_extremes_not_its_poles() {
    // A hose through points that swing in y, whose control points
    // overshoot them.
    let points = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.2),
        Point::new(1.2, -0.1, 0.4),
        Point::new(2.0, 1.05, 0.1),
        Point::new(3.0, 0.0, -0.3),
    ];
    let hose = BSplineCurve::interpolate(&points, 3, PointParametrization::Chordal).unwrap();
    let loose = hose.bounding_box().unwrap();
    let tight = hose.tight_bounding_box(TOL).unwrap();
    let sampled = sampled_box(|t| hose.point_at(t), hose.bounds(), 200_000);
    assert_tight(tight, sampled, 1e-9);
    assert!(loose.diagonal().magnitude() > tight.diagonal().magnitude() + 0.1);
}

#[test]
fn a_rational_circle_is_boxed_as_the_circle() {
    let circle = Circle::new(Point::new(1.0, 2.0, 3.0), Vector::new(1.0, 1.0, 1.0), 2.0).unwrap();
    let spline = circle.to_bspline();
    let tight = spline.tight_bounding_box(TOL).unwrap();
    let exact = circle.bounding_box().unwrap();
    assert_tight(tight, exact, 0.0);
    // A range of it that cuts its knot spans.
    let (first, last) = (0.3, 0.3 + FRAC_PI_2);
    let arc = TrimmedCurve::new(Curve::BSpline(spline), first, last).unwrap();
    let sampled = sampled_box(|t| arc.point_at(t), (first, last), 200_000);
    assert_tight(arc.tight_bounding_box(TOL).unwrap(), sampled, 1e-9);
}

#[test]
fn a_plane_curve_is_boxed_in_the_plane() {
    let points = [
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(2.0, -1.0),
        Point2D::new(3.0, 0.5),
    ];
    let curve = BSplineCurve2D::interpolate(&points, 3, PointParametrization::Centripetal).unwrap();
    let tight = curve.tight_bounding_box(TOL).unwrap();
    let (first, last) = curve.bounds();
    let sampled = BoundingBox2D::from_points(
        (0..=200_000).map(|k| curve.point_at(first + (last - first) * (k as f64 / 200_000.0))),
    )
    .unwrap();
    let embedded = |b: BoundingBox2D| {
        BoundingBox::new(
            Point::new(b.min_corner().x, b.min_corner().y, 0.0),
            Point::new(b.max_corner().x, b.max_corner().y, 0.0),
        )
    };
    assert_tight(embedded(tight), embedded(sampled), 1e-9);
}

#[test]
fn a_spline_surface_is_boxed_to_its_extremes() {
    // A bicubic patch whose interior poles bulge far past the surface.
    let poles = (0..4)
        .map(|i| {
            (0..4)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    let bulge = if (1..3).contains(&i) && (1..3).contains(&j) {
                        3.0
                    } else {
                        0.0
                    };
                    Point::new(x, y, bulge + 0.3 * (x - y))
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
    let tight = surface.tight_bounding_box(TOL).unwrap();
    let count = 600;
    let sampled = BoundingBox::from_points(
        (0..=count)
            .flat_map(|i| {
                (0..=count).map(move |j| (i as f64 / count as f64, j as f64 / count as f64))
            })
            .map(|(u, v)| surface.point_at(u, v)),
    )
    .unwrap();
    assert_tight(tight, sampled, 1e-5);
    assert!(surface.bounding_box().unwrap().max_corner().z > tight.max_corner().z + 1.0);
}

#[test]
fn a_revolved_circle_is_boxed_as_its_torus() {
    // A circle of radius 0.5 at distance 2 from a tilted axis: the torus
    // about that axis.
    let direction = Vector::new(1.0, 2.0, 2.0) * (1.0 / 3.0);
    let center = Point::new(1.0, -1.0, 0.5);
    let axis = Axis::new(center, direction).unwrap();
    let across = direction.cross(Vector::Z).normalized().unwrap();
    let profile = Circle::new(center + across * 2.0, direction.cross(across), 0.5).unwrap();
    let surface = RevolutionSurface::new(Curve::Circle(profile), axis).unwrap();
    let torus = Torus::new(center, direction, 2.0, 0.5).unwrap();
    let tight = surface.tight_bounding_box(TOL).unwrap();
    assert_tight(tight, torus.bounding_box().unwrap(), 0.0);
    assert!(surface.bounding_box().unwrap().diagonal().magnitude() > tight.diagonal().magnitude());
}

#[test]
fn an_elementary_box_is_already_tight() {
    let circle = Circle::new(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0), 2.0).unwrap();
    assert_eq!(circle.tight_bounding_box(TOL), circle.bounding_box());
    let torus = Torus::new(Point::ORIGIN, Vector::X, 3.0, 1.0).unwrap();
    assert_eq!(torus.tight_bounding_box(TOL), torus.bounding_box());
}
