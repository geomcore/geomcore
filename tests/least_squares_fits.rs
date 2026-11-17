//! The least-squares fits of planes, circles, spheres and cylinders: exact
//! and noisy points, invariance under a similarity, the orientation each
//! fit takes from the points' order, and the configurations it refuses.

use geomcore::prelude::*;
use geomcore::{
    Axis, Circle, Circle2D, Cylinder, FitError, Plane, Point, Point2D, Sphere, Tolerance,
    Transform, Vector, Vector2D,
};
use std::f64::consts::{PI, TAU};

const TOL: Tolerance = Tolerance::DEFAULT;

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

/// A unit vector across `normal`, and the one across both.
fn across(normal: Vector) -> (Vector, Vector) {
    let hint = if normal.x.abs() < 0.9 {
        Vector::X
    } else {
        Vector::Y
    };
    let first = normal.cross(hint).normalized().unwrap();
    (first, normal.normalized().unwrap().cross(first))
}

/// `count` points on the circle of `center`, `normal` and `radius`, over
/// `span` radians from angle 0, counterclockwise about `normal`.
fn circle_points(
    center: Point,
    normal: Vector,
    radius: f64,
    span: f64,
    count: usize,
) -> Vec<Point> {
    let (u, v) = across(normal);
    (0..count)
        .map(|k| {
            let angle = span * k as f64 / count as f64;
            center + u * (radius * angle.cos()) + v * (radius * angle.sin())
        })
        .collect()
}

/// Points on the cylinder of `axis` and `radius`: `rings` rings of
/// `per_ring` points over `span` radians, from height 0 to `length`.
fn cylinder_points(
    axis: Axis,
    radius: f64,
    span: f64,
    length: f64,
    rings: usize,
    per_ring: usize,
) -> Vec<Point> {
    let direction = axis.direction().vector();
    let (u, v) = across(direction);
    let mut points = Vec::new();
    for ring in 0..rings {
        let height = length * ring as f64 / (rings - 1).max(1) as f64;
        for k in 0..per_ring {
            let angle = span * (k as f64 + 0.5 * (ring % 2) as f64) / per_ring as f64;
            points.push(
                axis.origin()
                    + direction * height
                    + u * (radius * angle.cos())
                    + v * (radius * angle.sin()),
            );
        }
    }
    points
}

fn assert_parallel(a: Vector, b: Vector, bound: f64) {
    let sine = a
        .normalized()
        .unwrap()
        .cross(b.normalized().unwrap())
        .magnitude();
    assert!(sine <= bound, "{a:?} and {b:?} are {sine} apart");
}

// --- Planes -----------------------------------------------------------------

#[test]
fn a_plane_through_exact_points_has_zero_residuals() {
    let normal = Vector::new(1.0, -2.0, 0.5);
    let (u, v) = across(normal);
    let origin = Point::new(3.0, -1.0, 7.0);
    let points: Vec<Point> = (0..20)
        .map(|k| origin + u * (k as f64 * 0.37).sin() * 4.0 + v * (k as f64 * 0.61).cos() * 2.0)
        .collect();
    let fit = Plane::fit(&points, TOL).unwrap();
    assert!(fit.maximum_residual() < 1e-13, "{}", fit.maximum_residual());
    assert_parallel(fit.shape().axis().direction().vector(), normal, 1e-14);
    assert!(fit.shape().signed_distance(origin).abs() < 1e-13);
}

#[test]
fn a_plane_through_noisy_points_reports_the_noise() {
    let mut noise = Noise(7);
    let normal = Vector::new(0.2, 0.3, 1.0).normalized().unwrap();
    let (u, v) = across(normal);
    let points: Vec<Point> = (0..400)
        .map(|_| {
            Point::ORIGIN
                + u * noise.next() * 10.0
                + v * noise.next() * 10.0
                + normal * noise.next() * 1e-3
        })
        .collect();
    let fit = Plane::fit(&points, TOL).unwrap();
    // Uniform noise over [-1e-3, 1e-3] has a root mean square of 1e-3/√3.
    let expected = 1e-3 / 3f64.sqrt();
    assert!((fit.root_mean_square_residual() / expected - 1.0).abs() < 0.1);
    assert!(fit.maximum_residual() <= 1.1e-3 && fit.maximum_residual() > 0.9e-3);
    assert_parallel(fit.shape().axis().direction().vector(), normal, 1e-4);
}

#[test]
fn a_plane_is_oriented_by_the_points_order() {
    let square = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    ];
    let fit = Plane::fit(&square, TOL).unwrap();
    assert_eq!(fit.shape().frame().z_direction(), Vector::Z);
    let reversed: Vec<Point> = square.iter().rev().copied().collect();
    let fit = Plane::fit(&reversed, TOL).unwrap();
    assert_eq!(fit.shape().frame().z_direction(), -Vector::Z);
}

// --- Circles in the plane ---------------------------------------------------

/// `count` points on the circle of `center` and `radius` in the plane,
/// over `span` radians from `start`.
fn circle_2d_points(
    center: Point2D,
    radius: f64,
    start: f64,
    span: f64,
    count: usize,
) -> Vec<Point2D> {
    (0..count)
        .map(|k| {
            let angle = start + span * k as f64 / (count - 1) as f64;
            center + Vector2D::new(angle.cos(), angle.sin()) * radius
        })
        .collect()
}

#[test]
fn a_circle_in_the_plane_through_exact_points_is_recovered() {
    let center = Point2D::new(-4.0, 2.5);
    for (span, count) in [(TAU * 0.9, 12), (PI / 3.0, 7), (PI / 18.0, 5)] {
        let points = circle_2d_points(center, 3.0, 0.4, span, count);
        let fit = Circle2D::fit(&points, TOL).unwrap();
        let circle = fit.shape();
        assert!(
            circle.center().distance(center) < 1e-10,
            "{span}: {circle:?}"
        );
        assert!((circle.radius() - 3.0).abs() < 1e-10, "{span}: {circle:?}");
        assert!(
            fit.maximum_residual() < 1e-12,
            "{span}: {}",
            fit.maximum_residual()
        );
    }
}

#[test]
fn three_points_give_their_circumcircle() {
    let points = [
        Point2D::new(0.0, 0.0),
        Point2D::new(4.0, 0.0),
        Point2D::new(0.0, 3.0),
    ];
    let fit = Circle2D::fit(&points, TOL).unwrap();
    assert!(fit.shape().center().distance(Point2D::new(2.0, 1.5)) < 1e-14);
    assert!((fit.shape().radius() - 2.5).abs() < 1e-14);
    assert!(fit.maximum_residual() < 1e-14);
}

#[test]
fn a_circle_in_the_plane_through_noisy_points_reports_the_noise() {
    let mut noise = Noise(11);
    let center = Point2D::new(1.0, 1.0);
    let points: Vec<Point2D> = circle_2d_points(center, 5.0, 0.0, PI, 200)
        .into_iter()
        .map(|p| {
            let radial = (p - center).normalized().unwrap();
            p + radial * (noise.next() * 1e-2)
        })
        .collect();
    let fit = Circle2D::fit(&points, TOL).unwrap();
    assert!(fit.shape().center().distance(center) < 3e-3);
    assert!((fit.shape().radius() - 5.0).abs() < 2e-3);
    let expected = 1e-2 / 3f64.sqrt();
    assert!((fit.root_mean_square_residual() / expected - 1.0).abs() < 0.15);
    assert!(fit.maximum_residual() <= 1.1e-2);
}

#[test]
fn a_circle_in_the_plane_runs_through_its_points_in_order() {
    let points = circle_2d_points(Point2D::new(2.0, 0.0), 1.0, 0.5, 2.0, 6);
    let fit = Circle2D::fit(&points, TOL).unwrap();
    let circle = fit.shape();
    assert!(circle.is_counterclockwise());
    assert!(circle.parameter_of(points[0], TOL).abs() < 1e-12);
    assert!((circle.parameter_of(points[5], TOL) - 2.0).abs() < 1e-12);

    let reversed: Vec<Point2D> = points.iter().rev().copied().collect();
    let fit = Circle2D::fit(&reversed, TOL).unwrap();
    assert!(!fit.shape().is_counterclockwise());
    assert!(fit.shape().parameter_of(reversed[0], TOL).abs() < 1e-12);
    assert!((fit.shape().parameter_of(reversed[5], TOL) - 2.0).abs() < 1e-12);
}

#[test]
fn a_circle_in_the_plane_moves_and_scales_with_its_points() {
    let mut noise = Noise(3);
    let points: Vec<Point2D> = circle_2d_points(Point2D::new(0.5, -0.5), 2.0, 0.0, 4.0, 30)
        .into_iter()
        .map(|p| p + Vector2D::new(noise.next(), noise.next()) * 1e-3)
        .collect();
    let reference = Circle2D::fit(&points, TOL).unwrap();
    let (sine, cosine) = 0.7f64.sin_cos();
    for (scale, offset) in [(1.0, 1e3), (1.0, 1e7), (1024.0, 0.0), (1.0 / 1024.0, 5.0)] {
        let moved = |p: Point2D| {
            Point2D::new(
                scale * (cosine * p.x - sine * p.y) + offset,
                scale * (sine * p.x + cosine * p.y) - offset,
            )
        };
        let tol = Tolerance {
            confusion: TOL.confusion * scale,
            ..TOL
        };
        let moved_points: Vec<Point2D> = points.iter().map(|&p| moved(p)).collect();
        let fit = Circle2D::fit(&moved_points, tol).unwrap();
        let rounding = 1e-15 * (offset + scale * 3.0);
        assert!(
            fit.shape()
                .center()
                .distance(moved(reference.shape().center()))
                < 1e3 * rounding
        );
        assert!((fit.shape().radius() - scale * reference.shape().radius()).abs() < 1e3 * rounding);
        let ratio =
            fit.root_mean_square_residual() / (scale * reference.root_mean_square_residual());
        assert!(
            (ratio - 1.0).abs() < 1e-6 * (1.0 + offset / scale),
            "{scale} {offset}: {ratio}"
        );
        let ratio = fit.maximum_residual() / (scale * reference.maximum_residual());
        assert!(
            (ratio - 1.0).abs() < 1e-6 * (1.0 + offset / scale),
            "{scale} {offset}: {ratio}"
        );
    }
}

// --- Circles in space -------------------------------------------------------

#[test]
fn a_circle_in_space_through_exact_points_is_recovered() {
    let center = Point::new(1.0, -2.0, 3.0);
    let normal = Vector::new(0.3, -0.4, 0.8);
    for (span, count) in [(TAU, 9), (PI / 2.0, 6), (PI / 12.0, 4), (TAU, 3)] {
        let points = circle_points(center, normal, 2.5, span, count);
        let fit = Circle::fit(&points, TOL).unwrap();
        let circle = fit.shape();
        assert!(
            circle.center().distance(center) < 1e-10,
            "{span}: {circle:?}"
        );
        assert!((circle.radius() - 2.5).abs() < 1e-10, "{span}: {circle:?}");
        assert_parallel(circle.normal().vector(), normal, 1e-10);
        assert!(
            fit.maximum_residual() < 1e-12,
            "{span}: {}",
            fit.maximum_residual()
        );
    }
}

#[test]
fn a_circle_in_space_through_noisy_points_reports_the_noise() {
    let mut noise = Noise(5);
    let center = Point::new(10.0, 0.0, 0.0);
    let normal = Vector::new(1.0, 1.0, 1.0);
    let points: Vec<Point> = circle_points(center, normal, 4.0, 3.0, 300)
        .into_iter()
        .map(|p| p + noise.vector() * 1e-3)
        .collect();
    let fit = Circle::fit(&points, TOL).unwrap();
    assert!(fit.shape().center().distance(center) < 1e-3);
    assert!((fit.shape().radius() - 4.0).abs() < 1e-3);
    assert_parallel(fit.shape().normal().vector(), normal, 1e-3);
    // Two of the three components of the noise are across the circle.
    let expected = 1e-3 * (2.0f64 / 3.0).sqrt();
    assert!((fit.root_mean_square_residual() / expected - 1.0).abs() < 0.15);
    assert!(fit.maximum_residual() <= 3f64.sqrt() * 1e-3);
}

#[test]
fn a_circle_in_space_off_its_plane_fits_the_distance_in_space() {
    // Points alternately above and below a circle: the plane and the
    // circle stay, every point is 0.1 away.
    let points: Vec<Point> = circle_points(Point::ORIGIN, Vector::Z, 1.0, TAU, 12)
        .into_iter()
        .enumerate()
        .map(|(k, p)| p + Vector::Z * if k % 2 == 0 { 0.1 } else { -0.1 })
        .collect();
    let fit = Circle::fit(&points, TOL).unwrap();
    assert!(fit.shape().center().distance(Point::ORIGIN) < 1e-12);
    assert!((fit.shape().radius() - 1.0).abs() < 1e-12);
    assert!((fit.root_mean_square_residual() - 0.1).abs() < 1e-12);
    assert!((fit.maximum_residual() - 0.1).abs() < 1e-12);
}

#[test]
fn a_circle_in_space_runs_through_its_points_in_order() {
    let normal = Vector::new(0.0, 0.6, 0.8);
    let points = circle_points(Point::new(1.0, 1.0, 1.0), normal, 2.0, 2.5, 6);
    let fit = Circle::fit(&points, TOL).unwrap();
    assert!((fit.shape().normal().dot(normal) - 1.0).abs() < 1e-12);
    assert!(fit.shape().parameter_of(points[0], TOL).abs() < 1e-12);
    let reversed: Vec<Point> = points.iter().rev().copied().collect();
    let fit = Circle::fit(&reversed, TOL).unwrap();
    assert!((fit.shape().normal().dot(normal) + 1.0).abs() < 1e-12);
    assert!(fit.shape().parameter_of(reversed[0], TOL).abs() < 1e-12);
}

// --- Spheres ----------------------------------------------------------------

/// `count` points on the sphere of `center` and `radius`, within `cap`
/// radians of its north pole.
fn sphere_points(center: Point, radius: f64, cap: f64, count: usize) -> Vec<Point> {
    (0..count)
        .map(|k| {
            // A spiral over the cap.
            let fraction = (k as f64 + 0.5) / count as f64;
            let polar = cap * fraction.sqrt();
            let azimuth = k as f64 * 2.399_963;
            let (sine, cosine) = polar.sin_cos();
            center + Vector::new(sine * azimuth.cos(), sine * azimuth.sin(), cosine) * radius
        })
        .collect()
}

#[test]
fn a_sphere_through_exact_points_is_recovered() {
    let center = Point::new(-3.0, 4.0, 0.5);
    for (cap, count) in [(PI, 50), (PI / 6.0, 20), (PI / 2.0, 4)] {
        let points = sphere_points(center, 7.0, cap, count);
        let fit = Sphere::fit(&points, TOL).unwrap();
        assert!(
            fit.shape().center().distance(center) < 1e-10,
            "{cap}: {:?}",
            fit.shape()
        );
        assert!(
            (fit.shape().radius() - 7.0).abs() < 1e-10,
            "{cap}: {:?}",
            fit.shape()
        );
        assert!(
            fit.maximum_residual() < 1e-12,
            "{cap}: {}",
            fit.maximum_residual()
        );
    }
}

#[test]
fn a_sphere_through_noisy_points_reports_the_noise() {
    let mut noise = Noise(13);
    let center = Point::new(0.0, 0.0, 100.0);
    let points: Vec<Point> = sphere_points(center, 3.0, PI, 500)
        .into_iter()
        .map(|p| p + (p - center).normalized().unwrap() * (noise.next() * 1e-3))
        .collect();
    let fit = Sphere::fit(&points, TOL).unwrap();
    assert!(fit.shape().center().distance(center) < 1e-4);
    assert!((fit.shape().radius() - 3.0).abs() < 1e-4);
    let expected = 1e-3 / 3f64.sqrt();
    assert!((fit.root_mean_square_residual() / expected - 1.0).abs() < 0.1);
    assert!(fit.maximum_residual() <= 1.1e-3);
}

// --- Cylinders --------------------------------------------------------------

fn assert_cylinder(fit: &geomcore::CylinderFit, axis: Axis, radius: f64, bound: f64) {
    let cylinder = fit.shape();
    assert!((cylinder.radius() - radius).abs() < bound, "{cylinder:?}");
    assert_parallel(
        cylinder.axis().direction().vector(),
        axis.direction().vector(),
        bound,
    );
    assert!(
        axis.distance(cylinder.frame().origin()) < bound,
        "{cylinder:?}"
    );
}

#[test]
fn a_cylinder_through_exact_points_is_recovered() {
    let axis = Axis::new(Point::new(1.0, 2.0, -1.0), Vector::new(0.3, -0.5, 0.8)).unwrap();
    let cases = [
        // Whole rings along a length near the diameter.
        (TAU, 4.0, 4, 8),
        // A long thin tube.
        (TAU, 50.0, 6, 6),
        // A short wide ring band.
        (TAU, 0.2, 2, 16),
        // A quarter of the cylinder.
        (PI / 2.0, 3.0, 5, 6),
        // A sixth of it, longer than it is wide.
        (PI / 3.0, 6.0, 6, 5),
        // A sixth of it, wider than it is long.
        (PI / 3.0, 0.5, 3, 8),
    ];
    for (span, length, rings, per_ring) in cases {
        let points = cylinder_points(axis, 2.0, span, length, rings, per_ring);
        let fit = Cylinder::fit(&points, TOL).unwrap();
        assert!(
            fit.maximum_residual() < 1e-11,
            "{span} {length}: {} {:?}",
            fit.maximum_residual(),
            fit.shape()
        );
        assert_cylinder(&fit, axis, 2.0, 1e-9);
    }
}

#[test]
fn five_points_on_a_cylinder_fit_one_through_them() {
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let points = [
        Point::new(1.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.3),
        Point::new(-1.0, 0.0, 1.1),
        Point::new(0.0, -1.0, 0.4),
        Point::new(0.6, 0.8, 1.7),
    ];
    let fit = Cylinder::fit(&points, TOL).unwrap();
    assert!(fit.maximum_residual() < 1e-10, "{}", fit.maximum_residual());
    assert_cylinder(&fit, axis, 1.0, 1e-9);
}

#[test]
fn a_cylinder_through_noisy_points_reports_the_noise() {
    let mut noise = Noise(17);
    let axis = Axis::new(Point::new(5.0, 5.0, 5.0), Vector::new(1.0, 0.0, 1.0)).unwrap();
    for (span, length) in [(TAU, 10.0), (PI / 2.0, 5.0)] {
        let points: Vec<Point> = cylinder_points(axis, 3.0, span, length, 20, 20)
            .into_iter()
            .map(|p| {
                let radial = (p - axis.closest_point(p)).normalized().unwrap();
                p + radial * (noise.next() * 1e-3)
            })
            .collect();
        let fit = Cylinder::fit(&points, TOL).unwrap();
        assert_cylinder(&fit, axis, 3.0, 2e-3);
        let expected = 1e-3 / 3f64.sqrt();
        let rms = fit.root_mean_square_residual();
        assert!(
            (rms / expected - 1.0).abs() < 0.1,
            "{span}: {rms} {:?}",
            fit.shape()
        );
        // The noise, and the fit tilted by it towards the ends.
        assert!(
            fit.maximum_residual() <= 1.25e-3,
            "{span}: {}",
            fit.maximum_residual()
        );
    }
}

#[test]
fn a_cylinder_is_placed_and_oriented_by_its_points() {
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let points = cylinder_points(axis, 1.0, TAU, 2.0, 3, 6);
    let fit = Cylinder::fit(&points, TOL).unwrap();
    let frame = fit.shape().frame();
    // On the axis, level with the centroid; turning with the points.
    assert!(frame.origin().distance(Point::new(0.0, 0.0, 1.0)) < 1e-12);
    assert!((frame.z_direction().dot(Vector::Z) - 1.0).abs() < 1e-12);
    // Towards the first point.
    let first = (points[0] - frame.origin())
        .component_perpendicular_to(frame.z_direction())
        .normalized()
        .unwrap();
    assert!((frame.x_direction().dot(first) - 1.0).abs() < 1e-12);
    let reversed: Vec<Point> = points.iter().rev().copied().collect();
    let fit = Cylinder::fit(&reversed, TOL).unwrap();
    assert!((fit.shape().frame().z_direction().dot(Vector::Z) + 1.0).abs() < 1e-12);
}

// --- Similarity invariance --------------------------------------------------

/// The similarities the fits are checked under: a power-of-two scale, a
/// rotation and a translation, alone and together.
fn similarities() -> Vec<Transform> {
    let turn = Transform::rotation(
        Axis::new(Point::new(1.0, 2.0, 3.0), Vector::new(1.0, -1.0, 2.0)).unwrap(),
        0.9,
    );
    vec![
        turn,
        Transform::translation(Vector::new(1e3, -1e3, 1e3)),
        Transform::translation(Vector::new(1e7, 1e7, -1e7)),
        Transform::scaling(Point::ORIGIN, 1024.0).then(turn),
        Transform::scaling(Point::ORIGIN, 1.0 / 1024.0),
    ]
}

/// Checks that `fit_moved` of the moved points is the moved fit of the
/// points, comparing with `same`, and that the residuals scale with them.
fn check_similarity_invariance<S>(
    points: &[Point],
    fit: impl Fn(&[Point], Tolerance) -> geomcore::Fit<S>,
    same: impl Fn(&S, &S, f64),
) where
    S: Transformable + Copy,
{
    let reference = fit(points, TOL);
    for transform in similarities() {
        let scale = transform.scale_factor();
        let tol = Tolerance {
            confusion: TOL.confusion * scale,
            ..TOL
        };
        let moved: Vec<Point> = points.iter().map(|&p| transform.apply_point(p)).collect();
        let fitted = fit(&moved, tol);
        let expected = reference.shape().transformed(transform).unwrap();
        // The rounding of the moved coordinates, magnified by the fit.
        let size = moved
            .iter()
            .map(|&p| (p - Point::ORIGIN).magnitude())
            .fold(0.0, f64::max);
        same(fitted.shape(), &expected, 1e-11 * size.max(scale));
        for (got, want) in [
            (
                fitted.root_mean_square_residual(),
                reference.root_mean_square_residual(),
            ),
            (fitted.maximum_residual(), reference.maximum_residual()),
        ] {
            assert!(
                (got - scale * want).abs() <= 1e-6 * scale * want + 1e-11 * size,
                "{got} against {} for {transform:?}",
                scale * want
            );
        }
    }
}

#[test]
fn a_plane_fit_is_similarity_invariant() {
    let mut noise = Noise(19);
    let points: Vec<Point> = (0..40)
        .map(|_| {
            Point::new(
                noise.next() * 3.0,
                noise.next() * 2.0,
                0.5 + noise.next() * 0.01,
            )
        })
        .collect();
    check_similarity_invariance(
        &points,
        |p, tol| Plane::fit(p, tol).unwrap(),
        |a, b, bound| {
            assert!(a.frame().origin().distance(b.frame().origin()) < bound);
            assert!((a.frame().z_direction().dot(b.frame().z_direction()) - 1.0).abs() < 1e-12);
        },
    );
}

#[test]
fn a_circle_fit_is_similarity_invariant() {
    let mut noise = Noise(23);
    let points: Vec<Point> = circle_points(
        Point::new(1.0, 0.0, 0.0),
        Vector::new(0.0, 1.0, 1.0),
        2.0,
        4.0,
        25,
    )
    .into_iter()
    .map(|p| p + noise.vector() * 1e-3)
    .collect();
    check_similarity_invariance(
        &points,
        |p, tol| Circle::fit(p, tol).unwrap(),
        |a, b, bound| {
            assert!(a.center().distance(b.center()) < bound, "{a:?} {b:?}");
            assert!((a.radius() - b.radius()).abs() < bound);
            assert!((a.normal().dot(b.normal()) - 1.0).abs() < 1e-11);
        },
    );
}

#[test]
fn a_sphere_fit_is_similarity_invariant() {
    let mut noise = Noise(29);
    let points: Vec<Point> = sphere_points(Point::new(0.0, 1.0, 2.0), 3.0, PI / 3.0, 40)
        .into_iter()
        .map(|p| p + noise.vector() * 1e-3)
        .collect();
    check_similarity_invariance(
        &points,
        |p, tol| Sphere::fit(p, tol).unwrap(),
        |a, b, bound| {
            assert!(a.center().distance(b.center()) < bound, "{a:?} {b:?}");
            assert!((a.radius() - b.radius()).abs() < bound);
        },
    );
}

#[test]
fn a_cylinder_fit_is_similarity_invariant() {
    let mut noise = Noise(31);
    let axis = Axis::new(Point::new(0.0, 0.0, 1.0), Vector::new(1.0, 2.0, 0.5)).unwrap();
    let points: Vec<Point> = cylinder_points(axis, 1.5, 2.0, 4.0, 6, 8)
        .into_iter()
        .map(|p| p + noise.vector() * 1e-3)
        .collect();
    check_similarity_invariance(
        &points,
        |p, tol| Cylinder::fit(p, tol).unwrap(),
        |a, b, bound| {
            assert!(
                a.frame().origin().distance(b.frame().origin()) < bound,
                "{a:?} {b:?}"
            );
            assert!((a.radius() - b.radius()).abs() < bound);
            assert!((a.axis().direction().dot(b.axis().direction()) - 1.0).abs() < 1e-11);
        },
    );
}

// --- Refusals ---------------------------------------------------------------

#[test]
fn each_fit_refuses_fewer_points_than_its_shape_needs() {
    let points = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
        Point::new(0.0, 0.0, 1.0),
    ];
    let too_few = |required, given| FitError::TooFewPoints { required, given };
    assert_eq!(Plane::fit(&points[..2], TOL), Err(too_few(3, 2)));
    assert_eq!(Circle::fit(&points[..2], TOL), Err(too_few(3, 2)));
    assert_eq!(
        Circle2D::fit(&[Point2D::ORIGIN, Point2D::new(1.0, 0.0)], TOL),
        Err(too_few(3, 2))
    );
    assert_eq!(Sphere::fit(&points[..3], TOL), Err(too_few(4, 3)));
    assert_eq!(Cylinder::fit(&points, TOL), Err(too_few(5, 4)));
    assert_eq!(Plane::fit(&[], TOL), Err(too_few(3, 0)));
    assert_eq!(
        too_few(5, 4).to_string(),
        "the fit needs at least 5 points, 4 were given"
    );
}

#[test]
fn each_fit_refuses_a_point_that_is_not_finite() {
    let mut points = circle_points(Point::ORIGIN, Vector::Z, 1.0, TAU, 6);
    points.push(Point::new(0.0, f64::NAN, 1.0));
    assert_eq!(Plane::fit(&points, TOL), Err(FitError::NonFinite));
    assert_eq!(Circle::fit(&points, TOL), Err(FitError::NonFinite));
    assert_eq!(Sphere::fit(&points, TOL), Err(FitError::NonFinite));
    assert_eq!(Cylinder::fit(&points, TOL), Err(FitError::NonFinite));
    let plane_points = [
        Point2D::ORIGIN,
        Point2D::new(1.0, 0.0),
        Point2D::new(f64::INFINITY, 1.0),
    ];
    assert_eq!(Circle2D::fit(&plane_points, TOL), Err(FitError::NonFinite));
}

#[test]
fn each_fit_refuses_confused_points() {
    let points: Vec<Point> = (0..6)
        .map(|k| Point::new(1.0, 2.0, 3.0) + Vector::new(1.0, -1.0, 0.5) * (k as f64 * 1e-9))
        .collect();
    assert_eq!(Plane::fit(&points, TOL), Err(FitError::ConfusedPoints));
    assert_eq!(Circle::fit(&points, TOL), Err(FitError::ConfusedPoints));
    assert_eq!(Sphere::fit(&points, TOL), Err(FitError::ConfusedPoints));
    assert_eq!(Cylinder::fit(&points, TOL), Err(FitError::ConfusedPoints));
    let plane_points = [Point2D::new(5.0, 5.0); 4];
    assert_eq!(
        Circle2D::fit(&plane_points, TOL),
        Err(FitError::ConfusedPoints)
    );
}

#[test]
fn each_fit_refuses_collinear_points() {
    // Within the confusion tolerance of a line, not on it.
    let mut noise = Noise(37);
    let direction = Vector::new(1.0, 2.0, 2.0) / 3.0;
    let points: Vec<Point> = (0..10)
        .map(|k| Point::new(1e3, 0.0, 0.0) + direction * k as f64 + noise.vector() * 3e-8)
        .collect();
    assert_eq!(Plane::fit(&points, TOL), Err(FitError::CollinearPoints));
    assert_eq!(Circle::fit(&points, TOL), Err(FitError::CollinearPoints));
    assert_eq!(Sphere::fit(&points, TOL), Err(FitError::CollinearPoints));
    assert_eq!(Cylinder::fit(&points, TOL), Err(FitError::CollinearPoints));
    let plane_points: Vec<Point2D> = (0..10)
        .map(|k| Point2D::new(k as f64, 0.5 * k as f64 + noise.next() * 3e-8))
        .collect();
    assert_eq!(
        Circle2D::fit(&plane_points, TOL),
        Err(FitError::CollinearPoints)
    );

    // Bent a little further off the line, a circle fits them.
    let bent: Vec<Point2D> = (0..10)
        .map(|k| {
            let x = k as f64 - 4.5;
            Point2D::new(x, 1e-5 * x * x)
        })
        .collect();
    let fit = Circle2D::fit(&bent, TOL).unwrap();
    assert!(
        (fit.shape().radius() / 5e4 - 1.0).abs() < 1e-3,
        "{:?}",
        fit.shape()
    );
}

#[test]
fn a_sphere_and_a_cylinder_refuse_coplanar_points() {
    let mut noise = Noise(41);
    let normal = Vector::new(0.0, 0.6, 0.8);
    let points: Vec<Point> = circle_points(Point::new(3.0, 3.0, 3.0), normal, 2.0, TAU, 12)
        .into_iter()
        .map(|p| p + normal * (noise.next() * 3e-8))
        .collect();
    assert_eq!(Sphere::fit(&points, TOL), Err(FitError::CoplanarPoints));
    assert_eq!(Cylinder::fit(&points, TOL), Err(FitError::CoplanarPoints));
    // A circle and a plane fit them.
    assert!(Circle::fit(&points, TOL).unwrap().maximum_residual() < 1e-7);
    assert!(Plane::fit(&points, TOL).unwrap().maximum_residual() < 1e-7);
}

#[test]
fn the_refusals_follow_the_confusion_tolerance() {
    // Points 1e-5 off a plane are coplanar for a coarse tolerance only.
    let points: Vec<Point> = sphere_points(Point::ORIGIN, 1e3, 1e-4, 20);
    let coarse = Tolerance {
        confusion: 1e-3,
        ..TOL
    };
    assert_eq!(Sphere::fit(&points, coarse), Err(FitError::CoplanarPoints));
    assert!(Sphere::fit(&points, TOL).is_ok());
}
