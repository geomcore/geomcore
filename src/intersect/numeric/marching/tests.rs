use super::*;

#[test]
fn test_marching_plane_plane_line() {
    use crate::Plane;
    let tol = Tolerance::DEFAULT;
    let xy = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let zy = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    let curves = marched_curves(&xy, &zy, tol);
    assert_eq!(curves.len(), 1);
    // The traced line runs along Y through the origin.
    for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let p = curves[0].eval_point(t);
        assert!(p.x.abs() < 1e-6, "{p:?}");
        assert!(p.z.abs() < 1e-6, "{p:?}");
    }
}

#[test]
fn test_marching_plane_sphere_circle() {
    use crate::{Plane, Sphere};
    let tol = Tolerance::DEFAULT;
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let curves = marched_curves(&plane, &sphere, tol);
    assert_eq!(curves.len(), 1);
    // Approximate fit: residuals within marching scale (not confusion).
    for t in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0] {
        let p = curves[0].eval_point(t);
        assert!(plane.project_point(p, tol).distance < 1e-3, "{p:?}");
        assert!(sphere.project_point(p, tol).distance < 1e-3, "{p:?}");
    }
}

#[test]
fn test_marching_loop_is_traced_once() {
    use crate::{Plane, Sphere};
    let tol = Tolerance::DEFAULT;
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let curves = marched_curves(&plane, &sphere, tol);
    let curve = assert_one_curve_on_both(&curves, &plane, &sphere);
    let (t0, t1) = curve.bounds();
    let n = 2000;
    let length: f64 = (0..n)
        .map(|k| {
            let at = |k: usize| curve.eval_point(t0 + (t1 - t0) * k as f64 / n as f64);
            at(k).distance(at(k + 1))
        })
        .sum();
    let circumference = 4.0 * std::f64::consts::PI;
    assert!(
        (length - circumference).abs() < 1e-3 * circumference,
        "{length}"
    );
}

/// Asserts one traced curve whose samples lie on both surfaces.
fn assert_one_curve_on_both<A, B>(curves: &[BSplineCurve], a: &A, b: &B) -> BSplineCurve
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let tol = Tolerance::DEFAULT;
    assert_eq!(curves.len(), 1);
    let curve = curves[0].clone();
    let (t0, t1) = curve.bounds();
    for k in 0..=40 {
        let p = curve.eval_point(t0 + (t1 - t0) * k as f64 / 40.0);
        assert!(a.project_point(p, tol).distance < 1e-3, "{p:?}");
        assert!(b.project_point(p, tol).distance < 1e-3, "{p:?}");
    }
    curve
}

#[test]
fn test_marching_plane_sphere_off_grid_heights() {
    use crate::{Plane, Sphere};
    let tol = Tolerance::DEFAULT;
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    // No grid row of either surface lands on these circles.
    for h in [0.3, 0.8, -0.45, 0.123] {
        let plane = Plane::new(Point::new(0.0, 0.0, h), Vector::Z).unwrap();
        let curves = marched_curves(&plane, &sphere, tol);
        let curve = assert_one_curve_on_both(&curves, &plane, &sphere);
        let r = (1.0 - h * h).sqrt();
        for t in [0.0, 0.3, 0.6, 0.9] {
            let p = curve.eval_point(t);
            assert!((p.x.hypot(p.y) - r).abs() < 1e-3, "h = {h}: {p:?}");
        }
    }
}

#[test]
fn test_marching_cylinder_offset_sphere() {
    use crate::{Cylinder, Sphere};
    let tol = Tolerance::DEFAULT;
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let sphere = Sphere::new(Point::new(0.8, 0.0, 0.0), 1.78).unwrap();
    let curves = marched_curves(&cylinder, &sphere, tol);
    assert_one_curve_on_both(&curves, &cylinder, &sphere);
}

#[test]
fn test_marching_shallow_angle_spheres() {
    use crate::Sphere;
    let tol = Tolerance::DEFAULT;
    // Unit spheres whose centers are `d` apart meet at about `d`
    // radians, on the circle z = d/2.
    let a = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    for d in [0.05, 0.0175] {
        let b = Sphere::new(Point::new(0.0, 0.0, d), 1.0).unwrap();
        let curves = marched_curves(&a, &b, tol);
        let curve = assert_one_curve_on_both(&curves, &a, &b);
        let r = (1.0 - d * d / 4.0).sqrt();
        for t in [0.0, 0.3, 0.6, 0.9] {
            let p = curve.eval_point(t);
            assert!((p.z - d / 2.0).abs() < 1e-3, "d = {d}: {p:?}");
            assert!((p.x.hypot(p.y) - r).abs() < 1e-3, "d = {d}: {p:?}");
        }
    }
}

#[test]
fn test_marching_small_loop() {
    // A loop much smaller than the default step: the step follows the
    // curvature.
    use crate::{Plane, Sphere};
    let tol = Tolerance::DEFAULT;
    let plane = Plane::new(Point::new(0.0, 0.0, 0.003), Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 0.01).unwrap();
    let curves = marched_curves(&plane, &sphere, tol);
    let curve = assert_one_curve_on_both(&curves, &plane, &sphere);
    let r = (1e-4f64 - 9e-6).sqrt();
    for t in [0.0, 0.3, 0.6, 0.9] {
        let p = curve.eval_point(t);
        assert!((p.x.hypot(p.y) - r).abs() < 1e-4, "{p:?}");
    }
}

#[test]
fn test_marching_open_trace_is_not_closed() {
    use crate::Plane;
    let tol = Tolerance::DEFAULT;
    let xy = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let zy = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    let curves = marched_curves(&xy, &zy, tol);
    assert_eq!(curves.len(), 1);
    let (t0, t1) = curves[0].bounds();
    let gap = curves[0].eval_point(t0).distance(curves[0].eval_point(t1));
    assert!(gap > 10.0, "open line fitted as a loop: end gap {gap}");
}

/// Polyline length of a traced curve.
fn traced_length(curve: &BSplineCurve) -> f64 {
    let (t0, t1) = curve.bounds();
    let n = 4000;
    let at = |k: usize| curve.eval_point(t0 + (t1 - t0) * k as f64 / n as f64);
    (0..n).map(|k| at(k).distance(at(k + 1))).sum()
}

/// The torus of radii 3 and 1 about z, and the plane through its
/// center tilted `asin(1/3) + d` about y.
fn villarceau_pair(d: f64) -> (crate::Torus, crate::Plane) {
    let torus = crate::Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let a = (1.0f64 / 3.0).asin() + d;
    let plane = crate::Plane::new(Point::ORIGIN, Vector::new(-a.sin(), 0.0, a.cos())).unwrap();
    (torus, plane)
}

#[test]
fn test_marching_stops_at_a_self_crossing() {
    // Issue #7: at the Villarceau tilt the section is two circles of
    // radius 3 crossing where the plane touches the tube. A trace used
    // to turn onto the other circle there and cover an arc twice.
    let tol = Tolerance::DEFAULT;
    let (torus, plane) = villarceau_pair(0.0);
    let curves = marched_curves(&torus, &plane, tol);
    assert!(!curves.is_empty());
    for curve in &curves {
        let (t0, t1) = curve.bounds();
        for k in 0..=40 {
            let p = curve.eval_point(t0 + (t1 - t0) * k as f64 / 40.0);
            assert!(torus.project_point(p, tol).distance < 1e-3, "{p:?}");
            assert!(plane.project_point(p, tol).distance < 1e-3, "{p:?}");
        }
    }
    let total: f64 = curves.iter().map(traced_length).sum();
    let expected = 4.0 * std::f64::consts::PI * 3.0;
    assert!(
        (total - expected).abs() < 1e-3 * expected,
        "total length {total}, expected {expected}: {:?}",
        curves.iter().map(traced_length).collect::<Vec<_>>()
    );
}

#[test]
fn test_marching_near_the_villarceau_tilt_gives_two_loops() {
    let tol = Tolerance::DEFAULT;
    for d in [-1e-3, 1e-3] {
        let (torus, plane) = villarceau_pair(d);
        let curves = marched_curves(&torus, &plane, tol);
        assert_eq!(curves.len(), 2, "d = {d}");
        for curve in &curves {
            let (t0, t1) = curve.bounds();
            let gap = curve.eval_point(t0).distance(curve.eval_point(t1));
            assert!(gap < 1e-6, "d = {d}: open curve, end gap {gap}");
        }
    }
}

#[test]
fn test_segment_meets_box() {
    let unit = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let meets = |p: [f64; 3], q: [f64; 3]| {
        segment_meets_box(
            Point::new(p[0], p[1], p[2]),
            Point::new(q[0], q[1], q[2]),
            &unit,
        )
    };
    assert!(meets([-1.0, 0.5, 0.5], [2.0, 0.5, 0.5]));
    assert!(meets([0.2, 0.3, 0.4], [0.2, 0.3, 0.4]));
    assert!(meets([-1.0, -1.0, 0.5], [0.0, 0.0, 0.5]));
    assert!(!meets([-1.0, 0.5, 0.5], [-0.1, 0.5, 0.5]));
    assert!(!meets([-1.0, 0.0, 0.5], [0.0, -1.0, 0.5]));
    assert!(!meets([0.5, 0.5, 2.0], [0.5, 0.6, 2.0]));
}

/// Asserts every traced curve is closed and lies on both surfaces.
fn assert_closed_on_both<A, B>(curves: &[BSplineCurve], a: &A, b: &B)
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let tol = Tolerance::DEFAULT;
    for curve in curves {
        let (t0, t1) = curve.bounds();
        let gap = curve.eval_point(t0).distance(curve.eval_point(t1));
        assert!(gap < 1e-6, "open curve, end gap {gap}");
        for k in 0..=40 {
            let p = curve.eval_point(t0 + (t1 - t0) * k as f64 / 40.0);
            assert!(a.project_point(p, tol).distance < 1e-3, "{p:?}");
            assert!(b.project_point(p, tol).distance < 1e-3, "{p:?}");
        }
    }
}

#[test]
fn test_marching_finds_skew_cylinders_far_from_the_origin() {
    // The cylinders meet around x = 500, far beyond the grid window of
    // +-10 along the first cylinder's axis, in one loop about 510 long.
    use crate::Cylinder;
    let tol = Tolerance::DEFAULT;
    let along_x = Cylinder::new(Point::ORIGIN, Vector::X, 50.0).unwrap();
    let along_y = Cylinder::new(Point::new(500.0, 0.0, 30.0), Vector::Y, 50.0).unwrap();
    // The same second cylinder, its parameters measured from y = 300:
    // neither grid window holds a point of the intersection.
    let along_y_far = Cylinder::new(Point::new(500.0, 300.0, 30.0), Vector::Y, 50.0).unwrap();
    for curves in [
        marched_curves(&along_x, &along_y, tol),
        marched_curves(&along_y, &along_x, tol),
        marched_curves(&along_x, &along_y_far, tol),
    ] {
        assert!((1..=2).contains(&curves.len()), "{} curves", curves.len());
        assert_closed_on_both(&curves, &along_x, &along_y);
        for curve in &curves {
            let (t0, t1) = curve.bounds();
            for k in 0..=20 {
                let p = curve.eval_point(t0 + (t1 - t0) * k as f64 / 20.0);
                assert!((p.x - 500.0).abs() <= 50.0 + 1e-6, "{p:?}");
            }
        }
    }
}

#[test]
fn test_marching_steinmetz_cylinders_trace_both_ellipses() {
    // Equal radii, perpendicular axes crossing: two ellipses of
    // semi-axes 1 and sqrt(2), crossing where the cylinders touch, at
    // (0, 0, +-1). Traces stop at a crossing (no marching through
    // tangent points), so each ellipse comes as its two halves.
    use crate::Cylinder;
    let tol = Tolerance::DEFAULT;
    let along_x = Cylinder::new(Point::ORIGIN, Vector::X, 1.0).unwrap();
    let along_y = Cylinder::new(Point::ORIGIN, Vector::Y, 1.0).unwrap();
    let curves = marched_curves(&along_x, &along_y, tol);
    assert_eq!(curves.len(), 4);
    for curve in &curves {
        let (t0, t1) = curve.bounds();
        for end in [t0, t1] {
            let p = curve.eval_point(end);
            assert!((p.z.abs() - 1.0).abs() < 1e-6, "{p:?}");
        }
        for k in 0..=40 {
            let p = curve.eval_point(t0 + (t1 - t0) * k as f64 / 40.0);
            assert!(along_x.project_point(p, tol).distance < 1e-3, "{p:?}");
            assert!(along_y.project_point(p, tol).distance < 1e-3, "{p:?}");
        }
    }
    // Two ellipse perimeters (Ramanujan's approximation).
    let (a, b) = (2f64.sqrt(), 1.0);
    let h = ((a - b) / (a + b)).powi(2);
    let perimeter =
        std::f64::consts::PI * (a + b) * (1.0 + 3.0 * h / (10.0 + (4.0 - 3.0 * h).sqrt()));
    let total: f64 = curves.iter().map(traced_length).sum();
    assert!(
        (total - 2.0 * perimeter).abs() < 1e-3 * perimeter,
        "{total}"
    );
}

#[test]
fn test_marching_oblique_torus_plane() {
    use crate::{Plane, Torus};
    let tol = Tolerance::DEFAULT;
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let tilt: f64 = 0.6;
    let plane = Plane::new(
        Point::new(0.5, 0.0, 0.2),
        Vector::new(-tilt.sin(), 0.0, tilt.cos()),
    )
    .unwrap();
    let curves = marched_curves(&torus, &plane, tol);
    assert_eq!(curves.len(), 2);
    assert_closed_on_both(&curves, &torus, &plane);
}

/// Asserts every traced curve lies on both surfaces.
fn assert_on_both<A, B>(curves: &[BSplineCurve], a: &A, b: &B)
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let tol = Tolerance::DEFAULT;
    for curve in curves {
        let (t0, t1) = curve.bounds();
        for k in 0..=40 {
            let p = curve.eval_point(t0 + (t1 - t0) * k as f64 / 40.0);
            assert!(a.project_point(p, tol).distance < 1e-3, "{p:?}");
            assert!(b.project_point(p, tol).distance < 1e-3, "{p:?}");
        }
    }
}

/// A cubic wall waving between x = 0.1 and x = 0.4 over y in
/// `[-2.5, 2.5]`, extruded along z over `[-3, 3]`.
fn extruded_wall() -> crate::surfaces::BSplineSurface {
    let poles = (0..6)
        .map(|i| {
            let y = -2.5 + i as f64;
            let x = if i % 2 == 0 { 0.4 } else { 0.1 };
            vec![Point::new(x, y, -3.0), Point::new(x, y, 3.0)]
        })
        .collect();
    crate::surfaces::BSplineSurface::new(
        3,
        1,
        poles,
        vec![0.0, 1.0, 2.0, 3.0],
        vec![4, 1, 1, 4],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_cylinder_along_an_extruded_wall_keeps_the_wall_seeds() {
    // Every ruling runs along (or nearly along) the extrusion, so the
    // wall's own grid seeds find the two curves.
    use crate::{Cylinder, Intersect};
    let tol = Tolerance::DEFAULT;
    let wall = extruded_wall();
    let tilt = 3f64.to_radians();
    for axis in [Vector::Z, Vector::new(0.0, tilt.sin(), tilt.cos())] {
        let cylinder = Cylinder::new(Point::ORIGIN, axis, 1.0).unwrap();
        let result = cylinder.intersect(&wall, tol);
        let curves: Vec<BSplineCurve> = result
            .curves()
            .iter()
            .map(|c| c.curve.as_bspline().unwrap().clone())
            .collect();
        assert_eq!(curves.len(), 2, "axis {axis:?}");
        assert_on_both(&curves, &cylinder, &wall);
    }
}

#[test]
fn test_marching_parallel_cylinders_seed_from_circles() {
    // Every ruling is parallel to the other axis; a circle meets both
    // lines of the intersection.
    use crate::Cylinder;
    let tol = Tolerance::DEFAULT;
    let first = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let second = Cylinder::new(Point::new(0.5, 0.0, 0.0), Vector::Z, 1.0).unwrap();
    let curves = marched_curves(&first, &second, tol);
    assert_eq!(curves.len(), 2);
    assert_on_both(&curves, &first, &second);
}

#[test]
fn test_marching_sphere_between_torus_meridians() {
    // Spheres centred halfway between two seeding meridians of the
    // torus: the meridians on either side miss them, or lie wholly
    // inside the large one, so the grids find the curves.
    use crate::{Sphere, Torus};
    let tol = Tolerance::DEFAULT;
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let azimuth = 7.5f64.to_radians();
    for (distance, radius, count) in [(3.0, 1.2, 2), (4.0, 0.3, 1), (4.0, 0.15, 1), (2.0, 0.25, 1)]
    {
        let centre = Point::new(distance * azimuth.cos(), distance * azimuth.sin(), 0.0);
        let sphere = Sphere::new(centre, radius).unwrap();
        let curves = marched_curves(&sphere, &torus, tol);
        assert_eq!(curves.len(), count, "radius {radius} at {distance}");
        assert_closed_on_both(&curves, &sphere, &torus);
    }
}

#[test]
fn test_marching_against_a_torus_finds_what_the_grid_finds() {
    // Random spheres, cylinders and cones about a torus: wherever the
    // grid seeds find curves, so does the seeding.
    use crate::{Cone, Cylinder, Sphere, Torus};
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let mut state = 42u64;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    fn check<S: ParametricSurface>(surface: &S, torus: &crate::Torus, case: usize) {
        let tol = Tolerance::DEFAULT;
        if !grid_marching(surface, torus, tol).is_empty() {
            assert!(
                !marched_curves(surface, torus, tol).is_empty(),
                "case {case}"
            );
        }
    }
    for case in 0..24 {
        let azimuth = next() * std::f64::consts::TAU;
        let distance = 1.5 + 3.0 * next();
        let height = -1.5 + 3.0 * next();
        let centre = Point::new(distance * azimuth.cos(), distance * azimuth.sin(), height);
        let radius = 0.1 + 1.5 * next();
        let axis = Vector::new(next() - 0.5, next() - 0.5, next() - 0.5);
        check(&Sphere::new(centre, radius).unwrap(), &torus, case);
        check(
            &Cylinder::new(centre, axis, 0.1 + 0.6 * next()).unwrap(),
            &torus,
            case,
        );
        let (semi_angle, reference) = (0.1 + 0.5 * next(), 0.1 + 0.6 * next());
        check(
            &Cone::new(centre, axis, semi_angle, reference).unwrap(),
            &torus,
            case,
        );
    }
}

/// The straight wall along the diagonal of the xy plane over 280,
/// extruded along z over `[-1, 1]`, with one pole pushed off by `bump`
/// at `at` along it, its neighbours `width` apart, with the given knots
/// (uniform when `None`).
fn bumped_wall(
    bump: f64,
    width: f64,
    at: f64,
    knots: Option<&[f64]>,
) -> crate::surfaces::BSplineSurface {
    let r = std::f64::consts::FRAC_1_SQRT_2;
    let along = [
        0.0,
        40.0,
        80.0,
        120.0,
        at - 2.0 * width,
        at - width,
        at,
        at + width,
        at + 2.0 * width,
        200.0,
        240.0,
        280.0,
    ];
    let poles = along
        .iter()
        .enumerate()
        .map(|(i, &s)| {
            let off = if i == 6 { bump } else { 0.0 };
            let (x, y) = ((s + off) * r, (s - off) * r);
            vec![Point::new(x, y, -1.0), Point::new(x, y, 1.0)]
        })
        .collect::<Vec<_>>();
    let knots: Vec<f64> = knots.map_or_else(
        || (0..=poles.len() - 3).map(|k| k as f64).collect(),
        <[f64]>::to_vec,
    );
    let mut multiplicities = vec![1; knots.len()];
    multiplicities[0] = 4;
    *multiplicities.last_mut().unwrap() = 4;
    crate::surfaces::BSplineSurface::new(
        3,
        1,
        poles,
        knots,
        multiplicities,
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_marching_does_not_step_over_a_small_bump() {
    // A straight run with one small bump whose ends have the same
    // tangent: a long step across it would see no turn.
    use crate::Plane;
    let tol = Tolerance::DEFAULT;
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    // Chord-length knots, and a wide span before narrow ones at the
    // bump: a step must not cross the narrow spans in one go.
    let chord_length = [
        0.0, 40.0, 80.0, 120.0, 149.0, 149.5, 150.0, 150.5, 151.0, 280.0,
    ];
    let wide_then_narrow = [0.0, 1.0, 2.0, 3.0, 3.001, 3.002, 3.003, 3.004, 3.005, 4.0];
    for (bump, width, at, knots) in [
        (0.05, 0.25, 125.0, None),
        (0.02, 0.1, 150.0, None),
        (0.2, 0.5, 175.0, None),
        (0.02, 0.1, 150.0, Some(&chord_length[..])),
        (0.05, 0.25, 150.0, Some(&wide_then_narrow[..])),
        (0.2, 0.5, 150.0, Some(&wide_then_narrow[..])),
    ] {
        let wall = bumped_wall(bump, width, at, knots);
        let curves = marched_curves(&plane, &wall, tol);
        assert_eq!(curves.len(), 1, "bump {bump} at {at}");
        let curve = &curves[0];
        let knots = curve.knots();
        for span in knots.windows(2) {
            for k in 0..20 {
                let p = curve.eval_point(span[0] + (span[1] - span[0]) * k as f64 / 20.0);
                assert!(
                    wall.project_point(p, tol).distance < 1e-3,
                    "bump {bump}: {p:?}"
                );
                assert!(p.z.abs() < 1e-3, "{p:?}");
            }
        }
    }
}

/// A wall of degree 1 along x over `[0, 100]`, extruded along z over
/// `[-1, 1]`, at unit speed except over a tent of height `bump` on
/// `[at, at + 2 * half]`, whose three knots are `cluster` apart: the
/// speed jumps there by a factor `half / cluster`.
fn tent_wall(bump: f64, half: f64, at: f64, cluster: f64) -> crate::surfaces::BSplineSurface {
    let along = [0.0, at, at + half, at + 2.0 * half, 100.0];
    let poles = along
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let y = if i == 2 { bump } else { 0.0 };
            vec![Point::new(x, y, -1.0), Point::new(x, y, 1.0)]
        })
        .collect::<Vec<_>>();
    let knots = vec![
        0.0,
        at,
        at + 0.5 * cluster,
        at + cluster,
        100.0 - 2.0 * half + cluster,
    ];
    crate::surfaces::BSplineSurface::new(
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
    .unwrap()
}

#[test]
fn test_marching_does_not_skip_a_tent_in_a_knot_cluster() {
    // The tent's knots are so close in parameter that their spans look
    // tiny at the step's average speed, while they are 0.1 to 1 long
    // in space: a step must not cross them whole. The curve breaks at
    // the tent's three corners into four pieces, each on both surfaces
    // and between consecutive corners, so none is skipped or traced
    // twice.
    use crate::Plane;
    let tol = Tolerance::DEFAULT;
    let at = 50.0;
    for (cluster, bump, half, z0) in [
        (5e-8, 0.05, 0.05, 0.0),
        (5e-8, 0.2, 0.1, 0.0),
        (5e-8, 0.5, 0.5, 0.3),
        (1e-8, 0.05, 0.05, 0.3),
    ] {
        let wall = tent_wall(bump, half, at, cluster);
        let plane = Plane::new(Point::new(0.0, 0.0, z0), Vector::Z).unwrap();
        let curves = marched_curves(&plane, &wall, tol);
        assert_eq!(curves.len(), 4, "cluster {cluster:e} bump {bump}");
        // The distance to the wall, exactly: to the polyline through
        // its corners in the xy plane (its projection misreports inside
        // the tiny knot spans).
        let corners = [
            (0.0, 0.0),
            (at, 0.0),
            (at + half, bump),
            (at + 2.0 * half, 0.0),
            (100.0, 0.0),
        ];
        let wall_distance = |p: Point| {
            corners
                .windows(2)
                .map(|side| {
                    let ((ax, ay), (bx, by)) = (side[0], side[1]);
                    let (dx, dy) = (bx - ax, by - ay);
                    let along =
                        (((p.x - ax) * dx + (p.y - ay) * dy) / (dx * dx + dy * dy)).clamp(0.0, 1.0);
                    (p.x - ax - along * dx).hypot(p.y - ay - along * dy)
                })
                .fold(f64::INFINITY, f64::min)
        };
        let mut extents = Vec::new();
        for curve in &curves {
            let (mut lowest, mut highest) = (f64::INFINITY, f64::NEG_INFINITY);
            let knots = curve.knots();
            for span in knots.windows(2) {
                for k in 0..=20 {
                    let p = curve.eval_point(span[0] + (span[1] - span[0]) * k as f64 / 20.0);
                    lowest = lowest.min(p.x);
                    highest = highest.max(p.x);
                    assert!(
                        wall_distance(p) < 1e-3,
                        "cluster {cluster:e} bump {bump}: {p:?}"
                    );
                    assert!((p.z - z0).abs() < 1e-3, "{p:?}");
                }
            }
            extents.push((lowest, highest));
        }
        extents.sort_by(|a, b| a.0.total_cmp(&b.0));
        for ((lowest, highest), ends) in extents.iter().zip(corners.windows(2)) {
            assert!(
                (lowest - ends[0].0).abs() < 1e-3 && (highest - ends[1].0).abs() < 1e-3,
                "cluster {cluster:e} bump {bump}: piece [{lowest}, {highest}] for {ends:?}"
            );
        }
    }
}

#[test]
fn test_marching_stops_on_the_model_size_cube() {
    // A cylinder whose axis is parallel to a generator of the cone, off
    // the apex: the intersection runs to infinity.
    use crate::{Cone, Cylinder};
    let tol = Tolerance::DEFAULT;
    let cone = Cone::new(Point::ORIGIN, Vector::Z, 0.3, 1.0).unwrap();
    let generator = Vector::new(0.3f64.sin(), 0.0, 0.3f64.cos());
    let cylinder = Cylinder::new(Point::new(1.0, 0.0, 0.0), generator, 0.5).unwrap();
    let curves = marched_curves(&cone, &cylinder, tol);
    assert!(!curves.is_empty());
    let half = 0.5 * Tolerance::MODEL_SIZE;
    let on_cube = |p: Point| {
        [p.x, p.y, p.z]
            .iter()
            .any(|c| (c.abs() - half).abs() < 1e-3)
    };
    for curve in &curves {
        let (t0, t1) = curve.bounds();
        let ends = [curve.eval_point(t0), curve.eval_point(t1)];
        // Both branches of each curve run to infinity: every trace ends
        // on the cube, so none ran out of steps.
        assert!(ends.iter().all(|&p| on_cube(p)), "{ends:?}");
        assert!(curve.poles().len() < 2000, "{} poles", curve.poles().len());
    }
    assert_on_both(&curves, &cone, &cylinder);
}

#[test]
fn test_marching_finds_bounded_curves_beyond_the_model_size_cube() {
    // Two tori in general position cross in two loops, wherever they lie:
    // a bounded surface bounds the curves, so traces end on its box, not
    // on the cube around the origin (#138).
    use crate::Torus;
    let tol = Tolerance::DEFAULT;
    let pair = |offset: Vector| {
        let first = Torus::new(Point::ORIGIN + offset, Vector::Z, 3.0, 1.0).unwrap();
        let axis = Vector::new(1.0, 0.4, 0.3);
        let second = Torus::new(Point::new(2.0, 0.5, 0.2) + offset, axis, 2.5, 0.8).unwrap();
        (first, second)
    };
    let (first, second) = pair(Vector::ZERO);
    let near = marched_curves(&first, &second, tol);
    assert!(!near.is_empty());
    for offset in [1e5, 1e6] {
        let (first, second) = pair(Vector::new(offset, offset, offset));
        let far = marched_curves(&first, &second, tol);
        assert_eq!(far.len(), near.len(), "offset {offset:e}");
        assert_on_both(&far, &first, &second);
    }
}

#[test]
fn test_marching_miss_is_empty() {
    use crate::{Plane, Sphere};
    let tol = Tolerance::DEFAULT;
    let plane = Plane::new(Point::new(0.0, 0.0, 5.0), Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    assert!(marched_curves(&plane, &sphere, tol).is_empty());
}

fn cylinder_patch(u_range: (f64, f64)) -> crate::surfaces::BSplineSurface {
    crate::Cylinder::new(Point::ORIGIN, Vector::Z, 1.0)
        .unwrap()
        .to_bspline_patch(u_range, (0.0, 2.0))
        .unwrap()
}

#[test]
fn test_surface_and_its_bspline_form_are_coincident() {
    use crate::Sphere;
    let tol = Tolerance::DEFAULT;
    let sphere = Sphere::new(Point::new(0.3, 0.2, 0.5), 0.8).unwrap();
    let copy = sphere.to_bspline();
    for result in [
        marching_intersection(&sphere, &copy, tol),
        marching_intersection(&copy, &sphere, tol),
        marching_intersection(&copy, &copy, tol),
    ] {
        assert!(result.is_coincident() && !result.is_empty());
        assert!(result.curves().is_empty());
        assert_eq!(result.completeness(), &Completeness::Sampled);
    }
}

/// The region of an incomplete result with one region.
fn only_region(result: &SurfaceSurfaceIntersection) -> SurfaceSurfaceRegion {
    match result.completeness() {
        Completeness::Incomplete {
            reason: IncompleteReason::Unseparated,
            regions,
        } if regions.len() == 1 => regions[0],
        other => panic!("expected one unseparated region, got {other:?}"),
    }
}

#[test]
fn test_patch_lying_on_its_surface_is_incomplete_over_the_patch() {
    use crate::{Cylinder, Plane};
    use std::f64::consts::{PI, TAU};
    let tol = Tolerance::DEFAULT;
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let whole = cylinder_patch((0.0, TAU));
    let result = marching_intersection(&cylinder, &whole, tol);
    assert!(!result.is_coincident() && result.curves().is_empty());
    let region = only_region(&result);
    // The patch spans the whole turn: the cylinder's range is the period.
    assert_eq!(region.first_u_range, (0.0, TAU));
    assert_eq!(region.first_v_range, (0.0, 2.0));
    assert_eq!(region.second_u_range, whole.u_bounds());
    // Swapped, the patch is the first surface.
    let region = only_region(&marching_intersection(&whole, &cylinder, tol));
    assert_eq!(region.first_u_range, whole.u_bounds());
    assert_eq!(region.second_u_range, (0.0, TAU));
    // A half turn keeps to its range on the cylinder.
    let half = cylinder_patch((0.0, PI));
    let region = only_region(&marching_intersection(&cylinder, &half, tol));
    assert!(region.first_u_range.0.abs() < 1e-9 && (region.first_u_range.1 - PI).abs() < 1e-9);
    // A plane against its own bilinear patch.
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let patch = plane.to_bspline_patch((-1.0, 2.0), (0.0, 3.0)).unwrap();
    let region = only_region(&marching_intersection(&plane, &patch, tol));
    assert_eq!(
        (region.first_u_range, region.first_v_range),
        ((-1.0, 2.0), (0.0, 3.0))
    );
    assert_eq!(region.min_corner, Point::new(-1.0, 0.0, 0.0));
    assert_eq!(region.max_corner, Point::new(2.0, 3.0, 0.0));
}

#[test]
fn test_crossing_surfaces_are_not_taken_for_coincident() {
    use crate::Plane;
    let tol = Tolerance::DEFAULT;
    // Crossing the patch through its first sample row: the check stops
    // at the first sample off the plane.
    let plane = Plane::new(Point::new(1.0, 0.0, 0.0), Vector::new(1.0, 0.0, 1.0)).unwrap();
    let result = marching_intersection(&plane, &cylinder_patch((0.0, std::f64::consts::TAU)), tol);
    assert!(!result.is_coincident());
    assert_eq!(result.completeness(), &Completeness::Sampled);
    assert!(!result.curves().is_empty());
}

#[test]
fn test_plane_touching_a_cylinder_patch_traces_a_tangent_ruling() {
    use crate::Plane;
    let tol = Tolerance::DEFAULT;
    let patch = cylinder_patch((0.0, std::f64::consts::TAU));
    // Touching along x = 1, y = 0, and 1e-9 inside it, where the two
    // crossings lie 4.5e-5 apart at an angle of 4.5e-5 radians: a touch
    // within confusion.
    for offset in [1.0, 1.0 - 1e-9] {
        let plane = Plane::new(Point::new(offset, 0.0, 0.0), Vector::X).unwrap();
        let result = marching_intersection(&plane, &patch, tol);
        assert_eq!(result.curves().len(), 1);
        assert!(result.curves()[0].tangent);
        let Curve::BSpline(curve) = &result.curves()[0].curve else {
            panic!("a numeric curve is a B-spline");
        };
        let (t0, t1) = curve.bounds();
        let mut heights = [curve.eval_point(t0).z, curve.eval_point(t1).z];
        heights.sort_by(f64::total_cmp);
        assert!(heights[0].abs() < 1e-6 && (heights[1] - 2.0).abs() < 1e-6);
        for k in 0..=8 {
            let p = curve.eval_point(t0 + (t1 - t0) * k as f64 / 8.0);
            assert!((p.x - offset).abs() < 1e-7 && p.y.abs() < 1e-6, "{p:?}");
        }
    }
    // Crossing at 0.99, the two curves are not tangent.
    let plane = Plane::new(Point::new(0.99, 0.0, 0.0), Vector::X).unwrap();
    let result = marching_intersection(&plane, &patch, tol);
    assert_eq!(result.curves().len(), 2);
    assert!(result.curves().iter().all(|piece| !piece.tangent));
}

#[test]
fn test_sphere_inside_a_cylinder_patch_traces_the_touching_circle() {
    use crate::{Sphere, Torus};
    let tol = Tolerance::DEFAULT;
    let patch = cylinder_patch((0.0, std::f64::consts::TAU));
    let sphere = Sphere::new(Point::new(0.0, 0.0, 1.0), 1.0).unwrap();
    let torus = Torus::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 2.0, 1.0).unwrap();
    for result in [
        marching_intersection(&sphere, &patch, tol),
        marching_intersection(&torus, &patch, tol),
    ] {
        assert_eq!(result.curves().len(), 1);
        assert!(result.curves()[0].tangent);
        let Curve::BSpline(curve) = &result.curves()[0].curve else {
            panic!("a numeric curve is a B-spline");
        };
        let (t0, t1) = curve.bounds();
        assert!(curve.eval_point(t0).distance(curve.eval_point(t1)) < 1e-9);
        for k in 0..=16 {
            let p = curve.eval_point(t0 + (t1 - t0) * k as f64 / 16.0);
            let radius = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (radius - 1.0).abs() < 1e-6 && (p.z - 1.0).abs() < 1e-6,
                "{p:?}"
            );
        }
    }
}

#[test]
fn test_sphere_touching_a_cylinder_patch_at_a_point_is_a_contact_point() {
    use crate::Sphere;
    let tol = Tolerance::DEFAULT;
    let patch = cylinder_patch((0.0, std::f64::consts::TAU));
    let sphere = Sphere::new(Point::new(2.0, 0.0, 1.0), 1.0).unwrap();
    let result = marching_intersection(&sphere, &patch, tol);
    assert!(result.curves().is_empty());
    assert_eq!(result.contact_points().len(), 1);
    assert!(result.contact_points()[0].distance(Point::new(1.0, 0.0, 1.0)) < 1e-6);
    assert!(!result.is_empty());
}

#[test]
fn test_patches_overlapping_in_part_report_a_region_over_the_overlap() {
    let tol = Tolerance::DEFAULT;
    let first = cylinder_patch((0.0, 3.0));
    let second = cylinder_patch((1.5, 4.5));
    let result = marching_intersection(&first, &second, tol);
    assert!(!result.is_coincident());
    let Completeness::Incomplete { reason, regions } = result.completeness() else {
        panic!(
            "expected an incomplete result, got {:?}",
            result.completeness()
        );
    };
    assert_eq!(*reason, IncompleteReason::Unseparated);
    // The overlap is u in [1.5, 3] on both patches.
    assert!(
        regions
            .iter()
            .any(|region| region.first_u_range.0 <= 1.5 + 1e-6
                && region.first_u_range.1 >= 3.0 - 1e-6)
    );
    assert!(result.contact_points().is_empty());
}

/// `surface` moved by `offset`.
fn shifted(
    surface: &crate::surfaces::BSplineSurface,
    offset: Vector,
) -> crate::surfaces::BSplineSurface {
    surface.transformed(crate::Transform::translation(offset))
}

#[test]
fn test_near_copies_within_confusion_are_coincident() {
    use crate::{Cylinder, Sphere};
    let tol = Tolerance::DEFAULT;
    let sphere = Sphere::new(Point::new(0.3, 0.2, 0.5), 0.8).unwrap();
    let copy = sphere.to_bspline();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let patch = cylinder_patch((0.0, std::f64::consts::PI));
    for shift in [3e-8, 1e-7] {
        let offset = Vector::new(1.0, 0.5, 0.25).normalized().unwrap() * shift;
        let moved = shifted(&copy, offset);
        for result in [
            marching_intersection(&sphere, &moved, tol),
            marching_intersection(&moved, &sphere, tol),
            marching_intersection(&copy, &moved, tol),
        ] {
            assert!(result.is_coincident(), "shift {shift}");
            assert!(result.curves().is_empty());
        }
        let result = marching_intersection(&cylinder, &shifted(&patch, offset), tol);
        assert!(!result.is_coincident() && result.curves().is_empty());
        let region = only_region(&result);
        assert!(region.first_u_range.0.abs() < 1e-6);
        assert!((region.first_u_range.1 - std::f64::consts::PI).abs() < 1e-6);
    }
}

#[test]
fn test_copies_apart_beyond_confusion_are_not_coincident() {
    use crate::Sphere;
    let tol = Tolerance::DEFAULT;
    let sphere = Sphere::new(Point::new(0.3, 0.2, 0.5), 0.8).unwrap();
    let copy = sphere.to_bspline();
    // Equal spheres 1e-5 apart cross in one great circle.
    let offset = Vector::new(1.0, 0.5, 0.25).normalized().unwrap() * 1e-5;
    let result = marching_intersection(&sphere, &shifted(&copy, offset), tol);
    assert!(!result.is_coincident());
    assert_eq!(result.completeness(), &Completeness::Sampled);
    assert_eq!(result.curves().len(), 1);
}

/// A flat cubic patch over [0, 1] x [0, 1] in z = 0 with bumps of
/// `height`, alternately up and down, between the lines u = k / 6:
/// each bump is a single basis function of a uniform 24-span knot
/// vector, with zero value and slope on those lines.
fn bumpy_patch(height: f64) -> crate::surfaces::BSplineSurface {
    // A uniform cubic basis function peaks at 2/3 of its pole.
    let pole_height = 1.5 * height;
    let poles = (0..27)
        .map(|i| {
            let x = if i == 0 {
                0.0
            } else if i == 26 {
                1.0
            } else {
                (i as f64 - 1.0) / 24.0
            };
            let z = if i % 4 == 3 && i <= 23 {
                if (i / 4) % 2 == 0 {
                    pole_height
                } else {
                    -pole_height
                }
            } else {
                0.0
            };
            vec![Point::new(x, 0.0, z), Point::new(x, 1.0, z)]
        })
        .collect();
    let knots = (0..=24).map(|k| k as f64 / 24.0).collect();
    let mut multiplicities = vec![1; 25];
    (multiplicities[0], multiplicities[24]) = (4, 4);
    crate::surfaces::BSplineSurface::new(
        3,
        1,
        poles,
        knots,
        multiplicities,
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_bumps_between_the_sample_lines_are_not_coincident() {
    use crate::Plane;
    let tol = Tolerance::DEFAULT;
    let flat = Plane::new(Point::ORIGIN, Vector::Z)
        .unwrap()
        .to_bspline_patch((0.0, 1.0), (0.0, 1.0))
        .unwrap();
    for height in [0.1, 1e-3] {
        let bumpy = bumpy_patch(height);
        for k in 0..6 {
            let peak = bumpy.eval_point((k as f64 + 0.5) / 6.0, 0.5);
            assert!((peak.z.abs() - height).abs() < 1e-9, "{peak:?}");
        }
        for result in [
            marching_intersection(&flat, &bumpy, tol),
            marching_intersection(&bumpy, &flat, tol),
        ] {
            assert!(!result.is_coincident());
            // The tall bumps touch the flat patch along the seven
            // lines; the low ones lie flat enough along them that
            // marching may leave part unseparated, but never takes
            // the two for one surface.
            if height == 0.1 {
                assert_eq!(result.curves().len(), 7);
                assert_eq!(result.completeness(), &Completeness::Sampled);
            }
        }
    }
}
