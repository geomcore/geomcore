use super::*;
use crate::{Cylinder, Plane, Sphere, Torus};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The parameters of `point` on both surfaces of `pair`, converged onto
/// both from their projections.
fn seed_at<A: ParametricSurface, B: ParametricSurface>(
    pair: &SurfacePair<A, B>,
    point: Point,
) -> [f64; 4] {
    let (pa, pb) = (
        pair.first.project_point(point, TOL),
        pair.second.project_point(point, TOL),
    );
    pair.meet([pa.u, pa.v, pb.u, pb.v], None, TOL)
        .unwrap_or([pa.u, pa.v, pb.u, pb.v])
}

#[test]
fn test_a_crossing_closes_its_loop_on_the_seed() {
    let plane = Plane::new(Point::new(0.0, 0.0, 0.5), Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let pair = SurfacePair::new(&plane, &sphere);
    let seed = seed_at(&pair, Point::new(0.75f64.sqrt(), 0.0, 0.5));
    let Trace::Curve {
        points,
        tangents,
        stalls,
        touching,
    } = Tracer::new(pair, TOL).trace(seed)
    else {
        panic!("no curve traced");
    };
    assert!(!touching);
    assert!(stalls.is_empty());
    assert_eq!(points.first(), points.last());
    // One unit tangent per point, along the way the points run.
    assert_eq!(tangents.len(), points.len());
    for (k, t) in tangents.iter().enumerate().skip(1) {
        assert!((t.magnitude() - 1.0).abs() < 1e-9, "{t:?}");
        assert!(t.dot(points[k] - points[k - 1]) > 0.0);
    }
    for p in &points {
        assert!((p.z - 0.5).abs() <= TOL.confusion, "{p:?}");
        assert!(
            (p.x.hypot(p.y) - 0.75f64.sqrt()).abs() <= TOL.confusion,
            "{p:?}"
        );
    }
    // The steps follow the curvature: a chord deviation of
    // `CHORD_DEVIATION * confusion` on a circle of radius `r` takes steps
    // of `sqrt(8 r deviation)`.
    let radius = 0.75f64.sqrt();
    let step = (8.0 * radius * CHORD_DEVIATION * TOL.confusion).sqrt();
    let steps = std::f64::consts::TAU * radius / step;
    assert!(
        (points.len() as f64 - steps).abs() < 0.1 * steps,
        "{}",
        points.len()
    );
}

#[test]
fn test_a_trace_ends_on_the_parameter_bounds() {
    // A unit square of the plane `z = 0` against the plane `x = 0.5`: the
    // segment from `y = 0` to `y = 1`.
    let square = Plane::new(Point::ORIGIN, Vector::Z)
        .unwrap()
        .to_bspline_patch((0.0, 1.0), (0.0, 1.0))
        .unwrap();
    let wall = Plane::new(Point::new(0.5, 0.0, 0.0), Vector::X).unwrap();
    let pair = SurfacePair::new(&square, &wall);
    let seed = seed_at(&pair, Point::new(0.5, 0.5, 0.0));
    let Trace::Curve {
        points, touching, ..
    } = Tracer::new(pair, TOL).trace(seed)
    else {
        panic!("no curve traced");
    };
    assert!(!touching);
    let ends = [points[0], points[points.len() - 1]].map(|p| p.y);
    let (low, high) = (ends[0].min(ends[1]), ends[0].max(ends[1]));
    assert!(
        low.abs() <= TOL.confusion && (high - 1.0).abs() <= TOL.confusion,
        "{ends:?}"
    );
    for p in &points {
        assert!((p.x - 0.5).abs() <= TOL.confusion && p.z.abs() <= TOL.confusion);
    }
}

#[test]
fn test_a_touch_along_unbounded_surfaces_is_traced_to_the_model_cube() {
    let plane = Plane::new(Point::new(1.0, 0.0, 0.0), Vector::X).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let pair = SurfacePair::new(&plane, &cylinder);
    let seed = seed_at(&pair, Point::new(1.0, 0.0, 0.0));
    let Trace::Curve {
        points, touching, ..
    } = Tracer::new(pair, TOL).trace(seed)
    else {
        panic!("no curve traced");
    };
    assert!(touching);
    for p in &points {
        assert!(
            (p.x - 1.0).abs() <= TOL.confusion && p.y.abs() < 1e-3,
            "{p:?}"
        );
    }
    let ends = [points[0], points[points.len() - 1]].map(|p| p.z.abs());
    for end in ends {
        assert!((end - 0.5 * Tolerance::MODEL_SIZE).abs() < 1e-6, "{ends:?}");
    }
}

#[test]
fn test_a_touch_at_a_point_is_a_touch() {
    // The plane `x = 3` touches the outer equator of the torus of radii 2
    // and 1 at `(3, 0, 0)` only.
    let plane = Plane::new(Point::new(3.0, 0.0, 0.0), Vector::X).unwrap();
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 2.0, 1.0).unwrap();
    let pair = SurfacePair::new(&plane, &torus);
    let seed = seed_at(&pair, Point::new(3.0, 0.0, 0.0));
    assert!(matches!(Tracer::new(pair, TOL).trace(seed), Trace::Touch));
}

#[test]
fn test_a_crossing_without_a_direction_is_lost() {
    // The plane `z = 1` touches the unit sphere at its pole, where the
    // sphere's parametrization is singular: no touch is read there, and
    // the normals, both along z, give no direction to cross along.
    let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let pair = SurfacePair::new(&plane, &sphere);
    let seed = seed_at(&pair, Point::new(0.0, 0.0, 1.0));
    assert!(matches!(Tracer::new(pair, TOL).trace(seed), Trace::Lost));
}

#[test]
fn test_knot_counts_run_on_across_periods() {
    let knots = [0.0, 1.0, 2.0, 4.0];
    assert_eq!(knots_below(&knots, false, 1.0, false), 1);
    assert_eq!(knots_below(&knots, false, 1.0, true), 2);
    // Periodic: three knots per period of 4.
    assert_eq!(knots_below(&knots, true, 5.0, false), 4);
    assert_eq!(knots_below(&knots, true, -1.0, true), 0);
    assert_eq!(knots_below(&knots, true, -2.0, false), -1);
    assert_eq!(knot(&knots, true, 4), 5.0);
    assert_eq!(knot(&knots, true, -1), -2.0);
    assert_eq!(knot(&knots, false, -1), f64::NEG_INFINITY);
    assert_eq!(knot(&knots, false, 4), f64::INFINITY);
}
