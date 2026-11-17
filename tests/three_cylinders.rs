//! Three cylinders of equal radius along the X, Y and Z axes, intersected
//! by every route that reaches their common points.
//!
//! Each pair of the cylinders meets in two ellipses crossing at two
//! singular points. The curves of one pair cross those of another, and
//! cross the third cylinder, at the eight points `(±r/√2, ±r/√2, ±r/√2)`
//! on all three cylinders. The curve–curve route (the curves of two pairs
//! against each other) and the curve–surface route (the curves of a pair
//! against the third cylinder) must find the same eight points: exactly,
//! moved and scaled, on finite faces and on the cylinders as B-spline
//! patches, and, where the answer has no closed form, with one cylinder a
//! little wider, tilted or off the common centre.

use geomcore::prelude::*;
use geomcore::{
    Axis, Completeness, Curve, CurveCurveIntersection, CurveSurfaceIntersection, Cylinder,
    IntersectionMethod, Point, Surface, SurfaceSurfaceIntersection, Tolerance, Transform,
    TrimmedCurve, TrimmedSurface, Vector,
};
use std::f64::consts::{FRAC_1_SQRT_2, TAU};

/// The pairs of cylinders by their places in `[x, y, z]`, each with the
/// place of the third cylinder.
const PAIRS: [(usize, usize, usize); 3] = [(0, 1, 2), (1, 2, 0), (0, 2, 1)];

/// The pairings of pairs, by their places in `PAIRS`, whose curves the
/// curve–curve route crosses.
const PAIRINGS: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];

const PAIR_NAMES: [&str; 3] = ["xy", "yz", "xz"];
const AXIS_NAMES: [&str; 3] = ["x", "y", "z"];

/// The height of the finite cylinders: `v` in `[-HALF_HEIGHT, HALF_HEIGHT]`.
const HALF_HEIGHT: f64 = 2.0;

/// The cylinders along the X, Y and Z axes through the origin, of radius 1.
fn unit_cylinders() -> [Cylinder; 3] {
    [Vector::X, Vector::Y, Vector::Z]
        .map(|direction| Cylinder::new(Point::ORIGIN, direction, 1.0).unwrap())
}

/// The eight points on all three unit cylinders.
fn unit_common_points() -> Vec<Point> {
    let mut points = Vec::new();
    for x in [-1.0, 1.0] {
        for y in [-1.0, 1.0] {
            for z in [-1.0, 1.0] {
                points.push(Point::new(
                    x * FRAC_1_SQRT_2,
                    y * FRAC_1_SQRT_2,
                    z * FRAC_1_SQRT_2,
                ));
            }
        }
    }
    points
}

/// The true distance from `point` to `cylinder`, from its axis: not
/// snapped to zero within the confusion as a projection is.
fn distance_to_cylinder(point: Point, cylinder: &Cylinder) -> f64 {
    (cylinder.axis().distance(point) - cylinder.radius()).abs()
}

/// The largest distance from `point` to the three cylinders.
fn distance_to_all(point: Point, cylinders: &[Cylinder; 3]) -> f64 {
    cylinders
        .iter()
        .map(|cylinder| distance_to_cylinder(point, cylinder))
        .fold(0.0, f64::max)
}

/// The confusion of `Tolerance::DEFAULT` scaled by `scale`.
fn scaled_tolerance(scale: f64) -> Tolerance {
    let reference = Tolerance::DEFAULT;
    Tolerance::new(
        reference.confusion * scale,
        reference.angular,
        reference.parametric,
    )
    .unwrap()
}

/// The intersection of each pair of `surfaces`, in the order of `PAIRS`.
fn pair_intersections<S>(surfaces: &[S; 3], tol: Tolerance) -> [SurfaceSurfaceIntersection; 3]
where
    S: Intersect<S, Output = SurfaceSurfaceIntersection>,
{
    PAIRS.map(|(first, second, _)| surfaces[first].intersect(&surfaces[second], tol))
}

/// The curves of `result`, each trimmed to its parameter range.
fn trimmed_curves(result: &SurfaceSurfaceIntersection) -> Vec<TrimmedCurve> {
    result
        .curves()
        .iter()
        .map(|piece| {
            let (start, end) = piece.parameter_range;
            TrimmedCurve::new(piece.curve.clone(), start, end).unwrap()
        })
        .collect()
}

/// The points a route found, every hit kept: a point found twice is there
/// twice.
struct Route {
    name: String,
    points: Vec<Point>,
}

/// The six routes to the common points: the curves of each pairing of
/// pairs against each other, and the curves of each pair against the
/// third surface.
fn routes<S>(
    surfaces: &[S; 3],
    pairs: &[SurfaceSurfaceIntersection; 3],
    tol: Tolerance,
) -> Vec<Route>
where
    TrimmedCurve: Intersect<S, Output = CurveSurfaceIntersection>
        + Intersect<TrimmedCurve, Output = CurveCurveIntersection>,
{
    let curves = pairs.each_ref().map(trimmed_curves);
    let mut routes = Vec::new();
    for (first, second) in PAIRINGS {
        let mut points = Vec::new();
        for a in &curves[first] {
            for b in &curves[second] {
                points.extend(a.intersect(b, tol).hits().iter().map(|hit| hit.point));
            }
        }
        routes.push(Route {
            name: format!(
                "{} curves x {} curves",
                PAIR_NAMES[first], PAIR_NAMES[second]
            ),
            points,
        });
    }
    for (pair, (_, _, third)) in PAIRS.into_iter().enumerate() {
        let mut points = Vec::new();
        for curve in &curves[pair] {
            let result = curve.intersect(&surfaces[third], tol);
            points.extend(result.hits().iter().map(|hit| hit.point));
        }
        routes.push(Route {
            name: format!(
                "{} curves x {} surface",
                PAIR_NAMES[pair], AXIS_NAMES[third]
            ),
            points,
        });
    }
    routes
}

/// Every curve of every pair, sampled at 50 parameters across its range,
/// lies within `limit` of both its cylinders.
fn assert_curves_on_their_cylinders(
    pairs: &[SurfaceSurfaceIntersection; 3],
    cylinders: &[Cylinder; 3],
    limit: f64,
) {
    for (pair, (first, second, _)) in PAIRS.into_iter().enumerate() {
        for piece in pairs[pair].curves() {
            let (start, end) = piece.parameter_range;
            for i in 0..50 {
                let point = piece
                    .curve
                    .point_at(start + (end - start) * f64::from(i) / 49.0);
                for place in [first, second] {
                    let off = distance_to_cylinder(point, &cylinders[place]);
                    assert!(
                        off <= limit,
                        "a {} curve is {off:e} off the {} cylinder (limit {limit:e})",
                        PAIR_NAMES[pair],
                        AXIS_NAMES[place]
                    );
                }
            }
        }
    }
}

/// Each pair of exact cylinders meets in two ellipses and two contact
/// points, the singular points where the ellipses cross.
fn assert_pairs_meet_in_two_ellipses(pairs: &[SurfaceSurfaceIntersection; 3]) {
    for (pair, result) in pairs.iter().enumerate() {
        let name = PAIR_NAMES[pair];
        assert_eq!(result.method(), IntersectionMethod::Exact, "{name}");
        assert_eq!(result.completeness(), &Completeness::Complete, "{name}");
        assert_eq!(result.curves().len(), 2, "{name}");
        assert_eq!(result.contact_points().len(), 2, "{name}");
        for piece in result.curves() {
            assert!(
                matches!(piece.curve, Curve::Ellipse(_)),
                "{name}: {:?}",
                piece.curve
            );
        }
    }
}

/// How far a pair's curve may stray from its cylinders, a point found from
/// its expected place, and a point found from the three cylinders.
struct Limits {
    curve_from_cylinders: f64,
    point_from_expected: f64,
    point_from_cylinders: f64,
}

/// Every route finds each of `expected` exactly once and nothing else;
/// every point found lies near all three cylinders.
fn assert_routes_find(
    routes: &[Route],
    expected: &[Point],
    cylinders: &[Cylinder; 3],
    limits: &Limits,
) {
    let limit = limits.point_from_expected;
    for route in routes {
        assert_eq!(
            route.points.len(),
            expected.len(),
            "{}: {:?}",
            route.name,
            route.points
        );
        for target in expected {
            let near = route
                .points
                .iter()
                .filter(|point| point.distance(*target) <= limit)
                .count();
            assert_eq!(
                near, 1,
                "{}: {near} points within {limit:e} of {target:?}",
                route.name
            );
        }
        for point in &route.points {
            let off = distance_to_all(*point, cylinders);
            assert!(
                off <= limits.point_from_cylinders,
                "{}: {point:?} is {off:e} off a cylinder (limit {:e})",
                route.name,
                limits.point_from_cylinders
            );
        }
    }
}

/// The pairs of `surfaces`, checked by `check_pairs`, keep their curves
/// on the `cylinders`, and every route finds `expected`.
fn check_known_points<S>(
    surfaces: &[S; 3],
    cylinders: &[Cylinder; 3],
    expected: &[Point],
    tol: Tolerance,
    limits: &Limits,
    check_pairs: impl Fn(&[SurfaceSurfaceIntersection; 3]),
) where
    S: Intersect<S, Output = SurfaceSurfaceIntersection>,
    TrimmedCurve: Intersect<S, Output = CurveSurfaceIntersection>
        + Intersect<TrimmedCurve, Output = CurveCurveIntersection>,
{
    let pairs = pair_intersections(surfaces, tol);
    check_pairs(&pairs);
    assert_curves_on_their_cylinders(&pairs, cylinders, limits.curve_from_cylinders);
    let routes = routes(surfaces, &pairs, tol);
    assert_routes_find(&routes, expected, cylinders, limits);
}

/// The surfaces exactly as the cylinders.
fn whole(cylinders: &[Cylinder; 3]) -> [Surface; 3] {
    cylinders.map(Surface::from)
}

/// The cylinders as finite faces: a whole turn, `v` in
/// `[-HALF_HEIGHT, HALF_HEIGHT]`.
fn faces(cylinders: &[Cylinder; 3]) -> [TrimmedSurface; 3] {
    cylinders.map(|cylinder| {
        Surface::from(cylinder)
            .trimmed((0.0, TAU), (-HALF_HEIGHT, HALF_HEIGHT))
            .unwrap()
    })
}

/// The finite faces as exact rational B-spline patches, intersected by
/// marching.
fn patches(cylinders: &[Cylinder; 3]) -> [Surface; 3] {
    cylinders.map(|cylinder| {
        Surface::from(
            cylinder
                .to_bspline_patch((0.0, TAU), (-HALF_HEIGHT, HALF_HEIGHT))
                .unwrap(),
        )
    })
}

/// How far a point found from a closed form may stray from its expected
/// place, relative to the size of the setup.
const CLOSED_FORM_POINT_LIMIT: f64 = 1e-12;

/// The exact cylinders moved by `transform`: the limits scale with its
/// factor, with room for rounding coordinates as large as the moved
/// points'.
fn check_exact_under(transform: Transform) {
    let scale = transform.scale_factor();
    let cylinders = unit_cylinders().map(|cylinder| cylinder.transformed(transform).unwrap());
    let expected: Vec<Point> = unit_common_points()
        .into_iter()
        .map(|point| transform.apply_point(point))
        .collect();
    let magnitude = expected
        .iter()
        .flat_map(|point| [point.x.abs(), point.y.abs(), point.z.abs()])
        .fold(0.0, f64::max);
    let rounding = 16.0 * magnitude * f64::EPSILON;
    let tol = scaled_tolerance(scale);
    let limits = Limits {
        curve_from_cylinders: tol.confusion + rounding,
        point_from_expected: CLOSED_FORM_POINT_LIMIT * scale + rounding,
        point_from_cylinders: tol.confusion + rounding,
    };
    check_known_points(
        &whole(&cylinders),
        &cylinders,
        &expected,
        tol,
        &limits,
        assert_pairs_meet_in_two_ellipses,
    );
}

#[test]
fn test_unit_cylinders_meet_at_the_eight_points_by_every_route() {
    check_exact_under(Transform::IDENTITY);
}

#[test]
fn test_the_eight_points_follow_the_cylinders_under_a_similarity() {
    let skew_axis = Axis::new(Point::ORIGIN, Vector::new(1.0, 2.0, 3.0)).unwrap();
    let other_axis = Axis::new(Point::ORIGIN, Vector::new(1.0, -1.0, 2.0)).unwrap();
    for transform in [
        Transform::rotation(skew_axis, 0.7),
        Transform::translation(Vector::new(1e3, -2e3, 5e2)),
        Transform::translation(Vector::new(1e7, 3e6, -2e7)),
        Transform::scaling(Point::ORIGIN, 2f64.powi(20)),
        Transform::scaling(Point::ORIGIN, 2f64.powi(-20)),
        Transform::rotation(other_axis, 1.1)
            .then(Transform::translation(Vector::new(1e5, 2e5, 3e5)))
            .then(Transform::scaling(Point::ORIGIN, 1e3)),
    ] {
        check_exact_under(transform);
    }
}

#[test]
fn test_finite_cylinders_meet_at_the_eight_points_by_every_route() {
    let cylinders = unit_cylinders();
    let tol = Tolerance::DEFAULT;
    let limits = Limits {
        curve_from_cylinders: tol.confusion,
        point_from_expected: CLOSED_FORM_POINT_LIMIT,
        point_from_cylinders: tol.confusion,
    };
    check_known_points(
        &faces(&cylinders),
        &cylinders,
        &unit_common_points(),
        tol,
        &limits,
        assert_pairs_meet_in_two_ellipses,
    );
}

/// Each pair of patches is marched: its curves are sampled, the march
/// splitting the ellipses where they cross.
fn assert_pairs_marched(pairs: &[SurfaceSurfaceIntersection; 3]) {
    for (pair, result) in pairs.iter().enumerate() {
        let name = PAIR_NAMES[pair];
        assert_eq!(result.method(), IntersectionMethod::Numeric, "{name}");
        assert_eq!(result.completeness(), &Completeness::Sampled, "{name}");
        assert!(!result.curves().is_empty(), "{name}");
    }
}

#[test]
fn test_cylinders_as_bspline_patches_meet_at_the_eight_points() {
    // The march keeps to the surfaces within the confusion, and so do the
    // points found on its curves.
    let cylinders = unit_cylinders();
    let tol = Tolerance::DEFAULT;
    let limits = Limits {
        curve_from_cylinders: tol.confusion,
        point_from_expected: tol.confusion,
        point_from_cylinders: tol.confusion,
    };
    check_known_points(
        &patches(&cylinders),
        &cylinders,
        &unit_common_points(),
        tol,
        &limits,
        assert_pairs_marched,
    );
}

/// The unit cylinders with the Y one of radius `radius`, tilted by `tilt`
/// toward X and moved by `offset` along Z.
fn perturbed_cylinders(radius: f64, tilt: f64, offset: f64) -> [Cylinder; 3] {
    let [x, _, z] = unit_cylinders();
    let y = Cylinder::new(
        Point::new(0.0, 0.0, offset),
        Vector::new(tilt.sin(), tilt.cos(), 0.0),
        radius,
    )
    .unwrap();
    [x, y, z]
}

/// How many of `points` stand more than `merge` from every earlier one.
fn distinct_count(points: &[Point], merge: f64) -> usize {
    points
        .iter()
        .enumerate()
        .filter(|(i, point)| {
            points[..*i]
                .iter()
                .all(|other| point.distance(*other) > merge)
        })
        .count()
}

/// Where the common points have no closed form, the routes on `surfaces`
/// agree: each finds eight distinct points, every point has a partner
/// within the confusion in every other route, and lies within the
/// confusion of all three `cylinders`. The pairs' curves keep within the
/// confusion of their cylinders.
fn check_routes_agree<S>(surfaces: &[S; 3], cylinders: &[Cylinder; 3])
where
    S: Intersect<S, Output = SurfaceSurfaceIntersection>,
    TrimmedCurve: Intersect<S, Output = CurveSurfaceIntersection>
        + Intersect<TrimmedCurve, Output = CurveCurveIntersection>,
{
    let tol = Tolerance::DEFAULT;
    let pairs = pair_intersections(surfaces, tol);
    assert_curves_on_their_cylinders(&pairs, cylinders, tol.confusion);
    let routes = routes(surfaces, &pairs, tol);
    for route in &routes {
        assert_eq!(route.points.len(), 8, "{}: {:?}", route.name, route.points);
        assert_eq!(
            distinct_count(&route.points, tol.confusion),
            8,
            "{}: {:?}",
            route.name,
            route.points
        );
        for point in &route.points {
            let off = distance_to_all(*point, cylinders);
            assert!(
                off <= tol.confusion,
                "{}: {point:?} is {off:e} off a cylinder",
                route.name
            );
            for other in &routes {
                let gap = other
                    .points
                    .iter()
                    .map(|candidate| point.distance(*candidate))
                    .fold(f64::INFINITY, f64::min);
                assert!(
                    gap <= tol.confusion,
                    "{}: {point:?} is {gap:e} from the nearest point of {}",
                    route.name,
                    other.name
                );
            }
        }
    }
}

#[test]
fn test_routes_agree_with_the_y_cylinder_a_little_wider() {
    for widening in [1e-3, 1e-6, 1e-9, 1e-12] {
        let cylinders = perturbed_cylinders(1.0 + widening, 0.0, 0.0);
        check_routes_agree(&whole(&cylinders), &cylinders);
    }
}

#[test]
fn test_routes_agree_with_the_y_cylinder_a_little_tilted() {
    for tilt in [1e-3, 1e-6, 1e-9] {
        let cylinders = perturbed_cylinders(1.0, tilt, 0.0);
        check_routes_agree(&whole(&cylinders), &cylinders);
    }
}

#[test]
fn test_routes_agree_with_the_y_cylinder_a_little_off_centre() {
    for offset in [1e-3, 1e-6, 1e-9] {
        let cylinders = perturbed_cylinders(1.0, 0.0, offset);
        check_routes_agree(&whole(&cylinders), &cylinders);
    }
}

#[test]
fn test_routes_agree_on_perturbed_bspline_patches() {
    for cylinders in [
        perturbed_cylinders(1.0 + 1e-6, 0.0, 0.0),
        perturbed_cylinders(1.0, 1e-6, 0.0),
    ] {
        check_routes_agree(&patches(&cylinders), &cylinders);
    }
}
