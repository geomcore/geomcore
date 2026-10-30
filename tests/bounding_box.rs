//! Bounding boxes through the public interface: every box holds the
//! shape, sampled densely, and the box of an elementary shape is reached by
//! the shape on every face.

mod common;

use common::samples::trimmed;
use geomcore::prelude::*;
use geomcore::{
    BSplineCurve, BSplineSurface, BoundingBox, Circle, Cone, Curve, Cylinder, ElementaryCurve,
    ElementarySurface, Ellipse, Hyperbola, Line, Parabola, Plane, Point, Sphere, Surface, Torus,
    TrimmedCurve, Vector,
};
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI, TAU};

// ---------------------------------------------------------------- checker

fn coordinate(p: Point, k: usize) -> f64 {
    [p.x, p.y, p.z][k]
}

/// The largest value of `f` over `[a, b]`: dense samples, then a
/// golden-section refinement around the best one.
fn maximum(a: f64, b: f64, f: impl Fn(f64) -> f64) -> f64 {
    const SAMPLES: usize = 96;
    let step = (b - a) / SAMPLES as f64;
    let mut best = (f(a), a);
    for i in 1..=SAMPLES {
        let u = if i == SAMPLES { b } else { a + step * i as f64 };
        let value = f(u);
        if value > best.0 {
            best = (value, u);
        }
    }
    let (mut lo, mut hi) = ((best.1 - step).max(a), (best.1 + step).min(b));
    let ratio = 0.5 * (5.0_f64.sqrt() - 1.0);
    for _ in 0..60 {
        let (x1, x2) = (hi - ratio * (hi - lo), lo + ratio * (hi - lo));
        if f(x1) > f(x2) {
            hi = x2;
        } else {
            lo = x1;
        }
    }
    best.0.max(f(0.5 * (lo + hi)))
}

/// The extremes `(min, max)` of each coordinate of `point` over `[a, b]`.
fn extremes(a: f64, b: f64, point: impl Fn(f64) -> Point) -> [(f64, f64); 3] {
    std::array::from_fn(|k| {
        let low = -maximum(a, b, |t| -coordinate(point(t), k));
        let high = maximum(a, b, |t| coordinate(point(t), k));
        (low, high)
    })
}

/// The extremes of each coordinate of a surface over its parameter box.
fn surface_extremes<S: ParametricSurface>(
    surface: &S,
    (u0, u1): (f64, f64),
    (v0, v1): (f64, f64),
) -> [(f64, f64); 3] {
    std::array::from_fn(|k| {
        let value = |sign: f64| {
            maximum(u0, u1, |u| {
                maximum(v0, v1, |v| sign * coordinate(surface.eval_point(u, v), k))
            })
        };
        (-value(-1.0), value(1.0))
    })
}

fn scale(b: BoundingBox) -> f64 {
    1.0 + b.min().distance(Point::ORIGIN) + b.max().distance(Point::ORIGIN)
}

/// Asserts that the box holds every extreme and, if `tight`, that the
/// extremes reach every face.
fn check(b: BoundingBox, extremes: [(f64, f64); 3], tight: bool, what: &str) {
    let rounding = 1e-12 * scale(b);
    for (k, (low, high)) in extremes.into_iter().enumerate() {
        let (min, max) = (coordinate(b.min(), k), coordinate(b.max(), k));
        assert!(min <= low + rounding, "{what}: axis {k}, {min} > {low}");
        assert!(high <= max + rounding, "{what}: axis {k}, {high} > {max}");
        if tight {
            let reach = 1e-9 * scale(b);
            assert!(low - min <= reach, "{what}: axis {k}, min {min} vs {low}");
            assert!(max - high <= reach, "{what}: axis {k}, max {max} vs {high}");
        }
    }
}

/// The box of `shape` through the implementation for references.
fn through_reference<T: HasBoundingBox>(shape: &T) -> Option<BoundingBox> {
    <&T as HasBoundingBox>::bounding_box(&shape)
}

fn check_curve<C: ParametricCurve + HasBoundingBox>(curve: &C, tight: bool, what: &str) {
    let (first, last) = curve.bounds();
    let b = curve.bounding_box().expect("bounded curve");
    check(
        b,
        extremes(first, last, |t| curve.eval_point(t)),
        tight,
        what,
    );
}

/// Tilted directions: none along an axis, one along each.
fn normals() -> Vec<Vector> {
    vec![
        Vector::new(0.3, -0.5, 0.8),
        Vector::new(-1.0, 2.0, 0.25),
        Vector::new(0.0, 0.6, -0.8),
        Vector::X,
        Vector::Y,
        Vector::Z,
    ]
}

/// A direction orthogonal to `n`.
fn across(n: Vector) -> Vector {
    let other = if n.x.abs() < 0.9 {
        Vector::X
    } else {
        Vector::Y
    };
    n.cross(other)
}

/// Ranges of an angle: arcs on either side of the seam, across it, a
/// whole turn and a sliver.
const ARCS: [(f64, f64); 7] = [
    (0.0, FRAC_PI_2),
    (0.3, 2.9),
    (-1.0, 1.0),
    (5.0, 7.5),
    (1.0, 1.0 + TAU),
    (2.0, 2.001),
    (-4.0, -0.5),
];

// ---------------------------------------------------------------- elementary curves

#[test]
fn test_circles_and_arcs_are_tight() {
    let center = Point::new(1.0, -2.0, 3.0);
    for n in normals() {
        let circle = Circle::new(center, n, 2.5).unwrap();
        check_curve(&circle, true, "circle");
        for (first, last) in ARCS {
            let arc = TrimmedCurve::new(circle, first, last).unwrap();
            check_curve(&arc, true, &format!("arc {first}..{last} about {n:?}"));
        }
    }
}

#[test]
fn test_quarter_arc_box_is_its_extreme_points() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let arc = TrimmedCurve::new(circle, 0.25, FRAC_PI_2 + 0.25).unwrap();
    let b = arc.bounding_box().unwrap();
    // The arc passes through the top of the circle but not its sides.
    assert_eq!(b.max().y, 1.0);
    assert_eq!(b.min().x, arc.end_point().x);
    assert_eq!(b.max().x, arc.start_point().x);
    assert_eq!(b.min().y, arc.start_point().y);
}

#[test]
fn test_ellipses_and_arcs_are_tight() {
    let center = Point::new(-3.0, 0.5, 2.0);
    for n in normals() {
        let ellipse = Ellipse::new(center, n, across(n), 4.0, 1.5).unwrap();
        check_curve(&ellipse, true, "ellipse");
        for (first, last) in ARCS {
            let arc = TrimmedCurve::new(ellipse, first, last).unwrap();
            check_curve(&arc, true, &format!("elliptic arc {first}..{last}"));
        }
    }
}

#[test]
fn test_trimmed_parabolas_are_tight() {
    let apex = Point::new(0.5, 1.0, -1.0);
    for n in normals() {
        let parabola = Parabola::new(apex, n, across(n), 0.75).unwrap();
        for (first, last) in [(-3.0, 2.0), (0.5, 4.0), (-1e-3, 1e-3), (-10.0, -2.0)] {
            let piece = TrimmedCurve::new(parabola, first, last).unwrap();
            check_curve(&piece, true, &format!("parabola {first}..{last}"));
        }
    }
}

#[test]
fn test_trimmed_hyperbolas_are_tight() {
    let center = Point::new(2.0, 0.0, 1.0);
    for n in normals() {
        let hyperbola = Hyperbola::new(center, n, across(n), 2.0, 0.7).unwrap();
        for (first, last) in [(-1.5, 2.0), (0.2, 1.0), (-3.0, -0.1)] {
            let piece = TrimmedCurve::new(hyperbola, first, last).unwrap();
            check_curve(&piece, true, &format!("hyperbola {first}..{last}"));
        }
    }
}

#[test]
fn test_segments_are_their_end_points() {
    let (a, b) = (Point::new(1.0, 5.0, -2.0), Point::new(-3.0, 6.0, 0.0));
    let segment = TrimmedCurve::segment(a, b).unwrap();
    assert_eq!(segment.bounding_box(), Some(BoundingBox::new(a, b)));
    let line = Line::new(a, Vector::new(1.0, 1.0, 1.0)).unwrap();
    let piece = TrimmedCurve::new(line, -2.0, 3.0).unwrap();
    check_curve(&piece, true, "line piece");
}

#[test]
fn test_unbounded_curves_have_no_box() {
    let n = Vector::new(0.3, -0.5, 0.8);
    let curves: [ElementaryCurve; 3] = [
        Line::new(Point::ORIGIN, n).unwrap().into(),
        Parabola::new(Point::ORIGIN, n, across(n), 1.0)
            .unwrap()
            .into(),
        Hyperbola::new(Point::ORIGIN, n, across(n), 2.0, 1.0)
            .unwrap()
            .into(),
    ];
    for curve in curves {
        assert_eq!(curve.bounding_box(), None);
        assert_eq!(Curve::from(curve).bounding_box(), None);
    }
}

#[test]
fn test_run_time_curves_agree_with_concrete_ones() {
    let n = Vector::new(-1.0, 2.0, 0.25);
    let circle = Circle::new(Point::new(1.0, 2.0, 3.0), n, 2.0).unwrap();
    let ellipse = Ellipse::new(Point::ORIGIN, n, across(n), 3.0, 1.0).unwrap();
    let parabola = Parabola::new(Point::ORIGIN, n, across(n), 0.5).unwrap();
    let hyperbola = Hyperbola::new(Point::ORIGIN, n, across(n), 1.0, 2.0).unwrap();
    let line = Line::new(Point::ORIGIN, n).unwrap();
    let curves: [ElementaryCurve; 5] = [
        circle.into(),
        ellipse.into(),
        parabola.into(),
        hyperbola.into(),
        line.into(),
    ];
    for curve in curves {
        let wrapped = Curve::from(curve);
        assert_eq!(wrapped.bounding_box(), curve.bounding_box());
        let (first, last) = (0.5, 2.0);
        let piece = TrimmedCurve::new(curve, first, last).unwrap();
        let wrapped_piece = TrimmedCurve::new(wrapped, first, last).unwrap();
        assert!(piece.bounding_box().is_some());
        assert_eq!(piece.bounding_box(), wrapped_piece.bounding_box());
        assert_eq!(through_reference(&piece), piece.bounding_box());
    }
    assert_eq!(
        TrimmedCurve::new(circle, 0.5, 2.0).unwrap().bounding_box(),
        TrimmedCurve::new(Curve::from(circle), 0.5, 2.0)
            .unwrap()
            .bounding_box()
    );
    assert_eq!(through_reference(&circle), circle.bounding_box());
}

// ---------------------------------------------------------------- B-spline curves

fn cubic() -> BSplineCurve {
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, -1.0),
        Point::new(2.0, -1.0, 1.0),
        Point::new(3.0, 3.0, 0.5),
        Point::new(4.0, -2.0, 2.0),
        Point::new(5.0, 0.5, -0.5),
    ];
    BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.5, 3.0], vec![4, 1, 1, 4], false).unwrap()
}

/// A rational quadratic: the unit circle as four 90° arcs.
fn rational_circle() -> BSplineCurve {
    let s = FRAC_1_SQRT_2;
    let poles = vec![
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
        Point::new(-1.0, 1.0, 0.0),
        Point::new(-1.0, 0.0, 0.0),
        Point::new(-1.0, -1.0, 0.0),
        Point::new(0.0, -1.0, 0.0),
        Point::new(1.0, -1.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
    ];
    let weights = vec![1.0, s, 1.0, s, 1.0, s, 1.0, s, 1.0];
    BSplineCurve::new_rational(
        2,
        poles,
        weights,
        vec![0.0, 1.0, 2.0, 3.0, 4.0],
        vec![3, 2, 2, 2, 3],
        false,
    )
    .unwrap()
}

fn rational_cubic() -> BSplineCurve {
    let curve = cubic();
    let weights = vec![1.0, 4.0, 0.3, 2.0, 0.5, 1.0];
    BSplineCurve::new_rational(
        3,
        curve.poles().to_vec(),
        weights,
        curve.knots().to_vec(),
        curve.multiplicities().to_vec(),
        false,
    )
    .unwrap()
}

fn periodic() -> BSplineCurve {
    let poles = vec![
        Point::new(1.0, 0.0, 0.5),
        Point::new(0.0, 2.0, -0.5),
        Point::new(-1.5, 0.0, 1.0),
        Point::new(0.0, -1.0, 0.0),
        Point::new(2.0, -2.0, -1.0),
    ];
    BSplineCurve::new(
        3,
        poles,
        vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        vec![1; 6],
        true,
    )
    .unwrap()
}

fn bsplines() -> Vec<(&'static str, BSplineCurve)> {
    vec![
        ("cubic", cubic()),
        ("rational circle", rational_circle()),
        ("rational cubic", rational_cubic()),
        ("periodic", periodic()),
    ]
}

#[test]
fn test_bspline_boxes_hold_the_curve() {
    for (name, curve) in bsplines() {
        check_curve(&curve, false, name);
        let pole_box = BoundingBox::from_points(curve.poles().iter().copied()).unwrap();
        assert_eq!(curve.bounding_box(), Some(pole_box), "{name}");
        assert_eq!(Curve::from(curve.clone()).bounding_box(), Some(pole_box));
    }
}

#[test]
fn test_trimmed_bspline_boxes_hold_the_piece_and_reach_its_ends() {
    for (name, curve) in bsplines() {
        let (start, end) = curve.bounds();
        let length = end - start;
        let mut ranges = vec![
            (start, end),
            (start + 0.1 * length, start + 0.45 * length),
            (start + 0.37 * length, start + 0.38 * length),
            (start + 0.05 * length, start + 0.95 * length),
        ];
        if curve.is_periodic() {
            ranges.push((end - 0.3 * length, end + 0.2 * length));
            ranges.push((start - 2.0 * length, start - 1.1 * length));
        }
        let whole = curve.bounding_box().unwrap();
        for (first, last) in ranges {
            let what = format!("{name} over {first}..{last}");
            let piece = TrimmedCurve::new(&curve, first, last).unwrap();
            check_curve(&piece, false, &what);
            let b = piece.bounding_box().unwrap();
            let rounding = 1e-12 * scale(b);
            assert!(
                whole.enlarged(rounding).encloses(b.min())
                    && whole.enlarged(rounding).encloses(b.max()),
                "{what}: beyond the pole box"
            );
            for end_point in [piece.start_point(), piece.end_point()] {
                assert!(b.enlarged(rounding).encloses(end_point), "{what}");
            }
            let owned = TrimmedCurve::new(curve.clone(), first, last).unwrap();
            assert_eq!(owned.bounding_box(), Some(b));
            let wrapped = TrimmedCurve::new(Curve::from(curve.clone()), first, last).unwrap();
            assert_eq!(wrapped.bounding_box(), Some(b));
        }
    }
}

#[test]
fn test_trimmed_bspline_box_is_tighter_than_the_pole_box() {
    // The rational circle's poles reach the corners of the square around
    // it; a piece's Bézier control points stay near the piece.
    let curve = rational_circle();
    let piece = TrimmedCurve::new(&curve, 0.2, 0.8).unwrap();
    let b = piece.bounding_box().unwrap();
    assert!(b.max().y < 1.0 && b.min().x > 0.0, "{b:?}");
}

#[test]
fn test_trimmed_linear_bspline_box_is_exact() {
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(2.0, 2.0, 0.0),
        Point::new(4.0, 0.0, 1.0),
    ];
    let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
    let piece = TrimmedCurve::new(&curve, 0.5, 1.5).unwrap();
    check_curve(&piece, true, "polyline");
}

// ---------------------------------------------------------------- surfaces

#[test]
fn test_spheres_and_tori_are_tight() {
    let center = Point::new(1.0, -1.0, 2.0);
    let sphere = Sphere::new(center, 3.0).unwrap();
    let b = sphere.bounding_box().unwrap();
    check(
        b,
        surface_extremes(&sphere, sphere.u_bounds(), sphere.v_bounds()),
        true,
        "sphere",
    );
    for n in normals() {
        for (major, minor) in [(5.0, 1.5), (1.0, 2.0), (0.0, 1.0), (2.0, 0.0)] {
            let torus = Torus::new(center, n, major, minor).unwrap();
            let b = torus.bounding_box().unwrap();
            let extremes = surface_extremes(&torus, torus.u_bounds(), torus.v_bounds());
            check(
                b,
                extremes,
                true,
                &format!("torus {major}, {minor} about {n:?}"),
            );
        }
    }
}

#[test]
fn test_unbounded_surfaces_have_no_box() {
    let n = Vector::new(0.3, -0.5, 0.8);
    let surfaces: [ElementarySurface; 3] = [
        Plane::new(Point::ORIGIN, n).unwrap().into(),
        Cylinder::new(Point::ORIGIN, n, 1.0).unwrap().into(),
        Cone::new(Point::ORIGIN, n, 0.3, 1.0).unwrap().into(),
    ];
    for surface in surfaces {
        assert_eq!(surface.bounding_box(), None);
        assert_eq!(Surface::from(surface).bounding_box(), None);
    }
}

#[test]
fn test_run_time_surfaces_agree_with_concrete_ones() {
    let sphere = Sphere::new(Point::new(0.0, 1.0, 0.0), 2.0).unwrap();
    let torus = Torus::new(Point::ORIGIN, Vector::new(1.0, 1.0, 0.0), 3.0, 1.0).unwrap();
    for surface in [
        ElementarySurface::from(sphere),
        ElementarySurface::from(torus),
    ] {
        assert!(surface.bounding_box().is_some());
        assert_eq!(
            Surface::from(surface).bounding_box(),
            surface.bounding_box()
        );
        assert_eq!(through_reference(&surface), surface.bounding_box());
    }
}

fn bspline_surface(rational: bool) -> BSplineSurface {
    let poles: Vec<Vec<Point>> = (0..4)
        .map(|i| {
            (0..3)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    Point::new(x, y, (x * 1.3 + y * 0.7).sin() * 2.0)
                })
                .collect()
        })
        .collect();
    let (u_knots, u_mults) = (vec![0.0, 0.4, 1.0], vec![3, 1, 3]);
    let (v_knots, v_mults) = (vec![0.0, 1.0], vec![3, 3]);
    if rational {
        let weights = vec![
            vec![1.0, 2.0, 1.0],
            vec![0.5, 3.0, 1.0],
            vec![1.0, 0.2, 4.0],
            vec![1.0, 1.0, 1.0],
        ];
        BSplineSurface::new_rational(
            2, 2, poles, weights, u_knots, u_mults, v_knots, v_mults, false, false,
        )
        .unwrap()
    } else {
        BSplineSurface::new(
            2, 2, poles, u_knots, u_mults, v_knots, v_mults, false, false,
        )
        .unwrap()
    }
}

#[test]
fn test_bspline_surface_boxes_hold_the_surface() {
    for rational in [false, true] {
        let surface = bspline_surface(rational);
        let b = surface.bounding_box().unwrap();
        let (u, v) = (surface.u_bounds(), surface.v_bounds());
        check(
            b,
            surface_extremes(&surface, u, v),
            false,
            "B-spline surface",
        );
        let corners = [
            surface.eval_point(u.0, v.0),
            surface.eval_point(u.1, v.1),
            surface.eval_point(u.0, v.1),
        ];
        assert!(corners.iter().all(|&p| b.encloses(p)));
        assert_eq!(Surface::from(surface.clone()).bounding_box(), Some(b));
    }
}

// ---------------------------------------------------------------- every sample

fn is_bounded((first, last): (f64, f64)) -> bool {
    first.is_finite() && last.is_finite()
}

#[test]
fn test_every_sample_box_holds_the_shape() {
    for curve in common::samples::curves() {
        let what = format!("{curve:?}");
        let bounded = is_bounded(curve.bounds());
        assert_eq!(curve.bounding_box().is_some(), bounded, "{what}");
        if bounded {
            check_curve(&curve, false, &what);
        }
        check_curve(&trimmed(&curve), false, &format!("trimmed {what}"));
        if let Some(elementary) = curve.as_elementary() {
            assert_eq!(elementary.bounding_box().is_some(), bounded, "{what}");
            if bounded {
                check_curve(&elementary, false, &what);
            }
        }
    }
    for surface in common::samples::surfaces() {
        let what = format!("{surface:?}");
        let (u, v) = (surface.u_bounds(), surface.v_bounds());
        let bounded = is_bounded(u) && is_bounded(v);
        assert_eq!(surface.bounding_box().is_some(), bounded, "{what}");
        if let Some(b) = surface.bounding_box() {
            check(b, surface_extremes(&surface, u, v), false, &what);
        }
        if let Some(elementary) = surface.as_elementary() {
            assert_eq!(elementary.bounding_box(), surface.bounding_box(), "{what}");
        }
    }
}

// ---------------------------------------------------------------- culling

#[test]
fn test_boxes_cull_apart_shapes_and_keep_meeting_ones() {
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let circle = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 0.5).unwrap();
    let far = TrimmedCurve::new(
        Circle::new(Point::new(4.0, 0.0, 0.0), Vector::Z, 1.0).unwrap(),
        PI - 0.5,
        PI + 0.5,
    )
    .unwrap();
    let s = sphere.bounding_box().unwrap();
    assert!(s.overlaps(circle.bounding_box().unwrap()));
    let far = far.bounding_box().unwrap();
    assert!(!s.overlaps(far));
    // The arc's leftmost point is at x = 3, a distance of 2 from the sphere.
    assert!(!s.enlarged(1.9).overlaps(far));
    assert!(s.enlarged(1.0).overlaps(far.enlarged(1.0)));
}

// ---------------------------------------------------------------- 2D

mod planar {
    use geomcore::prelude::*;
    use geomcore::{
        BSplineCurve2D, BoundingBox2D, Circle2D, Curve2D, Ellipse2D, Hyperbola2D, Line2D,
        Parabola2D, Point2D, TrimmedCurve2D, Vector2D,
    };
    use std::f64::consts::{FRAC_PI_2, TAU};

    /// The extremes `(min, max)` of each coordinate of the curve `point`
    /// over `[a, b]`, from dense samples.
    fn sampled(a: f64, b: f64, point: impl Fn(f64) -> Point2D) -> [(f64, f64); 2] {
        let points: Vec<Point2D> = (0..=4000)
            .map(|i| point(a + (b - a) * i as f64 / 4000.0))
            .collect();
        [
            points
                .iter()
                .map(|p| p.x)
                .fold((f64::MAX, f64::MIN), minmax),
            points
                .iter()
                .map(|p| p.y)
                .fold((f64::MAX, f64::MIN), minmax),
        ]
    }

    fn minmax((lo, hi): (f64, f64), x: f64) -> (f64, f64) {
        (lo.min(x), hi.max(x))
    }

    /// The box holds every sample, and when `tight` the samples come within
    /// `slack` of every side.
    #[track_caller]
    fn check(b: BoundingBox2D, samples: [(f64, f64); 2], tight: bool) {
        let (min, max) = ([b.min().x, b.min().y], [b.max().x, b.max().y]);
        for k in 0..2 {
            assert!(min[k] <= samples[k].0 + 1e-12, "axis {k} low side");
            assert!(max[k] >= samples[k].1 - 1e-12, "axis {k} high side");
            if tight {
                assert!(min[k] >= samples[k].0 - 1e-5, "axis {k} low side is loose");
                assert!(max[k] <= samples[k].1 + 1e-5, "axis {k} high side is loose");
            }
        }
    }

    #[test]
    fn test_circles_and_arcs_are_tight() {
        let frames = [
            Circle2D::new(Point2D::new(1.0, -2.0), 1.5).unwrap(),
            Circle2D::from_frame(
                geomcore::Frame2D::new(
                    Point2D::ORIGIN,
                    Vector2D::new(3.0, 4.0),
                    -Vector2D::new(-4.0, 3.0),
                )
                .unwrap(),
                2.0,
            )
            .unwrap(),
        ];
        for circle in frames {
            check(
                circle.bounding_box().unwrap(),
                sampled(0.0, TAU, |t| circle.eval_point(t)),
                true,
            );
            for (first, last) in [
                (0.0, FRAC_PI_2),
                (-1.0, 2.0),
                (2.5, 2.5 + 5.0),
                (1.0, 1.0 + TAU),
            ] {
                let arc = TrimmedCurve2D::new(circle, first, last).unwrap();
                check(
                    arc.bounding_box().unwrap(),
                    sampled(first, last, |t| circle.eval_point(t)),
                    true,
                );
            }
        }
    }

    #[test]
    fn test_ellipses_and_arcs_are_tight() {
        let frames = [
            Ellipse2D::from_frame(geomcore::Frame2D::WORLD, 3.0, 1.5).unwrap(),
            Ellipse2D::from_frame(
                geomcore::Frame2D::new(
                    Point2D::new(1.0, -2.0),
                    Vector2D::new(0.6, 0.8),
                    Vector2D::new(-0.8, 0.6),
                )
                .unwrap(),
                2.5,
                1.25,
            )
            .unwrap(),
            Ellipse2D::from_frame(
                geomcore::Frame2D::new(
                    Point2D::new(-1.0, 0.5),
                    Vector2D::new(0.6, -0.8),
                    Vector2D::new(-0.8, -0.6),
                )
                .unwrap(),
                4.0,
                2.0,
            )
            .unwrap(),
        ];
        for ellipse in frames {
            check(
                ellipse.bounding_box().unwrap(),
                sampled(0.0, TAU, |t| ellipse.eval_point(t)),
                true,
            );
            for (first, last) in [
                (0.0, FRAC_PI_2),
                (-1.0, 2.0),
                (2.5, 2.5 + 5.0),
                (1.0, 1.0 + TAU),
            ] {
                let arc = TrimmedCurve2D::new(ellipse, first, last).unwrap();
                check(
                    arc.bounding_box().unwrap(),
                    sampled(first, last, |t| ellipse.eval_point(t)),
                    true,
                );
            }
            let wrapped = Curve2D::from(ellipse);
            assert_eq!(wrapped.bounding_box(), ellipse.bounding_box());
            let arc = TrimmedCurve2D::new(ellipse, 0.0, 2.0).unwrap();
            let boxed = TrimmedCurve2D::new(wrapped, 0.0, 2.0).unwrap();
            assert_eq!(boxed.bounding_box(), arc.bounding_box());
        }
    }

    #[test]
    fn test_parabolas_have_no_box_and_their_pieces_are_tight() {
        let parabolas = [
            Parabola2D::from_frame(geomcore::Frame2D::WORLD, 0.75).unwrap(),
            Parabola2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 1.5).unwrap(),
            Parabola2D::from_frame(
                geomcore::Frame2D::new(
                    Point2D::new(-1.0, 0.5),
                    Vector2D::new(0.6, -0.8),
                    Vector2D::new(-0.8, -0.6),
                )
                .unwrap(),
                0.5,
            )
            .unwrap(),
            Parabola2D::new(Point2D::ORIGIN, Vector2D::Y, 1.0).unwrap(),
        ];
        for parabola in parabolas {
            assert_eq!(parabola.bounding_box(), None);
            assert_eq!(Curve2D::from(parabola).bounding_box(), None);
            for (first, last) in [(-3.0, 2.0), (0.5, 4.0), (-1e-3, 1e-3), (-10.0, -2.0)] {
                let piece = TrimmedCurve2D::new(parabola, first, last).unwrap();
                check(
                    piece.bounding_box().unwrap(),
                    sampled(first, last, |t| parabola.eval_point(t)),
                    true,
                );
                let boxed = TrimmedCurve2D::new(Curve2D::from(parabola), first, last).unwrap();
                assert_eq!(boxed.bounding_box(), piece.bounding_box());
            }
        }
    }

    #[test]
    fn test_hyperbolas_have_no_box_and_their_pieces_are_tight() {
        let hyperbolas = [
            Hyperbola2D::from_frame(geomcore::Frame2D::WORLD, 2.0, 0.75).unwrap(),
            Hyperbola2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 1.5, 2.5).unwrap(),
            Hyperbola2D::from_frame(
                geomcore::Frame2D::new(
                    Point2D::new(-1.0, 0.5),
                    Vector2D::new(0.6, -0.8),
                    Vector2D::new(-0.8, -0.6),
                )
                .unwrap(),
                0.8,
                0.5,
            )
            .unwrap(),
            Hyperbola2D::new(Point2D::ORIGIN, Vector2D::Y, 1.0, 2.0).unwrap(),
        ];
        for hyperbola in hyperbolas {
            assert_eq!(hyperbola.bounding_box(), None);
            assert_eq!(Curve2D::from(hyperbola).bounding_box(), None);
            for (first, last) in [(-3.0, 2.0), (0.5, 4.0), (-1e-3, 1e-3), (-4.0, -0.2)] {
                let piece = TrimmedCurve2D::new(hyperbola, first, last).unwrap();
                check(
                    piece.bounding_box().unwrap(),
                    sampled(first, last, |t| hyperbola.eval_point(t)),
                    true,
                );
                let boxed = TrimmedCurve2D::new(Curve2D::from(hyperbola), first, last).unwrap();
                assert_eq!(boxed.bounding_box(), piece.bounding_box());
            }
        }
    }

    #[test]
    fn test_segments_are_their_end_points_and_lines_have_no_box() {
        let (start, end) = (Point2D::new(2.0, 1.0), Point2D::new(-1.0, 3.0));
        let segment = TrimmedCurve2D::segment(start, end).unwrap();
        assert_eq!(segment.bounding_box(), Some(BoundingBox2D::new(start, end)));
        assert_eq!(
            Line2D::new(Point2D::ORIGIN, Vector2D::Y)
                .unwrap()
                .bounding_box(),
            None
        );
    }

    #[test]
    fn test_bspline_boxes_hold_the_curve_and_pieces() {
        let poles = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 3.0),
            Point2D::new(2.0, -2.0),
            Point2D::new(4.0, 1.0),
        ];
        let curve =
            BSplineCurve2D::new(3, poles.clone(), vec![0.0, 1.0], vec![4, 4], false).unwrap();
        let whole = curve.bounding_box().unwrap();
        assert_eq!(whole, BoundingBox2D::from_points(poles).unwrap());
        check(whole, sampled(0.0, 1.0, |t| curve.eval_point(t)), false);
        let piece = TrimmedCurve2D::new(&curve, 0.2, 0.6).unwrap();
        let b = piece.bounding_box().unwrap();
        check(b, sampled(0.2, 0.6, |t| curve.eval_point(t)), false);
        assert!(whole.encloses(b.min()) && whole.encloses(b.max()));
    }

    #[test]
    fn test_run_time_curves_agree_with_concrete_ones() {
        let circle = Circle2D::new(Point2D::new(1.0, 1.0), 2.0).unwrap();
        let wrapped = Curve2D::from(circle);
        assert_eq!(wrapped.bounding_box(), circle.bounding_box());
        let arc = TrimmedCurve2D::new(circle, 0.0, 2.0).unwrap();
        let boxed = TrimmedCurve2D::new(wrapped, 0.0, 2.0).unwrap();
        assert_eq!(boxed.bounding_box(), arc.bounding_box());
    }

    #[test]
    fn test_every_sample_box_holds_the_curve() {
        for curve in crate::common::samples::curves_2d() {
            let (first, last) = curve.bounds();
            let bounded = first.is_finite() && last.is_finite();
            assert_eq!(curve.bounding_box().is_some(), bounded, "{curve:?}");
            if let Some(b) = curve.bounding_box() {
                check(b, sampled(first, last, |t| curve.eval_point(t)), false);
            }
            let piece = crate::common::samples::trimmed_2d(&curve);
            let (first, last) = piece.bounds();
            check(
                piece.bounding_box().unwrap(),
                sampled(first, last, |t| curve.eval_point(t)),
                false,
            );
        }
    }

    #[test]
    fn test_boxes_cull_apart_shapes() {
        let a = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
        let b = Circle2D::new(Point2D::new(5.0, 0.0), 1.0).unwrap();
        let (a, b) = (a.bounding_box().unwrap(), b.bounding_box().unwrap());
        assert!(!a.overlaps(b));
        let c = Circle2D::new(Point2D::new(1.5, 0.0), 1.0).unwrap();
        assert!(a.overlaps(c.bounding_box().unwrap()));
    }
}
