//! Samples of every variant of every shape enum, for the tests that check
//! a property on every shape. `src/shape_enum/tests.rs` includes this file
//! and fails, naming the variant, when a variant of `Curve`,
//! `ElementaryCurve`, `Curve2D`, `Surface` or `ElementarySurface` has no
//! sample here.
//!
//! Every sample sits in a tilted placement and has a point with a
//! coordinate of magnitude 2 or more, so its image under a scale of 1e308
//! about the origin overflows.

use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve, BSplineCurve2D, BSplineSurface, BezierCurve, BezierCurve2D, BezierSurface,
    Circle, Circle2D, Cone, Curve, Curve2D, Cylinder, Ellipse, Ellipse2D, ExtrusionSurface, Frame,
    Hyperbola, Hyperbola2D, Line, Line2D, Parabola, Parabola2D, Plane, Point, Point2D,
    RevolutionSurface, Sphere, Surface, Torus, TrimmedCurve, TrimmedCurve2D, Vector, Vector2D,
};

/// The placement of the samples: no axis along a coordinate axis.
pub fn tilted_frame() -> Frame {
    Frame::new(
        Point::new(1.0, -2.0, 0.5),
        Vector::new(0.3, -0.4, 0.86),
        Vector::new(1.0, 0.2, 0.0),
    )
    .unwrap()
}

/// The poles of the sample B-spline and Bézier curves.
fn curve_poles() -> Vec<Point> {
    vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.5),
        Point::new(3.0, 1.0, -0.5),
        Point::new(4.0, 3.0, 1.0),
    ]
}

/// A rational periodic quadratic B-spline curve: a closed loop around the
/// origin.
fn periodic_curve() -> BSplineCurve {
    let poles = (0..6)
        .map(|i| {
            let angle = i as f64 * std::f64::consts::TAU / 6.0;
            Point::new(2.0 * angle.cos(), angle.sin(), 0.2 * (2.0 * angle).sin())
        })
        .collect();
    let weights = vec![1.0, 0.7, 1.4, 1.0, 0.8, 1.2];
    let knots = (0..7).map(|i| i as f64).collect();
    BSplineCurve::new_rational(2, poles, weights, knots, vec![1; 7], true).unwrap()
}

/// Curves of every `Curve` variant: one of each, plus a rational periodic
/// B-spline and a rational Bézier curve.
pub fn curves() -> Vec<Curve> {
    let frame = tilted_frame();
    let (origin, z, x) = (frame.origin(), frame.z_direction(), frame.x_direction());
    vec![
        Line::new(origin, Vector::new(0.3, 0.9, -0.2))
            .unwrap()
            .into(),
        Circle::from_frame(frame, 2.0).unwrap().into(),
        Ellipse::new(origin, z, x, 3.0, 1.5).unwrap().into(),
        Parabola::new(origin, z, x, 0.8).unwrap().into(),
        Hyperbola::new(origin, z, x, 2.0, 1.0).unwrap().into(),
        BSplineCurve::new(2, curve_poles(), vec![-1.0, 0.5, 4.0], vec![3, 1, 3], false)
            .unwrap()
            .into(),
        periodic_curve().into(),
        BezierCurve::new(3, curve_poles()).unwrap().into(),
        BezierCurve::new_rational(3, curve_poles(), vec![1.0, 0.5, 2.0, 1.0])
            .unwrap()
            .into(),
    ]
}

/// The poles of the sample B-spline and Bézier curves in the plane.
fn curve_2d_poles() -> Vec<Point2D> {
    vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 1.0),
        Point2D::new(4.0, 3.0),
    ]
}

/// Curves of every `Curve2D` variant: one of each, plus a rational
/// periodic B-spline.
pub fn curves_2d() -> Vec<Curve2D> {
    let center = Point2D::new(0.5, -1.0);
    let x = Vector2D::new(0.6, 0.8);
    let periodic = periodic_curve();
    let periodic_poles = periodic
        .poles()
        .iter()
        .map(|pole| Point2D::new(pole.x, pole.y))
        .collect();
    vec![
        Line2D::new(center, x).unwrap().into(),
        Circle2D::new(center, 2.0).unwrap().into(),
        Ellipse2D::new(center, x, 3.0, 1.5).unwrap().into(),
        Parabola2D::new(center, x, 0.8).unwrap().into(),
        Hyperbola2D::new(center, x, 2.0, 1.0).unwrap().into(),
        BSplineCurve2D::new(
            2,
            curve_2d_poles(),
            vec![-1.0, 0.5, 4.0],
            vec![3, 1, 3],
            false,
        )
        .unwrap()
        .into(),
        BSplineCurve2D::new_rational(
            2,
            periodic_poles,
            periodic.weights().unwrap().to_vec(),
            periodic.knots().to_vec(),
            periodic.multiplicities().to_vec(),
            true,
        )
        .unwrap()
        .into(),
        BezierCurve2D::new(3, curve_2d_poles()).unwrap().into(),
    ]
}

/// The 3 by 3 poles of the sample B-spline and Bézier surfaces.
fn surface_poles() -> Vec<Vec<Point>> {
    (0..3)
        .map(|i| {
            (0..3)
                .map(|j| Point::new(i as f64, j as f64, ((i + 1) * (j + 2)) as f64 * 0.1))
                .collect()
        })
        .collect()
}

/// Weights of the rational sample surfaces.
fn surface_weights() -> Vec<Vec<f64>> {
    vec![
        vec![1.0, 2.0, 1.0],
        vec![0.5, 3.0, 1.0],
        vec![1.0, 0.2, 4.0],
    ]
}

/// Surfaces of every `Surface` variant: one of each elementary one;
/// extrusions of a line, a circle and a B-spline profile; revolutions of
/// a circle and of a line skew to the axis; and plain and rational
/// B-spline and Bézier surfaces.
pub fn surfaces() -> Vec<Surface> {
    let frame = tilted_frame();
    let circle = Circle::from_frame(
        Frame::new(Point::new(3.0, 0.5, 1.0), Vector::Y, Vector::X).unwrap(),
        0.7,
    )
    .unwrap();
    let line = Line::new(frame.origin(), Vector::new(0.3, 0.9, -0.2)).unwrap();
    let bspline_profile =
        BSplineCurve::new(2, curve_poles(), vec![-1.0, 0.5, 4.0], vec![3, 1, 3], false).unwrap();
    let sweep = Vector::new(0.2, -0.1, 1.5);
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let bspline = |weights: Option<Vec<Vec<f64>>>| {
        let (u_knots, v_knots) = (vec![-1.0, 5.0], vec![-2.0, 2.0]);
        match weights {
            None => BSplineSurface::new(
                2,
                2,
                surface_poles(),
                u_knots,
                vec![3, 3],
                v_knots,
                vec![3, 3],
                false,
                false,
            ),
            Some(weights) => BSplineSurface::new_rational(
                2,
                2,
                surface_poles(),
                weights,
                u_knots,
                vec![3, 3],
                v_knots,
                vec![3, 3],
                false,
                false,
            ),
        }
        .unwrap()
    };
    vec![
        Plane::from_frame(frame).into(),
        Cylinder::from_frame(frame, 2.0).unwrap().into(),
        Cone::from_frame(frame, 0.4, 2.0).unwrap().into(),
        Sphere::from_frame(frame, 2.0).unwrap().into(),
        Torus::from_frame(frame, 5.0, 1.5).unwrap().into(),
        ExtrusionSurface::new(line.into(), sweep).unwrap().into(),
        ExtrusionSurface::new(Circle::from_frame(frame, 2.0).unwrap().into(), sweep)
            .unwrap()
            .into(),
        ExtrusionSurface::new(bspline_profile.into(), sweep)
            .unwrap()
            .into(),
        RevolutionSurface::new(circle.into(), axis).unwrap().into(),
        RevolutionSurface::new(line.into(), axis).unwrap().into(),
        bspline(None).into(),
        bspline(Some(surface_weights())).into(),
        BezierSurface::new(2, 2, surface_poles()).unwrap().into(),
        BezierSurface::new_rational(2, 2, surface_poles(), surface_weights())
            .unwrap()
            .into(),
    ]
}

/// The range of the trimmed sample of a curve with these bounds: around
/// parameter 0 on an unbounded curve, inside the bounds of a bounded one.
pub fn trimmed_range((first, last): (f64, f64)) -> (f64, f64) {
    if first.is_finite() {
        (first + 0.2 * (last - first), first + 0.9 * (last - first))
    } else {
        (-1.5, 0.7)
    }
}

/// A trimmed piece of `curve`, over [`trimmed_range`].
pub fn trimmed(curve: &Curve) -> TrimmedCurve<Curve> {
    let (first, last) = trimmed_range(curve.bounds());
    TrimmedCurve::new(curve.clone(), first, last).unwrap()
}

/// A trimmed piece of a plane curve, over [`trimmed_range`].
pub fn trimmed_2d(curve: &Curve2D) -> TrimmedCurve2D<Curve2D> {
    let (first, last) = trimmed_range(curve.bounds());
    TrimmedCurve2D::new(curve.clone(), first, last).unwrap()
}

/// Parameters inside `bounds`, away from the ends; around 0 where the
/// bounds are infinite.
pub fn parameters((first, last): (f64, f64)) -> [f64; 3] {
    if first.is_finite() && last.is_finite() {
        [0.1, 0.45, 0.8].map(|fraction| first + fraction * (last - first))
    } else {
        [-1.3, 0.4, 2.1]
    }
}

/// Parameter pairs inside a surface's bounds, from [`parameters`].
pub fn surface_parameters(surface: &impl ParametricSurface) -> Vec<(f64, f64)> {
    let us = parameters(surface.u_bounds());
    let vs = parameters(surface.v_bounds());
    us.iter()
        .flat_map(|&u| vs.iter().map(move |&v| (u, v)))
        .collect()
}

/// Whether a curve's parameter is a length, and so scales with a
/// transform.
pub fn curve_parameter_is_length(curve: &Curve) -> bool {
    matches!(curve, Curve::Line(_) | Curve::Parabola(_))
}

/// Which of a surface's `(u, v)` parameters are lengths, and so scale
/// with a transform; the others are angles or B-spline parameters.
pub fn surface_parameters_are_lengths(surface: &Surface) -> (bool, bool) {
    match surface {
        Surface::Plane(_) => (true, true),
        Surface::Cylinder(_) | Surface::Cone(_) => (false, true),
        Surface::Extrusion(extrusion) => (curve_parameter_is_length(extrusion.profile()), true),
        Surface::Revolution(revolution) => (false, curve_parameter_is_length(revolution.profile())),
        _ => (false, false),
    }
}
