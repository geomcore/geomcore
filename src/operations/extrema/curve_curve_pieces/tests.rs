use std::f64::consts::PI;

use super::*;
use crate::operations::extrema::curve_curve::sampled_curve_curve_extrema;
use crate::{
    Circle, Curve, CurveCurveExtremum, Ellipse, Line, Parabola, Point, TrimmedCurve, Vector,
    curve_curve_extrema,
};

const TOL: Tolerance = Tolerance::DEFAULT;

/// A cubic swinging over the plane `z = 0`, never meeting it.
fn skew() -> BSplineCurve {
    let poles = [
        [-3.0, -1.0, 1.0],
        [-1.0, 2.0, 0.5],
        [1.0, -2.0, 1.5],
        [3.0, 1.0, 1.0],
    ];
    let poles = poles.map(|[x, y, z]| Point::new(x, y, z)).to_vec();
    BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap()
}

/// Asserts the minima are complete, are local minima of the distance
/// within the curves' bounds, and include those of the sampled search,
/// and returns them.
#[track_caller]
fn assert_extrema<A: ParametricCurve, B: ParametricCurve>(a: &A, b: &B) -> Vec<CurveCurveExtremum> {
    let extrema = curve_curve_extrema(a, b, TOL);
    assert_eq!(extrema.completeness(), &Completeness::Complete);
    let found = extrema.into_minima();
    assert!(!found.is_empty());
    let ((a0, a1), (b0, b1)) = (a.bounds(), b.bounds());
    for e in &found {
        let (s, t) = (e.first_parameter, e.second_parameter);
        assert!((a0..=a1).contains(&s) && (b0..=b1).contains(&t), "{e:?}");
        assert!((a.point_at(s).distance(b.point_at(t)) - e.distance).abs() < 1e-12);
        for (ds, dt) in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, -1.0)] {
            for h in [-1e-4, 1e-4] {
                let s = (s + h * ds).clamp(a0, a1);
                let t = (t + h * dt).clamp(b0, b1);
                let d = a.point_at(s).distance(b.point_at(t));
                assert!(d >= e.distance - 1e-12, "{e:?} not a minimum");
            }
        }
    }
    for want in sampled_curve_curve_extrema(a, b, TOL).into_minima() {
        let (pa, pb) = (
            a.point_at(want.first_parameter),
            b.point_at(want.second_parameter),
        );
        assert!(
            found
                .iter()
                .any(|e| (e.distance - want.distance).abs() < 1e-9
                    && a.point_at(e.first_parameter).distance(pa) < 1e-6
                    && b.point_at(e.second_parameter).distance(pb) < 1e-6),
            "{want:?} missing from {found:#?}"
        );
    }
    found
}

#[test]
fn test_extrema_of_segments_and_bspline_curve() {
    let curve = skew();
    let line = Line::new(Point::new(0.0, 0.0, -1.0), Vector::new(1.0, 0.3, 0.2)).unwrap();
    let far = Line::new(line.point_at(-5e3), line.direction()).unwrap();
    let segment = TrimmedCurve::new(line, -3.0, 3.0).unwrap();
    assert_extrema(&curve, &segment);
    let near = assert_extrema(&segment, &curve);
    let shifted = assert_extrema(&TrimmedCurve::new(far, 4997.0, 5003.0).unwrap(), &curve);
    assert_eq!(near.len(), shifted.len());
    for (a, b) in near.iter().zip(&shifted) {
        assert!((a.distance - b.distance).abs() < 1e-9);
        assert!((a.first_parameter + 5e3 - b.first_parameter).abs() < 1e-6);
    }
    let short = TrimmedCurve::new(line, -2.0, 1.0).unwrap();
    assert_extrema(&curve, &short);
    assert_sampled_extrema(&line, &curve);
    assert_sampled_extrema(&Curve::Line(line), &Curve::BSpline(curve.clone()));
}

/// Asserts the minima of `a` and `b` are the grid's, sampled.
#[track_caller]
fn assert_sampled_extrema<A: ParametricCurve, B: ParametricCurve>(a: &A, b: &B) {
    let extrema = curve_curve_extrema(a, b, TOL);
    assert_eq!(extrema.completeness(), &Completeness::Sampled);
    let grid = sampled_curve_curve_extrema(a, b, TOL);
    assert!(!grid.is_empty());
    assert_eq!(extrema.into_minima(), grid.into_minima());
}

#[test]
fn test_extrema_of_curved_elementary_curves_stay_sampled() {
    let curve = skew();
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.5).unwrap();
    let ellipse = Ellipse::new(
        Point::new(0.2, 0.0, 0.0),
        Vector::Z,
        Vector::new(1.0, 0.2, 0.0),
        2.5,
        1.0,
    )
    .unwrap();
    let parabola = Parabola::new(Point::new(0.0, -1.0, 0.0), Vector::Z, Vector::Y, 0.5).unwrap();
    assert_sampled_extrema(&circle, &curve);
    assert_sampled_extrema(&curve, &circle);
    assert_sampled_extrema(&ellipse, &curve);
    assert_sampled_extrema(&TrimmedCurve::new(circle, 1.0, 4.0).unwrap(), &curve);
    assert_sampled_extrema(&TrimmedCurve::new(ellipse, 0.0, PI).unwrap(), &curve);
    assert_sampled_extrema(&parabola, &curve);
    assert_sampled_extrema(&TrimmedCurve::new(parabola, -3.0, 3.0).unwrap(), &curve);
    let (a, b) = (Curve::Ellipse(ellipse), Curve::BSpline(curve.clone()));
    assert_sampled_extrema(&a, &b);
    assert_sampled_extrema(&Curve::Circle(circle), &curve);
}
