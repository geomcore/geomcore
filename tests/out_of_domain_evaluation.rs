//! Evaluation past a shape's bounds (ADR 0037): a non-periodic B-spline
//! or Bézier direction clamps to the nearer bound, point and derivatives
//! alike; elementary and trimmed curves continue along their equation;
//! operations that build from a parameter refuse one outside the bounds.

use geomcore::prelude::*;
use geomcore::{
    BSplineSplitError, BezierCurve, BezierSurface, KnotInsertionError, Line, Point,
    SurfaceIsoError, TrimmedCurve, Vector,
};

fn cubic() -> BezierCurve {
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.0),
        Point::new(3.0, 2.0, 1.0),
        Point::new(4.0, 0.0, 1.0),
    ];
    BezierCurve::new(3, poles).unwrap()
}

/// The grid `(i, j, i * j)`: the saddle `(2u, 2v, 4uv)`.
fn saddle() -> BezierSurface {
    let poles = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| Point::new(i as f64, j as f64, (i * j) as f64))
                .collect()
        })
        .collect();
    BezierSurface::new(2, 2, poles).unwrap()
}

#[test]
fn a_spline_curve_clamps_to_the_nearer_bound() {
    let bezier = cubic();
    let spline = bezier.to_bspline();
    let curves: [&dyn ParametricCurve; 2] = [&bezier, &spline];
    for curve in curves {
        for (outside, bound) in [(5.0, 1.0), (-3.0, 0.0)] {
            assert_eq!(curve.point_at(outside), curve.point_at(bound));
            for order in 1..=4 {
                assert_eq!(
                    curve.derivative(outside, order),
                    curve.derivative(bound, order)
                );
            }
        }
    }
}

#[test]
fn a_spline_surface_clamps_each_direction_to_its_nearer_bound() {
    let surface = saddle();
    assert_eq!(surface.point_at(99.0, 99.0), surface.point_at(1.0, 1.0));
    assert_eq!(surface.point_at(-2.0, 0.5), surface.point_at(0.0, 0.5));
    for (du, dv) in [(1, 0), (0, 1), (1, 1), (2, 0)] {
        assert_eq!(
            surface.derivative(7.0, -7.0, du, dv),
            surface.derivative(1.0, 0.0, du, dv)
        );
    }
}

#[test]
fn a_line_and_a_trimmed_line_continue_past_their_bounds() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let trimmed = TrimmedCurve::new(line, 0.0, 1.0).unwrap();
    assert_eq!(line.point_at(5.0), Point::new(5.0, 0.0, 0.0));
    assert_eq!(trimmed.point_at(5.0), Point::new(5.0, 0.0, 0.0));
}

#[test]
fn a_trimmed_spline_evaluates_its_basis_which_clamps_at_its_own_bounds() {
    let spline = cubic().to_bspline();
    let trimmed = TrimmedCurve::new(spline.clone(), 0.2, 0.6).unwrap();
    assert_eq!(trimmed.point_at(0.8), spline.point_at(0.8));
    assert_eq!(trimmed.point_at(5.0), spline.point_at(1.0));
}

#[test]
fn building_from_a_parameter_outside_the_bounds_is_refused() {
    let spline = cubic().to_bspline();
    assert_eq!(
        spline.split_at(5.0).err(),
        Some(BSplineSplitError::OutOfRange)
    );
    assert_eq!(
        spline.insert_knot(5.0, 1).err(),
        Some(KnotInsertionError::KnotOutOfRange)
    );
    let surface = saddle();
    assert_eq!(surface.u_iso(99.0).err(), Some(SurfaceIsoError::OutOfRange));
    assert_eq!(surface.v_iso(-1.0).err(), Some(SurfaceIsoError::OutOfRange));
}
