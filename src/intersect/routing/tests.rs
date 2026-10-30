use super::*;
use std::f64::consts::{PI, TAU};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The speed of a unit-speed line or of a circle of radius 1.
fn unit_speed(_: f64) -> f64 {
    1.0
}

#[test]
fn test_in_range_plain() {
    assert_eq!(in_range(0.5, (0.0, 1.0), None, unit_speed, TOL), Some(0.5));
    assert_eq!(in_range(1.5, (0.0, 1.0), None, unit_speed, TOL), None);
    assert_eq!(
        in_range(1.0 + 1e-12, (0.0, 1.0), None, unit_speed, TOL),
        Some(1.0)
    );
    assert_eq!(
        in_range(-1e-12, (0.0, 1.0), None, unit_speed, TOL),
        Some(0.0)
    );
}

#[test]
fn test_in_range_measures_the_excess_in_space() {
    let fast = |_: f64| 10.0;
    assert_eq!(in_range(1.0 + 5e-9, (0.0, 1.0), None, fast, TOL), Some(1.0));
    assert_eq!(in_range(1.0 + 2e-8, (0.0, 1.0), None, fast, TOL), None);
    assert_eq!(
        in_range(1.0 + 5e-8, (0.0, 1.0), None, unit_speed, TOL),
        Some(1.0)
    );
}

#[test]
fn test_in_range_wraps_across_the_seam() {
    let shifted = in_range(6.0, (-1.0, 1.0), Some(TAU), unit_speed, TOL).unwrap();
    assert!((shifted - (6.0 - TAU)).abs() < 1e-15);
    assert_eq!(in_range(3.0, (-1.0, 1.0), Some(TAU), unit_speed, TOL), None);
    // Just past 2π against a range ending there: the end, not the start.
    assert_eq!(
        in_range(5e-8, (PI, TAU), Some(TAU), unit_speed, TOL),
        Some(TAU)
    );
    assert_eq!(
        in_range(TAU - 5e-8, (0.0, PI), Some(TAU), unit_speed, TOL),
        Some(0.0)
    );
}

#[test]
fn test_in_range_accepts_just_below_the_start_on_a_periodic_curve() {
    let found = in_range(-1.0 - 1e-12, (-1.0, 1.0), Some(TAU), unit_speed, TOL).unwrap();
    assert_eq!(found, -1.0);
    assert_eq!(
        in_range(-1.0 - 1e-3, (-1.0, 1.0), Some(TAU), unit_speed, TOL),
        None
    );
}

/// The speed of a unit-speed line or of a circle of radius 1.
const UNIT_SPEED: fn(f64) -> f64 = |_| 1.0;

#[test]
fn test_clip_range() {
    assert_eq!(
        *clip_range(
            (f64::NEG_INFINITY, f64::INFINITY),
            (2.0, 3.0),
            None,
            &UNIT_SPEED,
            TOL
        ),
        [(2.0, 3.0)]
    );
    assert_eq!(
        *clip_range((0.0, TAU), (-1.0, 1.0), Some(TAU), &UNIT_SPEED, TOL),
        [(-1.0, 1.0)]
    );
    assert!(clip_range((0.0, 1.0), (2.0, 3.0), None, &UNIT_SPEED, TOL).is_empty());
}

#[test]
fn test_clip_range_shifts_a_partial_overlap_across_the_seam() {
    let pieces = clip_range((6.0, 6.0 + 0.5), (-1.0, 1.0), Some(TAU), &UNIT_SPEED, TOL);
    assert_eq!(pieces.len(), 1);
    let (low, high) = pieces[0];
    assert!((low - (6.0 - TAU)).abs() < 1e-15);
    assert!((high - (6.5 - TAU)).abs() < 1e-15);
    assert!(clip_range((3.0, 4.0), (-1.0, 1.0), Some(TAU), &UNIT_SPEED, TOL).is_empty());
}

#[test]
fn test_clip_range_returns_both_pieces_across_the_seam() {
    // Arcs [0, 5] and [4, 7] of one circle share [4, 5] and [0, 7 - TAU].
    let pieces = clip_range((4.0, 7.0), (0.0, 5.0), Some(TAU), &UNIT_SPEED, TOL);
    assert_eq!(pieces.len(), 2);
    assert!((pieces[0].0 - 0.0).abs() < 1e-15 && (pieces[0].1 - (7.0 - TAU)).abs() < 1e-15);
    assert_eq!(pieces[1], (4.0, 5.0));
}

#[test]
fn test_clip_range_keeps_a_touch_as_a_single_point() {
    assert_eq!(
        *clip_range((1.0, 2.0), (0.0, 1.0), None, &UNIT_SPEED, TOL),
        [(1.0, 1.0)]
    );
    let seam = clip_range((0.0, 1.0), (5.0, TAU), Some(TAU), &UNIT_SPEED, TOL);
    assert_eq!(*seam, [(TAU, TAU)]);
    assert!(is_single_point(seam[0], &UNIT_SPEED, TOL));
}

#[test]
fn test_clip_range_drops_a_touch_repeating_the_end_of_a_piece() {
    // The copy shifted by a period touches the full range at TAU, the
    // same point as the start of the real piece.
    assert_eq!(
        *clip_range((0.0, 1.0), (0.0, TAU), Some(TAU), &UNIT_SPEED, TOL),
        [(0.0, 1.0)]
    );
}

#[test]
fn test_trimmed_operand_reports_basis_and_range() {
    let circle = crate::Circle::new(crate::Point::ORIGIN, crate::Vector::Z, 1.0).unwrap();
    let arc = crate::curves::TrimmedCurve::new(circle, -1.0, 1.0).unwrap();
    let (basis, range) = arc.basis_and_range();
    // A trimmed circle routes through the typed circle solver.
    assert_eq!(basis, circle);
    assert_eq!(range, Some((-1.0, 1.0)));
}

#[test]
fn test_trimmed_of_trimmed_keeps_the_outer_range() {
    let line = crate::Line::new(crate::Point::ORIGIN, crate::Vector::X).unwrap();
    let outer = TrimmedCurve::new(TrimmedCurve::new(line, 0.0, 10.0).unwrap(), 2.0, 3.0).unwrap();
    assert_eq!(outer.basis_and_range().1, Some((2.0, 3.0)));
    let spline_free: Curve = line.into();
    assert_eq!(spline_free.basis_and_range().1, None);
}
