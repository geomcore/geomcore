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

#[test]
fn test_representative_shifts_by_whole_periods_into_the_range() {
    assert_eq!(representative(0.5, (0.0, 1.0), None), Some(0.5));
    assert_eq!(representative(1.0 + 1e-12, (0.0, 1.0), None), None);
    assert_eq!(representative(-1e-12, (0.0, 1.0), None), None);
    let shifted = representative(6.0, (-1.0, 1.0), Some(TAU)).unwrap();
    assert!((shifted - (6.0 - TAU)).abs() < 1e-15);
    assert_eq!(representative(3.0, (-1.0, 1.0), Some(TAU)), None);
}

#[test]
fn test_nearest_end_is_measured_across_the_seam() {
    assert_eq!(nearest_end(1.5, (0.0, 1.0), None), 1.0);
    assert_eq!(nearest_end(-0.5, (0.0, 1.0), None), 0.0);
    // Just past 2π against a range ending there: the end, not the start
    // that clamping would give.
    assert_eq!(nearest_end(5e-8, (PI, TAU), Some(TAU)), TAU);
    assert_eq!(nearest_end(TAU - 5e-8, (0.0, PI), Some(TAU)), 0.0);
}

/// The bare `start + (u − start).rem_euclid(period)` that
/// [`wrap_into_period`] replaces.
fn bare(u: f64, start: f64, period: f64) -> f64 {
    start + (u - start).rem_euclid(period)
}

#[test]
fn test_wrap_into_period_stays_below_the_end() {
    // (u, start, period): rows at and around the seam, a start of −π, and
    // multiples of the period far from the start.
    let rows = [
        (-0.0, 0.0, TAU),
        (-1e-17, 0.0, TAU),
        (-1e-300, 0.0, TAU),
        (TAU, 0.0, TAU),
        (TAU - 1e-16, 0.0, TAU),
        (0.0, 0.0, TAU),
        (1.0 - 1e-16, 1.0, TAU),
        (-PI - 1e-16, -PI, TAU),
        (PI, -PI, TAU),
        (1e6 * TAU, 0.0, TAU),
        (-1e6 * TAU - 1e-12, 0.0, TAU),
        (3.0 - 1e-16, 3.0, 2.0),
    ];
    for (u, start, period) in rows {
        let w = wrap_into_period(u, start, period);
        assert!(
            start <= w && w < start + period,
            "{u} from {start} by {period}: {w}"
        );
    }
}

#[test]
fn test_wrap_into_period_folds_the_end_onto_the_start() {
    // The bare expression returns the end of the period here.
    assert_eq!(bare(-1e-17, 0.0, TAU), TAU);
    assert_eq!(wrap_into_period(-1e-17, 0.0, TAU), 0.0);
    assert_eq!(wrap_angle(-1e-17), 0.0);
    assert_eq!(wrap_angle(TAU), 0.0);
}

#[test]
fn test_wrap_into_period_keeps_every_other_value() {
    for k in 0..2000 {
        let u = -50.0 + 0.0517 * f64::from(k);
        for (start, period) in [(0.0, TAU), (-PI, TAU), (1.5, 2.0)] {
            let b = bare(u, start, period);
            if b < start + period {
                assert_eq!(wrap_into_period(u, start, period).to_bits(), b.to_bits());
            }
        }
    }
    assert!(wrap_into_period(f64::NAN, 0.0, TAU).is_nan());
    assert!(wrap_angle(f64::INFINITY).is_nan());
}

#[test]
fn test_checked_refuses_in_order() {
    let open = (f64::NEG_INFINITY, f64::INFINITY);
    let check = |range, bounds, period| ParameterRange::checked(range, bounds, period);
    assert_eq!(
        check((f64::NAN, f64::INFINITY), open, None),
        Err(RangeError::NotANumber)
    );
    assert_eq!(
        check((0.0, f64::INFINITY), open, None),
        Err(RangeError::Infinite)
    );
    assert_eq!(check((1.0, 1.0), open, None), Err(RangeError::Empty));
    assert_eq!(check((2.0, 1.0), open, Some(TAU)), Err(RangeError::Empty));
    assert_eq!(
        check((0.0, 1.5), (0.0, 1.0), None),
        Err(RangeError::OutsideBounds)
    );
    assert_eq!(
        check((-0.5, 0.5), (0.0, 1.0), None),
        Err(RangeError::OutsideBounds)
    );
    assert_eq!(
        check((0.0, TAU + 1e-9), (0.0, TAU), Some(TAU)),
        Err(RangeError::ExceedsPeriod)
    );
}

#[test]
fn test_checked_accepts_a_period_from_any_start() {
    let bounds = (0.0, TAU);
    // A periodic range ignores the bounds and may cross the seam.
    for (first, last) in [
        (-1.0, 1.0),
        (5.0, 7.0),
        (-PI, PI),
        (0.0, TAU),
        (10.0, 10.0 + TAU),
    ] {
        let range = ParameterRange::checked((first, last), bounds, Some(TAU)).unwrap();
        assert_eq!(
            (range.first, range.last, range.period),
            (first, last, Some(TAU))
        );
    }
    // A whole turn computed with rounding errors is still within one period.
    assert!(ParameterRange::checked((0.1, 0.1 + TAU * (1.0 + 5e-13)), bounds, Some(TAU)).is_ok());
    assert!(ParameterRange::checked((0.2, 0.7), (0.0, 1.0), None).is_ok());
}

#[test]
fn test_spans_whole_period_at_the_boundary() {
    // (length, by the rounding band, within the default `parametric = 1e-9`)
    let rows = [
        (TAU, true, true),
        (TAU * (1.0 - 2e-12), false, true),
        (TAU - 5e-10, false, true),
        (TAU - 2e-9, false, false),
    ];
    for start in [0.0, -PI, 2.5] {
        for (length, rounding, within) in rows {
            let range = ParameterRange {
                first: start,
                last: start + length,
                period: Some(TAU),
            };
            assert_eq!(range.spans_whole_period(), rounding, "{start} {length}");
            assert_eq!(
                range.spans_whole_period_within(TOL),
                within,
                "{start} {length}"
            );
        }
    }
    let open = ParameterRange {
        first: 0.0,
        last: 1e9,
        period: None,
    };
    assert!(!open.spans_whole_period());
    assert!(!open.spans_whole_period_within(TOL));
}

#[test]
fn test_clip_range() {
    assert_eq!(
        *clip_range(
            (f64::NEG_INFINITY, f64::INFINITY),
            (2.0, 3.0),
            None,
            &unit_speed,
            TOL
        ),
        [(2.0, 3.0)]
    );
    assert_eq!(
        *clip_range((0.0, TAU), (-1.0, 1.0), Some(TAU), &unit_speed, TOL),
        [(-1.0, 1.0)]
    );
    assert!(clip_range((0.0, 1.0), (2.0, 3.0), None, &unit_speed, TOL).is_empty());
}

#[test]
fn test_clip_range_shifts_a_partial_overlap_across_the_seam() {
    let pieces = clip_range((6.0, 6.0 + 0.5), (-1.0, 1.0), Some(TAU), &unit_speed, TOL);
    assert_eq!(pieces.len(), 1);
    let (low, high) = pieces[0];
    assert!((low - (6.0 - TAU)).abs() < 1e-15);
    assert!((high - (6.5 - TAU)).abs() < 1e-15);
    assert!(clip_range((3.0, 4.0), (-1.0, 1.0), Some(TAU), &unit_speed, TOL).is_empty());
}

#[test]
fn test_clip_range_returns_both_pieces_across_the_seam() {
    // Arcs [0, 5] and [4, 7] of one circle share [4, 5] and [0, 7 - TAU].
    let pieces = clip_range((4.0, 7.0), (0.0, 5.0), Some(TAU), &unit_speed, TOL);
    assert_eq!(pieces.len(), 2);
    assert!((pieces[0].0 - 0.0).abs() < 1e-15 && (pieces[0].1 - (7.0 - TAU)).abs() < 1e-15);
    assert_eq!(pieces[1], (4.0, 5.0));
}

#[test]
fn test_clip_range_keeps_a_touch_as_a_single_point() {
    assert_eq!(
        *clip_range((1.0, 2.0), (0.0, 1.0), None, &unit_speed, TOL),
        [(1.0, 1.0)]
    );
    let seam = clip_range((0.0, 1.0), (5.0, TAU), Some(TAU), &unit_speed, TOL);
    assert_eq!(*seam, [(TAU, TAU)]);
    assert!(is_single_point(seam[0], &unit_speed, TOL));
}

#[test]
fn test_clip_range_drops_a_touch_repeating_the_end_of_a_piece() {
    // The copy shifted by a period touches the full range at TAU, the
    // same point as the start of the real piece.
    assert_eq!(
        *clip_range((0.0, 1.0), (0.0, TAU), Some(TAU), &unit_speed, TOL),
        [(0.0, 1.0)]
    );
}
