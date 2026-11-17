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
