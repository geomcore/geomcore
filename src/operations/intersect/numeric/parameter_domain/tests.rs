use super::*;

/// The period `[1, 3)`.
fn period() -> ParameterDomain {
    ParameterDomain::Periodic {
        start: 1.0,
        period: 2.0,
    }
}

#[test]
fn test_start_and_end_of_a_period() {
    assert_eq!(period().start(), 1.0);
    assert_eq!(period().end(), 3.0);
}

#[test]
fn test_start_and_end_of_an_interval() {
    let interval = ParameterDomain::Interval(-1.0, 4.0);
    assert_eq!(interval.start(), -1.0);
    assert_eq!(interval.end(), 4.0);
}

#[test]
fn test_from_bounds_is_a_period_or_an_interval() {
    let periodic = ParameterDomain::from_bounds((1.0, 3.0), true);
    assert!(periodic.is_periodic());
    assert_eq!((periodic.start(), periodic.end()), (1.0, 3.0));
    let interval = ParameterDomain::from_bounds((1.0, 3.0), false);
    assert!(!interval.is_periodic());
    assert_eq!((interval.start(), interval.end()), (1.0, 3.0));
}

#[test]
fn test_turn_is_a_whole_period() {
    assert!(ParameterDomain::TURN.is_periodic());
    assert_eq!(ParameterDomain::TURN.start(), 0.0);
    assert_eq!(ParameterDomain::TURN.end(), TAU);
}

#[test]
fn test_slack_is_a_fraction_of_the_width() {
    assert_eq!(period().slack(), SLACK * 2.0);
    assert_eq!(ParameterDomain::Interval(-1.0, 4.0).slack(), SLACK * 5.0);
}

#[test]
fn test_admit_wraps_a_period() {
    assert_eq!(period().admit(3.5), Some(1.5));
    assert_eq!(period().admit(-0.5), Some(1.5));
    assert_eq!(period().admit(2.0), Some(2.0));
}

#[test]
fn test_admit_reads_a_value_within_slack_of_a_period_end_as_its_start() {
    assert_eq!(period().admit(3.0 - 1e-12), Some(1.0));
    let beyond = period().admit(3.0 - 1e-6).unwrap();
    assert!((beyond - (3.0 - 1e-6)).abs() < 1e-12);
}

#[test]
fn test_admit_keeps_an_interval_value_inside() {
    assert_eq!(ParameterDomain::Interval(0.0, 1.0).admit(0.5), Some(0.5));
}

#[test]
fn test_admit_clamps_within_slack_of_an_interval() {
    let interval = ParameterDomain::Interval(0.0, 1.0);
    assert_eq!(interval.admit(1.0 + 5e-11), Some(1.0));
    assert_eq!(interval.admit(-5e-11), Some(0.0));
}

#[test]
fn test_admit_rejects_beyond_the_slack_of_an_interval() {
    let interval = ParameterDomain::Interval(0.0, 1.0);
    assert_eq!(interval.admit(1.0 + 1e-9), None);
    assert_eq!(interval.admit(-1e-9), None);
}

#[test]
fn test_settle_wraps_a_period() {
    assert_eq!(period().settle(3.5), 1.5);
}

#[test]
fn test_settle_clamps_an_interval_beyond_its_slack() {
    let interval = ParameterDomain::Interval(0.0, 1.0);
    assert_eq!(interval.settle(2.0), 1.0);
    assert_eq!(interval.settle(-5.0), 0.0);
    assert_eq!(interval.settle(0.25), 0.25);
}

#[test]
fn test_within_slack_of_widens_the_pair_by_the_slack() {
    let interval = ParameterDomain::Interval(0.0, 1.0);
    assert!(interval.within_slack_of(1.0 + 5e-11, (0.0, 1.0)));
    assert!(!interval.within_slack_of(1.0 + 1e-9, (0.0, 1.0)));
    assert!(interval.within_slack_of(0.2 - 1e-12, (0.2, 0.4)));
    assert!(!interval.within_slack_of(0.1, (0.2, 0.4)));
}
