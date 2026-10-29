use super::*;

#[test]
fn test_snap_distance() {
    assert_eq!(snap_distance(1e-9, 1e-7), 0.0);
    assert_eq!(snap_distance(1e-7, 1e-7), 0.0);
    assert_eq!(snap_distance(1e-6, 1e-7), 1e-6);
}
