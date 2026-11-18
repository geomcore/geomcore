use super::*;

/// The lemma of the module: an arc of a unit circle and its chord, run at
/// constant speed over the same parameter, are never farther apart than
/// the sagitta, for arcs up to a half turn.
#[test]
fn an_arc_and_its_chord_at_the_same_parameter_are_within_the_sagitta() {
    for k in 1..=60 {
        let angle = std::f64::consts::PI * k as f64 / 60.0;
        let (start, end) = ((0.0f64, 1.0f64), (angle.sin(), angle.cos()));
        let bound = sagitta(angle);
        let farthest = (0..=200)
            .map(|i| {
                let t = i as f64 / 200.0;
                let (x, y) = ((angle * t).sin(), (angle * t).cos());
                let chord = (
                    start.0 + (end.0 - start.0) * t,
                    start.1 + (end.1 - start.1) * t,
                );
                (x - chord.0).hypot(y - chord.1)
            })
            .fold(0.0, f64::max);
        assert!(farthest <= bound * (1.0 + 1e-12), "angle {angle}");
        // Reached at the middle.
        assert!(farthest >= bound * (1.0 - 1e-12), "angle {angle}");
    }
}

#[test]
fn a_spindle_torus_crosses_its_axis_twice_a_turn() {
    let crossings = torus_axis_crossings(1.0, 2.0, (0.0, TAU));
    assert_eq!(crossings.len(), 2);
    for v in crossings {
        assert!((1.0 + 2.0 * v.cos()).abs() < 1e-15);
    }
    assert!(torus_axis_crossings(3.0, 1.0, (0.0, TAU)).is_empty());
    // Across the seam.
    assert_eq!(torus_axis_crossings(1.0, 2.0, (3.0, 3.0 + TAU)).len(), 2);
}
