//! Curves with no closed-form deviations (a quartic branch): equal
//! parameter steps, each step's deviations measured on the curve.

use super::DiscretizationCriteria;
use super::elementary::equal_steps;
use crate::Point;
use crate::geometry::curves::Run;

/// One run of equal steps over `[first, last]` of the curve evaluated by
/// `point_at_parameter`, the fewest whose measured deviations meet the
/// criteria: each step's chord deviation is the distance of its middle
/// point from its chord, and its angle deviation twice the angle between
/// its chord and the chord to its middle, both exact for a circle and
/// close for any step short against the curve's turning.
pub(crate) fn sampled_steps<F>(
    point_at_parameter: impl Fn(f64) -> Point,
    (first, last): (f64, f64),
    criteria: &DiscretizationCriteria,
    min_segments: usize,
    visit: F,
) where
    F: FnMut(Run, &dyn Fn(f64) -> Point),
{
    let deviations = |count: usize| {
        let step = (last - first) / count as f64;
        let (mut chord, mut angle) = (0.0_f64, 0.0_f64);
        let mut start = point_at_parameter(first);
        for k in 0..count {
            let from = first + step * k as f64;
            let end = point_at_parameter(if k + 1 == count { last } else { from + step });
            let middle = point_at_parameter(from + 0.5 * step);
            let (across, half) = (end - start, middle - start);
            let length = across.square_magnitude();
            let off = if length > 0.0 {
                half - across * (half.dot(across) / length)
            } else {
                half
            };
            chord = chord.max(off.magnitude());
            angle = angle.max(2.0 * half.cross(across).magnitude().atan2(half.dot(across)));
            start = end;
        }
        (chord, angle)
    };
    // Deviations fall as the square of the count: from those of a few
    // steps, about the count that fits.
    const PROBE: usize = 8;
    let (chord, _) = deviations(PROBE);
    let guess = PROBE as f64 * (chord / criteria.max_chord_deviation).sqrt();
    equal_steps(
        &point_at_parameter,
        (first, last),
        criteria,
        min_segments,
        guess,
        deviations,
        visit,
    );
}
