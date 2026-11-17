use super::*;

#[test]
fn a_linear_model_is_solved_in_one_step() {
    // y = 2x + 1 through exact samples.
    let samples: Vec<(f64, f64)> = (0..5).map(|i| (i as f64, 2.0 * i as f64 + 1.0)).collect();
    let fitted = levenberg_marquardt(
        [0.0, 0.0],
        samples.len(),
        |&[slope, intercept]: &[f64; 2], k| {
            let (x, y) = samples[k];
            [(slope * x + intercept - y, [x, 1.0])]
        },
        |&[slope, intercept], [ds, di]| [slope + ds, intercept + di],
    );
    assert!((fitted[0] - 2.0).abs() < 1e-12 && (fitted[1] - 1.0).abs() < 1e-12);
}

#[test]
fn an_exponential_decay_is_recovered_from_a_far_start() {
    // y = 3 exp(-0.7 x): nonlinear in the rate.
    let samples: Vec<(f64, f64)> = (0..12)
        .map(|i| {
            let x = i as f64 * 0.5;
            (x, 3.0 * (-0.7 * x).exp())
        })
        .collect();
    let fitted = levenberg_marquardt(
        [1.0, 0.1],
        samples.len(),
        |&[amplitude, rate]: &[f64; 2], k| {
            let (x, y) = samples[k];
            let decay = (-rate * x).exp();
            [(amplitude * decay - y, [decay, -amplitude * x * decay])]
        },
        |&[amplitude, rate], [da, dr]| [amplitude + da, rate + dr],
    );
    assert!((fitted[0] - 3.0).abs() < 1e-12, "{fitted:?}");
    assert!((fitted[1] - 0.7).abs() < 1e-12, "{fitted:?}");
}

#[test]
fn two_rows_per_item_count_both() {
    // The point nearest three points in the plane: their mean, each point
    // giving two residual rows.
    let points = [(0.0, 0.0), (3.0, 0.0), (0.0, 6.0)];
    let fitted = levenberg_marquardt(
        [10.0, -4.0],
        points.len(),
        |&[x, y]: &[f64; 2], k| {
            let (px, py) = points[k];
            [(x - px, [1.0, 0.0]), (y - py, [0.0, 1.0])]
        },
        |&[x, y], [dx, dy]| [x + dx, y + dy],
    );
    assert!(
        (fitted[0] - 1.0).abs() < 1e-12 && (fitted[1] - 2.0).abs() < 1e-12,
        "{fitted:?}"
    );
}
