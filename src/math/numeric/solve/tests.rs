use super::*;

#[test]
fn test_solve_2x2() {
    let x = solve_2x2([[2.0, 1.0], [1.0, 3.0]], [5.0, 6.0], 1e-12).unwrap();
    assert!((x[0] - 1.8).abs() < 1e-12);
    assert!((x[1] - 1.4).abs() < 1e-12);
    assert!(solve_2x2([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0], 1e-12).is_none());
}

#[test]
fn test_solve_3x3() {
    let x = solve_3x3(
        [[3.0, 2.0, -1.0], [2.0, -2.0, 4.0], [-1.0, 0.5, -1.0]],
        [1.0, -2.0, 0.0],
        1e-12,
    )
    .unwrap();
    assert!((x[0] - 1.0).abs() < 1e-9);
    assert!((x[1] + 2.0).abs() < 1e-9);
    assert!((x[2] + 2.0).abs() < 1e-9);
    assert!(solve_3x3([[1.0, 2.0, 3.0]; 3], [1.0, 2.0, 3.0], 1e-12).is_none());
}

#[test]
fn test_solve_4x4() {
    // Needs a row swap: the first pivot is zero.
    let a = [
        [0.0, 1.0, 2.0, 1.0],
        [1.0, 0.0, 1.0, 0.0],
        [2.0, 1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0, 3.0],
    ];
    let x = solve_4x4(a, [10.0, 4.0, 6.0, 16.0], 1e-12).unwrap();
    for (xi, want) in x.iter().zip([1.0, 0.0, 3.0, 4.0]) {
        assert!((xi - want).abs() < 1e-12, "{x:?}");
    }
    assert!(solve_4x4([[1.0, 2.0, 3.0, 4.0]; 4], [1.0; 4], 1e-12).is_none());
}

#[test]
fn test_newton_3d() {
    // Intersection of x^2 + y^2 = 2, y = x, z = 1 near (1, 1, 1).
    let x = newton_3d(
        |v| [v[0] * v[0] + v[1] * v[1] - 2.0, v[1] - v[0], v[2] - 1.0],
        |v| {
            [
                [2.0 * v[0], 2.0 * v[1], 0.0],
                [-1.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        },
        [1.2, 0.8, 0.5],
        1e-12,
        50,
    )
    .unwrap();
    assert!((x[0] - 1.0).abs() < 1e-9);
    assert!((x[1] - 1.0).abs() < 1e-9);
    assert!((x[2] - 1.0).abs() < 1e-9);
}

#[test]
fn test_polish_root_sqrt2() {
    let x = polish_root(|x| x * x - 2.0, |x| 2.0 * x, 1.4, 8);
    assert!((x - std::f64::consts::SQRT_2).abs() < 1e-15);
    // A flat function leaves the candidate unchanged.
    assert_eq!(polish_root(|_| 1.0, |_| 0.0, 0.5, 8), 0.5);
    // A step that would increase |f| is refused.
    assert_eq!(
        polish_root(|x: f64| x.atan(), |x| 1.0 / (1.0 + x * x), 3.0, 8),
        3.0
    );
}

#[test]
fn test_brent_minimum() {
    // min of (x - 2)^4 + (x - 2)^2 at x = 2.
    let x = brent_minimum(|x| (x - 2.0).powi(4) + (x - 2.0).powi(2), -5.0, 5.0, 1e-10);
    assert!((x - 2.0).abs() < 1e-7);
}
