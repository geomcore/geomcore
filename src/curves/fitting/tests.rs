use super::*;

/// Dense Gaussian elimination with partial pivoting, the reference.
fn dense_solve(mut a: Vec<Vec<f64>>, mut b: Vec<Point>) -> Vec<Point> {
    let n = b.len();
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .unwrap();
        a.swap(k, pivot);
        b.swap(k, pivot);
        for i in k + 1..n {
            let f = a[i][k] / a[k][k];
            let pivot_row = a[k].clone();
            for (x, y) in a[i][k..].iter_mut().zip(&pivot_row[k..]) {
                *x -= f * y;
            }
            let bk = b[k];
            b[i] = Point::new(b[i].x - f * bk.x, b[i].y - f * bk.y, b[i].z - f * bk.z);
        }
    }
    for i in (0..n).rev() {
        let mut x = b[i];
        for j in i + 1..n {
            x = Point::new(
                x.x - a[i][j] * b[j].x,
                x.y - a[i][j] * b[j].y,
                x.z - a[i][j] * b[j].z,
            );
        }
        b[i] = Point::new(x.x / a[i][i], x.y / a[i][i], x.z / a[i][i]);
    }
    b
}

/// The full basis matrix row at `t`, by the Cox-de Boor recurrence.
fn dense_row(flat: &[f64], p: usize, t: f64) -> Vec<f64> {
    let n = flat.len() - p - 1;
    let mut rows = BasisRows::new(flat, p);
    let first = rows.at(t);
    let mut row = vec![0.0; n];
    row[first..=first + p].copy_from_slice(&rows.values[..=p]);
    row
}

fn wavy_points(count: usize) -> Vec<Point> {
    (0..count)
        .map(|i| {
            let s = i as f64 / (count - 1) as f64;
            Point::new(
                3.0 * s,
                (7.0 * s).sin() + 0.3 * s * s,
                (3.0 * s).cos() * 0.5,
            )
        })
        .collect()
}

fn assert_points_close(a: &[Point], b: &[Point], tol: f64) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        assert!(x.distance(*y) <= tol, "{x:?} vs {y:?}");
    }
}

#[test]
fn test_basis_rows_sum_to_one_and_match_full_rows() {
    let points = wavy_points(17);
    let params = point_parameters(&points, PointParametrization::Chordal).unwrap();
    for p in 1..=5 {
        let flat = averaged_knots(&params, p);
        let mut rows = BasisRows::new(&flat, p);
        for &t in &params {
            rows.at(t);
            let sum: f64 = rows.values[..=p].iter().sum();
            assert!((sum - 1.0).abs() < 1e-14, "degree {p} at {t}: {sum}");
        }
    }
}

#[test]
fn test_averaged_knots_match_definition() {
    let params = [0.0, 0.1, 0.25, 0.3, 0.6, 0.8, 1.0];
    let p = 3;
    let flat = averaged_knots(&params, p);
    assert_eq!(flat.len(), params.len() + p + 1);
    for j in 1..params.len() - p {
        let mean = params[j..j + p].iter().sum::<f64>() / p as f64;
        assert!((flat[j + p] - mean).abs() < 1e-15);
    }
    assert_eq!(&flat[..=p], &[0.0; 4]);
    assert_eq!(&flat[flat.len() - p - 1..], &[1.0; 4]);
}

#[test]
fn test_interpolation_matches_dense_solve() {
    for (count, p) in [(2, 1), (5, 2), (9, 3), (40, 3), (23, 5)] {
        let points = wavy_points(count);
        let curve = BSplineCurve::interpolate(&points, p, PointParametrization::Chordal).unwrap();
        let params = point_parameters(&points, PointParametrization::Chordal).unwrap();
        let flat = averaged_knots(&params, p);
        let a = params.iter().map(|&t| dense_row(&flat, p, t)).collect();
        let expected = dense_solve(a, points.clone());
        assert_points_close(curve.poles(), &expected, 1e-11);
        for (q, &t) in points.iter().zip(&params) {
            assert!(curve.eval_point(t).distance(*q) < 1e-11);
        }
    }
}

#[test]
fn test_approximation_matches_dense_normal_equations() {
    for (count, p, m) in [(4, 1, 2), (12, 3, 7), (60, 3, 20), (31, 4, 9), (8, 2, 8)] {
        let points = wavy_points(count);
        let curve =
            BSplineCurve::approximate(&points, p, m, PointParametrization::Centripetal).unwrap();
        let params = point_parameters(&points, PointParametrization::Centripetal).unwrap();
        let flat = approximation_knots(&params, m, p);
        let rows: Vec<Vec<f64>> = params.iter().map(|&t| dense_row(&flat, p, t)).collect();
        let mut normal = vec![vec![0.0; m]; m];
        let mut rhs = vec![Point::ORIGIN; m];
        for (row, q) in rows.iter().zip(&points) {
            for i in 0..m {
                for j in 0..m {
                    normal[i][j] += row[i] * row[j];
                }
                rhs[i] = Point::new(
                    rhs[i].x + row[i] * q.x,
                    rhs[i].y + row[i] * q.y,
                    rhs[i].z + row[i] * q.z,
                );
            }
        }
        assert_points_close(curve.poles(), &dense_solve(normal, rhs), 1e-10);
    }
}

#[test]
fn test_approximation_rejects_poles_without_data() {
    // Two clusters a millionth wide and far apart leave a 12-pole fit
    // numerically undetermined: the basis functions barely reach the
    // points.
    let mut points: Vec<Point> = (0..6)
        .map(|i| Point::new(i as f64 * 1e-6, 0.0, 0.0))
        .collect();
    points.extend((0..6).map(|i| Point::new(10.0 + i as f64 * 1e-6, 1.0, 0.0)));
    assert_eq!(
        BSplineCurve::approximate(&points, 3, 12, PointParametrization::Chordal),
        Err(BSplineConstructionError::PoleCountMismatch)
    );
}

#[test]
fn test_approximation_within_keeps_points_close_with_few_poles() {
    let points = wavy_points(200);
    for p in 1..=4 {
        for tolerance in [1e-1, 1e-3, 1e-5] {
            let curve = BSplineCurve::approximate_within(
                &points,
                p,
                tolerance,
                PointParametrization::Chordal,
            )
            .unwrap();
            let params = point_parameters(&points, PointParametrization::Chordal).unwrap();
            for (q, &t) in points.iter().zip(&params) {
                assert!(curve.eval_point(t).distance(*q) <= tolerance, "degree {p}");
            }
            assert!(curve.poles().len() < points.len(), "degree {p}");
        }
    }
}

#[test]
fn test_approximation_within_gathers_knots_where_points_bend() {
    // Straight for x in [0, 2], then a tight wave over [2, 3], joined
    // smoothly enough for a cubic.
    let points: Vec<Point> = (0..300)
        .map(|i| {
            let x = 3.0 * i as f64 / 299.0;
            let t = (x - 2.0).max(0.0);
            let y = t.powi(4) * (12.0 * t).sin();
            Point::new(x, y, 0.0)
        })
        .collect();
    let curve =
        BSplineCurve::approximate_within(&points, 3, 1e-4, PointParametrization::Chordal).unwrap();
    // The wave, under a third of the length, holds most interior knots.
    let start = curve
        .project_point(Point::new(2.0, 0.0, 0.0), crate::Tolerance::DEFAULT)
        .parameter;
    let interior = &curve.knots()[1..curve.knots().len() - 1];
    let in_wave = interior.iter().filter(|&&u| u > start).count();
    assert!(
        4 * in_wave > 3 * interior.len(),
        "{in_wave} of {}",
        interior.len()
    );
}

#[test]
fn test_approximation_within_zero_deviation_interpolates() {
    let points = wavy_points(12);
    for par in [
        PointParametrization::Chordal,
        PointParametrization::Centripetal,
    ] {
        let within = BSplineCurve::approximate_within(&points, 3, 0.0, par).unwrap();
        assert_eq!(within, BSplineCurve::interpolate(&points, 3, par).unwrap());
        let nan = BSplineCurve::approximate_within(&points, 3, f64::NAN, par).unwrap();
        assert_eq!(nan, within);
    }
}

#[test]
fn test_approximation_within_rejects_bad_input() {
    let points = wavy_points(4);
    let fit = |points: &[Point], degree| {
        BSplineCurve::approximate_within(points, degree, 1e-3, PointParametrization::Chordal)
    };
    assert_eq!(
        fit(&points[..1], 1),
        Err(BSplineConstructionError::TooFewPoints)
    );
    assert_eq!(
        fit(&points, 0),
        Err(BSplineConstructionError::InvalidDegree)
    );
    assert_eq!(
        fit(&points, 4),
        Err(BSplineConstructionError::DegreeTooHigh)
    );
    let confused = [points[0], points[0], points[1]];
    assert_eq!(
        fit(&confused, 1),
        Err(BSplineConstructionError::ConfusedPoints)
    );
}

#[test]
fn test_fewest_poles_search_finds_the_threshold() {
    for threshold in 3..=40 {
        let mut attempts = 0;
        let found = fit_with_fewest_poles(3, 41, |count| {
            attempts += 1;
            (count >= threshold).then_some(count)
        });
        assert_eq!(found, (threshold < 41).then_some(threshold));
        assert!(attempts <= 12, "{threshold}: {attempts}");
    }
}

#[test]
fn test_interpolation_of_two_points_is_their_segment() {
    let points = [Point::new(1.0, 2.0, 3.0), Point::new(4.0, 6.0, 3.0)];
    let curve = BSplineCurve::interpolate(&points, 1, PointParametrization::Chordal).unwrap();
    assert_eq!(curve.poles(), &points);
    assert_eq!(curve.knots(), &[0.0, 1.0]);
    assert_eq!(curve.multiplicities(), &[2, 2]);
}
