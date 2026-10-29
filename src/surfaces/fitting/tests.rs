use super::*;

fn wavy_grid(rows: usize, cols: usize) -> Vec<Vec<Point>> {
    (0..rows)
        .map(|i| {
            (0..cols)
                .map(|j| {
                    let (x, y) = (i as f64 / (rows - 1) as f64, j as f64 / (cols - 1) as f64);
                    // Uneven spacing along both directions.
                    let (x, y) = (x + 0.1 * x * x, y - 0.05 * y * y);
                    Point::new(x, y, (3.0 * x).sin() * (2.0 * y).cos() + 0.2 * x * y)
                })
                .collect()
        })
        .collect()
}

/// Dense Gaussian elimination with partial pivoting.
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

/// The full basis row of `basis` at parameter `k`, `count` poles.
fn dense_row(basis: &BasisTable, k: usize, count: usize) -> Vec<f64> {
    let (first, values) = basis.at(k);
    let mut row = vec![0.0; count];
    row[first..first + values.len()].copy_from_slice(values);
    row
}

#[test]
fn test_approximation_matches_dense_tensor_least_squares() {
    for (rows, cols, p, q, mu, mv) in [(5, 6, 2, 1, 3, 4), (9, 7, 3, 3, 6, 5), (12, 12, 3, 2, 6, 6)]
    {
        let points = wavy_grid(rows, cols);
        let par = PointParametrization::Chordal;
        let surface = BSplineSurface::approximate(&points, p, q, mu, mv, par).unwrap();
        let grid = Grid::new(&points, par).unwrap();
        let u_basis = BasisTable::new(
            &approximation_knots(&grid.u_params, mu, p),
            p,
            &grid.u_params,
        );
        let v_basis = BasisTable::new(
            &approximation_knots(&grid.v_params, mv, q),
            q,
            &grid.v_params,
        );
        // One equation per grid point over the `mu * mv` poles.
        let mut normal = vec![vec![0.0; mu * mv]; mu * mv];
        let mut rhs = vec![Point::ORIGIN; mu * mv];
        for (i, line) in points.iter().enumerate() {
            let nu = dense_row(&u_basis, i, mu);
            for (j, &x) in line.iter().enumerate() {
                let nv = dense_row(&v_basis, j, mv);
                let row: Vec<f64> = nu
                    .iter()
                    .flat_map(|a| nv.iter().map(move |b| a * b))
                    .collect();
                for (r, &ar) in row.iter().enumerate() {
                    for (c, &ac) in row.iter().enumerate() {
                        normal[r][c] += ar * ac;
                    }
                    rhs[r] = Point::new(
                        rhs[r].x + ar * x.x,
                        rhs[r].y + ar * x.y,
                        rhs[r].z + ar * x.z,
                    );
                }
            }
        }
        let expected = dense_solve(normal, rhs);
        assert_eq!((surface.u_pole_count(), surface.v_pole_count()), (mu, mv));
        for (a, b) in surface.poles().iter().zip(&expected) {
            assert!(a.distance(*b) < 1e-10, "{a:?} vs {b:?}");
        }
    }
}

#[test]
fn test_approximation_within_keeps_grid_points_close() {
    let points = wavy_grid(25, 18);
    let par = PointParametrization::Centripetal;
    let grid = Grid::new(&points, par).unwrap();
    for (p, q) in [(1, 1), (2, 3), (3, 3)] {
        for tolerance in [1e-1, 1e-3, 1e-5, 0.0] {
            let surface =
                BSplineSurface::approximate_within(&points, p, q, tolerance, par).unwrap();
            for (i, row) in points.iter().enumerate() {
                for (j, x) in row.iter().enumerate() {
                    let at = surface.eval_point(grid.u_params[i], grid.v_params[j]);
                    assert!(
                        at.distance(*x) <= tolerance.max(1e-12),
                        "({p}, {q}) {tolerance}"
                    );
                }
            }
            if tolerance == 1e-1 {
                assert!(surface.u_pole_count() < 25 && surface.v_pole_count() < 18);
            }
            if tolerance == 0.0 {
                assert_eq!((surface.u_pole_count(), surface.v_pole_count()), (25, 18));
            }
        }
    }
}

#[test]
fn test_approximation_rejects_bad_grids() {
    let par = PointParametrization::Chordal;
    let points = wavy_grid(4, 5);
    let fixed =
        |points: &[Vec<Point>], p, mu| BSplineSurface::approximate(points, p, 1, mu, 3, par);
    assert_eq!(
        fixed(&points[..1], 1, 2),
        Err(BSplineConstructionError::TooFewPoints)
    );
    assert_eq!(
        fixed(&[], 1, 2),
        Err(BSplineConstructionError::TooFewPoints)
    );
    let mut ragged = points.clone();
    ragged[2].pop();
    assert_eq!(
        fixed(&ragged, 1, 2),
        Err(BSplineConstructionError::PoleCountMismatch)
    );
    assert_eq!(
        fixed(&points, 0, 2),
        Err(BSplineConstructionError::InvalidDegree)
    );
    assert_eq!(
        fixed(&points, 2, 2),
        Err(BSplineConstructionError::DegreeTooHigh)
    );
    assert_eq!(
        fixed(&points, 1, 5),
        Err(BSplineConstructionError::TooManyPoles)
    );
    let within = |p| BSplineSurface::approximate_within(&points, p, 1, 1e-3, par);
    assert_eq!(within(4), Err(BSplineConstructionError::DegreeTooHigh));
    assert!(within(3).is_ok());
    // Every row the same point: no `v` parameters.
    let collapsed = vec![vec![Point::ORIGIN; 3]; 3];
    assert_eq!(
        BSplineSurface::approximate(&collapsed, 1, 1, 2, 2, par),
        Err(BSplineConstructionError::ConfusedPoints)
    );
}

#[test]
fn test_parameters_skip_lines_without_length() {
    // A grid whose first row collapses to a point, as at a pole.
    let mut points = wavy_grid(4, 5);
    points[0] = vec![points[0][0]; 5];
    let grid = Grid::new(&points, PointParametrization::Chordal).unwrap();
    let rest = Grid::new(&points[1..], PointParametrization::Chordal).unwrap();
    for (a, b) in grid.v_params.iter().zip(&rest.v_params) {
        assert!((a - b).abs() < 1e-15);
    }
}
