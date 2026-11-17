use super::*;

/// The largest entry of `a v - λ v` over every eigenpair.
fn largest_defect<const N: usize>(
    a: [[f64; N]; N],
    values: [f64; N],
    vectors: [[f64; N]; N],
) -> f64 {
    let mut largest = 0.0f64;
    for (value, vector) in values.iter().zip(vectors) {
        for row in 0..N {
            let product: f64 = (0..N).map(|k| a[row][k] * vector[k]).sum();
            largest = largest.max((product - value * vector[row]).abs());
        }
    }
    largest
}

/// The largest entry of `vᵀ v - I`.
fn largest_departure_from_orthonormal<const N: usize>(vectors: [[f64; N]; N]) -> f64 {
    let mut largest = 0.0f64;
    for i in 0..N {
        for j in 0..N {
            let dot: f64 = (0..N).map(|k| vectors[i][k] * vectors[j][k]).sum();
            let expected = if i == j { 1.0 } else { 0.0 };
            largest = largest.max((dot - expected).abs());
        }
    }
    largest
}

#[test]
fn a_diagonal_matrix_is_its_own_decomposition_sorted() {
    let (values, vectors) = symmetric_eigen([[3.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 2.0]]);
    assert_eq!(values, [-1.0, 2.0, 3.0]);
    assert_eq!(vectors, [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]);
}

#[test]
fn a_full_matrix_decomposes_to_rounding() {
    let a = [[4.0, 1.0, -2.0], [1.0, 2.0, 0.5], [-2.0, 0.5, 3.0]];
    let (values, vectors) = symmetric_eigen(a);
    assert!(values[0] <= values[1] && values[1] <= values[2]);
    assert!(largest_defect(a, values, vectors) < 1e-14);
    assert!(largest_departure_from_orthonormal(vectors) < 1e-15);
    let trace: f64 = values.iter().sum();
    assert!((trace - 9.0).abs() < 1e-14);
}

#[test]
fn a_repeated_eigenvalue_keeps_orthonormal_vectors() {
    // The covariance of points spread evenly on a circle: two equal
    // eigenvalues, and a zero one across the circle's plane.
    let a = [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.0]];
    let rotation = [[0.6, 0.8, 0.0], [-0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
    let tilt = [[1.0, 0.0, 0.0], [0.0, 0.8, -0.6], [0.0, 0.6, 0.8]];
    let product = |x: [[f64; 3]; 3], y: [[f64; 3]; 3]| -> [[f64; 3]; 3] {
        std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| x[i][k] * y[k][j]).sum()))
    };
    let transposed = |x: [[f64; 3]; 3]| -> [[f64; 3]; 3] {
        std::array::from_fn(|i| std::array::from_fn(|j| x[j][i]))
    };
    let r = product(tilt, rotation);
    let rotated = product(product(r, a), transposed(r));
    let (values, vectors) = symmetric_eigen(rotated);
    assert!(values[0].abs() < 1e-16);
    assert!((values[1] - 0.5).abs() < 1e-15 && (values[2] - 0.5).abs() < 1e-15);
    assert!(largest_defect(rotated, values, vectors) < 1e-15);
    assert!(largest_departure_from_orthonormal(vectors) < 1e-15);
    // The smallest eigenvector is the rotated z axis, up to its sign.
    let normal = [r[0][2], r[1][2], r[2][2]];
    let alignment: f64 = (0..3).map(|k| normal[k] * vectors[0][k]).sum();
    assert!((alignment.abs() - 1.0).abs() < 1e-15);
}

#[test]
fn a_badly_scaled_matrix_keeps_its_small_eigenvector() {
    // Eigenvalues 1, 1e-4 and 1e-10, as the covariance of a nearly flat
    // point set: the smallest eigenvalue comes out to the rounding of the
    // largest, and its eigenvector to that over the gap above it.
    let a = [1.0, 1e-4, 1e-10];
    let c = std::f64::consts::FRAC_1_SQRT_2;
    let r = [[c, -c, 0.0], [0.5, 0.5, -c], [0.5, 0.5, c]];
    let rotated: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..3).map(|k| r[i][k] * a[k] * r[j][k]).sum())
    });
    let (values, vectors) = symmetric_eigen(rotated);
    assert!((values[2] - 1.0).abs() < 1e-15);
    assert!((values[1] - 1e-4).abs() < 1e-15);
    assert!((values[0] - 1e-10).abs() < 1e-15);
    let smallest = [r[0][2], r[1][2], r[2][2]];
    let alignment: f64 = (0..3).map(|k| smallest[k] * vectors[0][k]).sum();
    assert!((alignment.abs() - 1.0).abs() < 1e-12);
}

#[test]
fn a_two_by_two_matrix_decomposes() {
    let a = [[2.0, 1.0], [1.0, 2.0]];
    let (values, vectors) = symmetric_eigen(a);
    assert!((values[0] - 1.0).abs() < 1e-15 && (values[1] - 3.0).abs() < 1e-15);
    assert!(largest_defect(a, values, vectors) < 1e-15);
}
