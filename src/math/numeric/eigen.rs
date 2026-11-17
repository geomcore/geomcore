//! Eigenvalues and eigenvectors of a small symmetric matrix, by cyclic
//! Jacobi rotations.
//!
//! Each rotation zeroes one off-diagonal entry; a sweep over all of them
//! shrinks the off-diagonal part quadratically, so a 3 by 3 matrix settles
//! in a handful of sweeps. The eigenvectors come out orthonormal to
//! rounding, the smallest eigenvalue's included, which is what a fit by
//! orthogonal regression reads.

/// The most sweeps over the off-diagonal entries: Jacobi converges
/// quadratically, so a small matrix needs far fewer.
const MAX_SWEEPS: usize = 50;

/// The eigenvalues of the symmetric `matrix`, ascending, with a unit
/// eigenvector for each: `vectors[k]` goes with `values[k]`. Only the
/// upper triangle is read as the matrix's; the vectors are orthonormal.
#[allow(clippy::needless_range_loop)] // index arithmetic on a symmetric matrix
pub(crate) fn symmetric_eigen<const N: usize>(matrix: [[f64; N]; N]) -> ([f64; N], [[f64; N]; N]) {
    let mut a = matrix;
    for row in 0..N {
        for column in 0..row {
            a[row][column] = a[column][row];
        }
    }
    // The columns of `v` are the eigenvectors.
    let mut v = [[0.0; N]; N];
    for (k, row) in v.iter_mut().enumerate() {
        row[k] = 1.0;
    }
    for _ in 0..MAX_SWEEPS {
        let mut rotated = false;
        for p in 0..N {
            for q in p + 1..N {
                let off = a[p][q];
                let diagonal = a[p][p].abs() + a[q][q].abs();
                // An entry below the rounding of the diagonal it sits
                // between is zero.
                if diagonal + 100.0 * off.abs() == diagonal {
                    a[p][q] = 0.0;
                    a[q][p] = 0.0;
                    continue;
                }
                rotated = true;
                let theta = (a[q][q] - a[p][p]) / (2.0 * off);
                // The smaller root of t² + 2θt − 1 = 0, the rotation's
                // tangent; for a huge θ, its first-order form.
                let tangent = if theta.abs() > 1e150 {
                    0.5 / theta
                } else {
                    theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt())
                };
                let cosine = 1.0 / (tangent * tangent + 1.0).sqrt();
                let sine = tangent * cosine;
                for row in a.iter_mut() {
                    let (at_p, at_q) = (row[p], row[q]);
                    row[p] = cosine * at_p - sine * at_q;
                    row[q] = sine * at_p + cosine * at_q;
                }
                for column in 0..N {
                    let (at_p, at_q) = (a[p][column], a[q][column]);
                    a[p][column] = cosine * at_p - sine * at_q;
                    a[q][column] = sine * at_p + cosine * at_q;
                }
                for row in v.iter_mut() {
                    let (at_p, at_q) = (row[p], row[q]);
                    row[p] = cosine * at_p - sine * at_q;
                    row[q] = sine * at_p + cosine * at_q;
                }
            }
        }
        if !rotated {
            break;
        }
    }
    let mut order: [usize; N] = std::array::from_fn(|k| k);
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    let values = order.map(|k| a[k][k]);
    let vectors = order.map(|k| std::array::from_fn(|row| v[row][k]));
    (values, vectors)
}

#[cfg(test)]
mod tests;
