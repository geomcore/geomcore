use super::*;

#[test]
fn test_elevation_matrix_rows_are_convex_weights() {
    for (p, by) in [(1, 1), (2, 3), (3, 1), (5, 4), (12, 13)] {
        let w = p + 1;
        let mut m = vec![0.0; (p + by + 1) * w];
        elevation_matrix(p, by, &mut m);
        for row in m.chunks_exact(w) {
            let sum: f64 = row.iter().sum();
            assert!((sum - 1.0).abs() < 1e-14, "p {p} by {by}: {row:?}");
            assert!(row.iter().all(|&x| x >= 0.0));
        }
    }
}
