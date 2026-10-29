use super::*;
use crate::Point;
use crate::math::boxtree::boxes_apart;

/// The point at `f` in `[0, 1]` of the homogeneous Bézier `piece`.
fn bezier_point(piece: &[f64], f: f64) -> Point {
    let mut points: Vec<[f64; 4]> = piece.chunks(4).map(|c| [c[0], c[1], c[2], c[3]]).collect();
    for r in 1..points.len() {
        for i in 0..points.len() - r {
            points[i] =
                std::array::from_fn(|c| points[i][c] + f * (points[i + 1][c] - points[i][c]));
        }
    }
    let [x, y, z, w] = points[0];
    Point::new(x / w, y / w, z / w)
}

#[track_caller]
fn assert_pieces_follow(curve: &BSplineCurve, range: (f64, f64), breaks: &[f64]) {
    let pieces = CurvePieces::new(curve, range);
    assert_eq!(pieces.width, curve.degree() + 1);
    let ends: Vec<f64> = std::iter::once(range.0)
        .chain(pieces.ranges.iter().map(|r| r.1))
        .collect();
    assert_eq!(ends, breaks);
    for i in 0..pieces.len() {
        let (t0, t1) = pieces.ranges[i];
        for f in [0.0, 0.3, 0.5, 1.0] {
            let expected = curve.eval_point(t0 + f * (t1 - t0));
            assert!(bezier_point(pieces.piece(i), f).distance(expected) < 1e-12);
        }
    }
}

fn knotted(weights: Option<Vec<f64>>) -> BSplineCurve {
    let poles: Vec<Point> = (0..6)
        .map(|k| Point::new(k as f64, (k * k % 5) as f64, (k % 3) as f64))
        .collect();
    let (knots, mults) = (vec![0.0, 1.0, 1.5, 3.0], vec![4, 1, 1, 4]);
    match weights {
        None => BSplineCurve::new(3, poles, knots, mults, false),
        Some(w) => BSplineCurve::new_rational(3, poles, w, knots, mults, false),
    }
    .unwrap()
}

#[test]
fn test_curve_pieces_follow_the_curve() {
    let weights = vec![1.0, 0.5, 2.0, 1.5, 0.7, 1.0];
    for curve in [knotted(None), knotted(Some(weights))] {
        assert_pieces_follow(&curve, (0.0, 3.0), &[0.0, 1.0, 1.5, 3.0]);
        assert_pieces_follow(&curve, (0.2, 1.2), &[0.2, 1.0, 1.2]);
        assert_pieces_follow(&curve, (1.1, 1.4), &[1.1, 1.4]);
        assert_pieces_follow(&curve, (1.0, 1.5), &[1.0, 1.5]);
    }
}

#[test]
fn test_curve_pieces_cross_a_periodic_seam() {
    let poles = (0..6)
        .map(|k| {
            let angle = k as f64 * std::f64::consts::TAU / 6.0;
            Point::new(2.0 * angle.cos(), 2.0 * angle.sin(), (k % 2) as f64)
        })
        .collect();
    let knots = (0..=6).map(f64::from).collect();
    let ring = BSplineCurve::new(3, poles, knots, vec![1; 7], true).unwrap();
    assert_pieces_follow(&ring, (0.0, 6.0), &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    assert_pieces_follow(&ring, (4.5, 7.25), &[4.5, 5.0, 6.0, 7.0, 7.25]);
    assert_pieces_follow(&ring, (-1.5, 0.5), &[-1.5, -1.0, 0.0, 0.5]);
}

#[test]
fn test_pool_halves_into_two_slots() {
    let mut pool = Pool::new(4 * 3);
    let slot = pool.take();
    // A quadratic from (0, 0, 0) to (4, 0, 0) through control (2, 2, 0).
    pool.get_mut(slot)
        .copy_from_slice(&[0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 0.0, 1.0, 4.0, 0.0, 0.0, 1.0]);
    let upper = pool.halve(slot, (3, 1, 4), true);
    assert_ne!(upper, slot);
    assert_eq!(
        pool.get(slot),
        &[0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0]
    );
    assert_eq!(
        pool.get(upper),
        &[2.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 4.0, 0.0, 0.0, 1.0]
    );
    pool.release(slot);
    assert_eq!(pool.take(), slot);
}

const SQUARE: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
];

#[test]
fn test_separating_axis_along_the_patch_normal() {
    // The plane `z = x + y`, and points 0.5 and 0.6 above it.
    let tilted = SQUARE.map(|[x, y, _]| [x, y, x + y]);
    let normal = net_normal(&tilted, 2, 2);
    let above = [[0.9, 0.1, 1.5], [0.1, 0.9, 1.6]];
    let apart = |gap: f64| boxes_apart(&bounding_box(&above), &bounding_box(&tilted), gap);
    assert!(!apart(0.0));
    assert!(apart_along(normal, &above, &tilted, 0.25));
    assert!(!apart_along(normal, &above, &tilted, 0.3));
}

fn flatten(points: &[[f64; 3]]) -> Vec<f64> {
    points.iter().flatten().copied().collect()
}

fn cones(points: &[[f64; 3]]) -> Option<[Cone; 2]> {
    derivative_cones(&flatten(points), points, (2, 2, 3))
}

#[test]
fn test_normals_never_parallel() {
    let flat = cones(&SQUARE).unwrap();
    let tilted = cones(&SQUARE.map(|[x, y, _]| [x, y, x])).unwrap();
    let lifted = cones(&SQUARE.map(|[x, y, _]| [x, y, 1.0])).unwrap();
    // Crossing planes never have parallel normals; parallel ones do.
    assert!(normals_never_parallel(flat, tilted));
    assert!(normals_never_parallel(tilted, flat));
    assert!(!normals_never_parallel(flat, lifted));
    // A slightly warped square: never parallel to the tilted plane, but
    // parallel at one point to a shallow one.
    let warped = cones(&[
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.1],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, -0.1],
    ])
    .unwrap();
    let shallow = cones(&SQUARE.map(|[x, y, _]| [x, y, -0.05 * x])).unwrap();
    assert!(normals_never_parallel(warped, tilted));
    assert!(!normals_never_parallel(warped, shallow));
    // The rational square's cones come from interval bounds.
    let rational: Vec<f64> = SQUARE
        .iter()
        .zip([1.0, 1.2, 0.9, 1.1])
        .flat_map(|([x, y, z], w)| [x * w, y * w, z * w, w])
        .collect();
    let weighted = derivative_cones(&rational, &SQUARE, (2, 2, 4)).unwrap();
    assert!(normals_never_parallel(weighted, tilted));
    assert!(!normals_never_parallel(weighted, lifted));
}

#[test]
fn test_derivative_cones_reject_coincident_points() {
    let collapsed = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
    ];
    assert!(cones(&collapsed).is_none());
    let line = SQUARE.map(|[x, _, _]| [x, 0.0, 0.0]);
    assert!(cones(&line).is_none());
}

#[test]
fn test_stretch_and_thickness() {
    let long_in_s = SQUARE.map(|[x, y, z]| [4.0 * x, y, z]);
    assert!(stretches_along_rows(&long_in_s, 2, 2));
    let long_in_t = SQUARE.map(|[x, y, z]| [x, 4.0 * y, z]);
    assert!(!stretches_along_rows(&long_in_t, 2, 2));
    let tilted = SQUARE.map(|[x, y, _]| [x, y, x + y]);
    assert!(thickness(&tilted, net_normal(&tilted, 2, 2)) < 1e-15);
    assert_eq!(thickness(&tilted, [0.0, 0.0, 3.0]), 2.0);
    assert_eq!(thickness(&tilted, [0.0; 3]), f64::INFINITY);
}

#[test]
fn test_meets_at_most_once() {
    let polynomial: Vec<f64> = SQUARE.iter().flatten().copied().collect();
    // The same square, rational: weighted corners, still flat.
    let rational: Vec<f64> = SQUARE
        .iter()
        .zip([1.0, 1.2, 0.9, 1.1])
        .flat_map(|([x, y, z], w)| [x * w, y * w, z * w, w])
        .collect();
    for (net, h) in [(polynomial, 3), (rational, 4)] {
        let once = |curve: &[[f64; 3]]| meets_at_most_once(curve, &net, &SQUARE, (2, 2, h));
        // Steep through the square: once, even curving a little.
        assert!(once(&[[0.5, 0.5, -1.0], [0.6, 0.4, 0.0], [0.5, 0.5, 1.0]]));
        // Bending back, or lying in its plane: maybe more.
        assert!(!once(&[
            [0.5, 0.5, -1.0],
            [0.5, 0.5, 1.0],
            [0.6, 0.5, -1.0]
        ]));
        assert!(!once(&[[0.0, 0.5, 0.0], [1.0, 0.5, 0.0]]));
    }
    // A line meets a plane once, however shallow; a polynomial net's
    // differences say so, the rational one's interval bounds need not.
    let net: Vec<f64> = SQUARE.iter().flatten().copied().collect();
    assert!(meets_at_most_once(
        &[[0.0, 0.5, -0.01], [1.0, 0.5, 0.01]],
        &net,
        &SQUARE,
        (2, 2, 3)
    ));
}
