use super::*;

#[test]
fn test_nearest_first_visits_by_box_distance_and_prunes() {
    // A 1 x 4 row of unit boxes along x.
    let leaves: Vec<Aabb> = (0..4)
        .map(|i| {
            let x = 2.0 * i as f64;
            [x, 0.0, 0.0, x + 1.0, 1.0, 1.0]
        })
        .collect();
    let tree = BoxTree::new(&leaves, 4);
    let mut seen = Vec::new();
    // The query sits right of the last box; each visited leaf reports
    // its own box distance as the best, so only the nearest is visited.
    tree.nearest_first(Point::new(10.0, 0.5, 0.5), |k| {
        seen.push(k);
        box_distance(Point::new(10.0, 0.5, 0.5), &leaves[k])
    });
    assert_eq!(seen, vec![3]);
    // Without pruning every leaf comes, nearest first.
    let mut all = Vec::new();
    tree.nearest_first(Point::new(10.0, 0.5, 0.5), |k| {
        all.push(k);
        f64::INFINITY
    });
    assert_eq!(all, vec![3, 2, 1, 0]);
}

#[test]
fn test_overlapping_visits_only_boxes_within_the_gap() {
    // A 2 x 3 grid of unit boxes, a step of 2 apart.
    let leaves: Vec<Aabb> = (0..6)
        .map(|k| {
            let (x, y) = (2.0 * (k / 3) as f64, 2.0 * (k % 3) as f64);
            [x, y, 0.0, x + 1.0, y + 1.0, 1.0]
        })
        .collect();
    let tree = BoxTree::new(&leaves, 3);
    let query = [1.2, 1.2, 0.5, 1.8, 2.5, 0.6];
    let mut seen = Vec::new();
    tree.overlapping(&query, 0.0, |k| seen.push(k));
    assert!(seen.is_empty());
    tree.overlapping(&query, 0.25, |k| seen.push(k));
    seen.sort();
    assert_eq!(seen, vec![0, 1, 3, 4]);
    let mut all = Vec::new();
    tree.overlapping(tree.bounds(), 0.0, |k| all.push(k));
    all.sort();
    assert_eq!(all, (0..6).collect::<Vec<_>>());
}

#[test]
fn test_box_distance() {
    let b = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    assert_eq!(box_distance(Point::new(0.5, 0.5, 0.5), &b), 0.0);
    assert_eq!(box_distance(Point::new(4.0, 5.0, 0.5), &b), 5.0);
}
