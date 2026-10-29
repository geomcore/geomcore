use super::*;

/// `lo + 1.0 * (hi - lo)` can round one step past `hi`: a projection
/// onto the upper edge must still land inside the parameter range.
#[test]
fn test_projection_on_the_upper_edge_stays_in_the_parameter_range() {
    let poles = vec![
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
        vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0)],
    ];
    let surface = BSplineSurface::new(
        1,
        1,
        poles,
        vec![0.03, 0.3],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap();
    let projection = surface.project_point(Point::new(2.0, 0.5, 0.0), Tolerance::DEFAULT);
    assert_eq!(projection.distance, 1.0, "{projection:?}");
    assert!(projection.u <= 0.3, "{projection:?}");
    assert!((projection.v - 0.5).abs() < 1e-12, "{projection:?}");
}

/// A flat wall whose first u span is 1e-9 wide. Points on it, including
/// inside that span, project to themselves; the narrow u direction must
/// not stall the much wider v direction or stop before the step moves a
/// visible distance in space.
#[test]
fn test_projection_inside_a_narrow_u_span() {
    let xs = [
        0.0, 40.0, 80.0, 120.0, 149.5, 149.75, 150.0, 150.25, 150.5, 200.0, 240.0, 280.0,
    ];
    let poles = xs
        .iter()
        .map(|&x| vec![Point::new(x, 0.0, -1.0), Point::new(x, 0.0, 1.0)])
        .collect::<Vec<_>>();
    let knots = vec![0.0, 1e-9, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    let surface = BSplineSurface::new(
        3,
        1,
        poles,
        knots,
        vec![4, 1, 1, 1, 1, 1, 1, 1, 1, 4],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap();
    let tol = Tolerance::DEFAULT;
    for u in [0.0, 1e-10, 5e-10, 1e-9, 0.5] {
        let point = surface.eval_point(u, 0.5);
        let projection = surface.project_point(point, tol);
        assert!(
            projection.distance <= tol.confusion,
            "u = {u:e}: {projection:?}"
        );
        assert!(
            projection.point.distance(point) <= tol.confusion,
            "u = {u:e}: {projection:?}"
        );
        assert!(
            (projection.v - 0.5).abs() <= 1e-9,
            "u = {u:e}: {projection:?}"
        );
    }
}

// ---- bilinear patch by hand ----

fn bilinear() -> BSplineSurface {
    // Degree 1x1 over the unit square:
    //   P00=(0,0,0) P01=(0,1,0)
    //   P10=(1,0,0) P11=(1,1,1)
    // S(u,v) = (u, v, u*v).
    let poles = vec![
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
        vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    ];
    BSplineSurface::new(
        1,
        1,
        poles,
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_accessors() {
    let s = bilinear();
    assert_eq!(s.u_degree(), 1);
    assert_eq!(s.v_degree(), 1);
    assert!(!s.is_u_periodic());
    assert!(!s.is_v_periodic());
    assert!(!s.is_rational());
}

#[test]
fn test_bilinear_point_by_hand() {
    let s = bilinear();
    assert_eq!(s.eval_point(0.5, 0.5), Point::new(0.5, 0.5, 0.25));
    assert_eq!(s.eval_point(0.0, 0.0), Point::new(0.0, 0.0, 0.0));
    assert_eq!(s.eval_point(1.0, 1.0), Point::new(1.0, 1.0, 1.0));
    assert_eq!(s.eval_point(0.3, 0.7), Point::new(0.3, 0.7, 0.3 * 0.7));
}

#[test]
fn test_bilinear_derivatives_by_hand() {
    let s = bilinear();
    // Su = (1, 0, v), Sv = (0, 1, u).
    assert_eq!(
        s.eval_derivative(0.3, 0.7, 1, 0),
        Vector::new(1.0, 0.0, 0.7)
    );
    assert_eq!(
        s.eval_derivative(0.3, 0.7, 0, 1),
        Vector::new(0.0, 1.0, 0.3)
    );
}

#[test]
fn test_eval_points_matches_mapped_eval_point() {
    let s = bilinear();
    let uvs = [(0.0, 0.0), (0.2, 0.4), (0.5, 0.5), (1.0, 1.0)];
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| s.eval_point(u, v)).collect();
    assert_eq!(s.eval_points(&uvs), expected);
}

#[test]
#[should_panic(expected = "du + dv must be >= 1")]
fn test_eval_derivative_zero_order_panics() {
    bilinear().eval_derivative(0.5, 0.5, 0, 0);
}

#[test]
fn test_eval_derivative_second_order_bilinear() {
    // S = (u, v, u*v): Suu = Svv = 0, Suv = (0, 0, 1).
    let s = bilinear();
    assert_eq!(
        s.eval_derivative(0.5, 0.5, 2, 0),
        Vector::new(0.0, 0.0, 0.0)
    );
    assert_eq!(
        s.eval_derivative(0.5, 0.5, 1, 1),
        Vector::new(0.0, 0.0, 1.0)
    );
    assert_eq!(
        s.eval_derivative(0.5, 0.5, 0, 2),
        Vector::new(0.0, 0.0, 0.0)
    );
}

#[test]
#[should_panic(expected = "du + dv must be <= 2")]
fn test_eval_derivative_third_order_panics() {
    bilinear().eval_derivative(0.5, 0.5, 2, 1);
}

// ---- bounds / periods ----

#[test]
fn test_bounds_clamped() {
    let s = bilinear();
    assert_eq!(s.u_bounds(), (0.0, 1.0));
    assert_eq!(s.v_bounds(), (0.0, 1.0));
    assert_eq!(s.u_period(), None);
    assert_eq!(s.v_period(), None);
}

// ---- validation error variants ----

fn ok_grid() -> Vec<Vec<Point>> {
    vec![
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
        vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    ]
}

#[test]
fn test_empty_grid_is_pole_count_mismatch() {
    let err = BSplineSurface::new(
        1,
        1,
        vec![],
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap_err();
    assert_eq!(err, BSplineConstructionError::PoleCountMismatch);
}

#[test]
fn test_ragged_grid_is_pole_count_mismatch() {
    let poles = vec![
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
        vec![Point::new(1.0, 0.0, 0.0)], // short row
    ];
    let err = BSplineSurface::new(
        1,
        1,
        poles,
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap_err();
    assert_eq!(err, BSplineConstructionError::PoleCountMismatch);
}

#[test]
fn test_invalid_u_degree() {
    let err = BSplineSurface::new(
        0,
        1,
        ok_grid(),
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap_err();
    assert_eq!(err, BSplineConstructionError::InvalidDegree);
}

#[test]
fn test_u_pole_count_mismatch_from_knots() {
    // Degree 1: expected u rows = sum(u_mults) - degree - 1. Knots
    // [0,1,2] mults [2,1,2] (valid for degree 1) expect 3 rows, but the
    // grid has only 2 => pole-count mismatch in the u direction.
    let err = BSplineSurface::new(
        1,
        1,
        ok_grid(),
        vec![0.0, 1.0, 2.0],
        vec![2, 1, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap_err();
    assert_eq!(err, BSplineConstructionError::PoleCountMismatch);
}

#[test]
fn test_v_knots_not_increasing() {
    let err = BSplineSurface::new(
        1,
        1,
        ok_grid(),
        vec![0.0, 1.0],
        vec![2, 2],
        vec![1.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap_err();
    assert_eq!(err, BSplineConstructionError::KnotsNotIncreasing);
}

#[test]
fn test_weight_count_mismatch_ragged() {
    let err = BSplineSurface::new_rational(
        1,
        1,
        ok_grid(),
        vec![vec![1.0, 1.0], vec![1.0]], // ragged weights
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap_err();
    assert_eq!(err, BSplineConstructionError::WeightCountMismatch);
}

#[test]
fn test_weight_count_mismatch_wrong_rows() {
    let err = BSplineSurface::new_rational(
        1,
        1,
        ok_grid(),
        vec![vec![1.0, 1.0]], // one row instead of two
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap_err();
    assert_eq!(err, BSplineConstructionError::WeightCountMismatch);
}

#[test]
fn test_non_positive_weight() {
    let err = BSplineSurface::new_rational(
        1,
        1,
        ok_grid(),
        vec![vec![1.0, 0.0], vec![1.0, 1.0]],
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap_err();
    assert_eq!(err, BSplineConstructionError::NonPositiveWeight);

    let err = BSplineSurface::new_rational(
        1,
        1,
        ok_grid(),
        vec![vec![1.0, 1.0], vec![-2.0, 1.0]],
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap_err();
    assert_eq!(err, BSplineConstructionError::NonPositiveWeight);
}

// ---- periodic seam through the public type ----

fn periodic_tube() -> BSplineSurface {
    // Fixture uperiodic_tube: u_degree 2 periodic (6 pole rows), v_degree 1
    // clamped (2 columns).
    let poles = vec![
        vec![Point::new(2.0, 0.0, 0.0), Point::new(2.2, 0.0, 3.0)],
        vec![
            Point::new(1.0000000000000002, 1.7320508075688772, 0.0),
            Point::new(1.1000000000000003, 1.905255888325765, 3.0),
        ],
        vec![
            Point::new(-0.9999999999999996, 1.7320508075688774, 0.0),
            Point::new(-1.0999999999999996, 1.9052558883257653, 3.0),
        ],
        vec![
            Point::new(-2.0, 2.4492935982947064e-16, 0.0),
            Point::new(-2.2, 2.6942229581241775e-16, 3.0),
        ],
        vec![
            Point::new(-1.0000000000000009, -1.7320508075688767, 0.0),
            Point::new(-1.100000000000001, -1.9052558883257646, 3.0),
        ],
        vec![
            Point::new(1.0000000000000002, -1.7320508075688772, 0.0),
            Point::new(1.1000000000000003, -1.905255888325765, 3.0),
        ],
    ];
    BSplineSurface::new(
        2,
        1,
        poles,
        vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        vec![1, 1, 1, 1, 1, 1, 1],
        vec![0.0, 1.0],
        vec![2, 2],
        true,
        false,
    )
    .unwrap()
}

#[test]
fn test_periodic_bounds_and_period() {
    let s = periodic_tube();
    assert_eq!(s.u_bounds(), (0.0, 6.0));
    assert_eq!(s.v_bounds(), (0.0, 1.0));
    assert_eq!(s.u_period(), Some(6.0));
    assert_eq!(s.v_period(), None);
}

#[test]
fn test_periodic_tube_is_u_closed_not_v() {
    let s = periodic_tube();
    assert!(s.is_u_periodic());
    assert!(s.is_u_closed(Tolerance::DEFAULT));
    assert!(s.is_u_closed(Tolerance::DEFAULT));
    assert!(!s.is_v_closed(Tolerance::DEFAULT));
    assert!(!s.is_v_closed(Tolerance::DEFAULT));
}

#[test]
fn test_periodic_seam_matches_across_period() {
    let s = periodic_tube();
    // u = 0 and u = period (6) must give the same point and u-derivative.
    for v in [0.0, 0.3, 0.7, 1.0] {
        let p0 = s.eval_point(0.0, v);
        let p1 = s.eval_point(6.0, v);
        assert!((p0.x - p1.x).abs() < 1e-9, "x seam v={v}");
        assert!((p0.y - p1.y).abs() < 1e-9, "y seam v={v}");
        assert!((p0.z - p1.z).abs() < 1e-9, "z seam v={v}");

        let d0 = s.eval_derivative(0.0, v, 1, 0);
        let d1 = s.eval_derivative(6.0, v, 1, 0);
        assert!((d0.x - d1.x).abs() < 1e-9, "du.x seam v={v}");
        assert!((d0.y - d1.y).abs() < 1e-9, "du.y seam v={v}");
        assert!((d0.z - d1.z).abs() < 1e-9, "du.z seam v={v}");
    }
}

#[test]
fn test_periodic_out_of_window_u_matches_golden() {
    let s = periodic_tube();
    // Golden from fixtures/surfaces_bspline.json (uperiodic_tube), sample
    // with out-of-window u = -0.8.
    let p = s.eval_point(-0.8, 0.6);
    assert!((p.x - 1.7596000000000003).abs() < 1e-7);
    assert!((p.y - -0.5507921568069025).abs() < 1e-7);
    assert!((p.z - 1.7999999999999998).abs() < 1e-7);
}

fn bilinear_patch() -> BSplineSurface {
    let poles = vec![
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 2.0, 0.0)],
        vec![Point::new(2.0, 0.0, 0.0), Point::new(2.0, 2.0, 0.0)],
    ];
    BSplineSurface::new(
        1,
        1,
        poles,
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_bspline_surface_project_point_planar() {
    let s = bilinear_patch();
    let tol = Tolerance::DEFAULT;
    let proj = s.project_point(Point::new(1.0, 1.0, 1.0), tol);
    assert!((proj.u - 0.5).abs() < 1e-9);
    assert!((proj.v - 0.5).abs() < 1e-9);
    assert!((proj.distance - 1.0).abs() < 1e-9);
    // Beyond the edge: boundary foot wins.
    let edge = s.project_point(Point::new(3.0, 1.0, 0.0), tol);
    assert!((edge.u - 1.0).abs() < 1e-9);
    assert!((edge.v - 0.5).abs() < 1e-9);
    assert!((edge.distance - 1.0).abs() < 1e-9);
}

#[test]
fn test_bspline_surface_contains_planar() {
    let s = bilinear_patch();
    let tol = Tolerance::DEFAULT;
    assert!(s.contains(Point::new(1.0, 1.0, 0.0), tol));
    assert!(s.contains(s.eval_point(0.25, 0.75), tol));
    assert!(!s.contains(Point::new(1.0, 1.0, 1.0), tol));
}

#[test]
fn test_bspline_surface_project_point_tube() {
    let s = periodic_tube();
    let tol = Tolerance::DEFAULT;
    // On-surface point projects to (near-)zero distance.
    let on = s.eval_point(0.5, 0.5);
    assert_eq!(s.project_point(on, tol).distance, 0.0);
    assert!(s.contains(on, tol));
    // Off-surface: beat against dense brute-force sampling.
    let query = Point::new(0.0, 0.0, 5.0);
    let proj = s.project_point(query, tol);
    let (u0, u1) = s.u_bounds();
    let (v0, v1) = s.v_bounds();
    let mut best = f64::INFINITY;
    for j in 0..=60 {
        for i in 0..=60 {
            let d = s
                .eval_point(
                    u0 + (u1 - u0) * i as f64 / 60.0,
                    v0 + (v1 - v0) * j as f64 / 60.0,
                )
                .distance(query);
            best = best.min(d);
        }
    }
    assert!(proj.distance <= best + 1e-9);
    assert!(best - proj.distance < 1e-2);
}

// ---- knot insertion ----

use crate::BSplineCurve;

/// A flat row-major net as rows of `row_length`.
fn rows<T: Clone>(net: &[T], row_length: usize) -> Vec<Vec<T>> {
    net.chunks(row_length).map(<[T]>::to_vec).collect()
}

/// Cubic in `u` (periodic: six uniform rows; else clamped with a
/// double knot), quadratic in `v` (periodic nonuniform or clamped),
/// rational if `rational`.
fn insertion_surface(rational: bool, u_periodic: bool, v_periodic: bool) -> BSplineSurface {
    let (u_knots, u_mults) = if u_periodic {
        ((0..7).map(f64::from).collect(), vec![1; 7])
    } else {
        (vec![0.0, 0.8, 2.0, 3.0], vec![4, 1, 2, 4])
    };
    let (v_knots, v_mults) = if v_periodic {
        (vec![0.0, 1.0, 2.5, 3.0, 4.0], vec![1; 5])
    } else {
        (vec![0.0, 1.0, 2.0], vec![3, 1, 3])
    };
    let (n_u, n_v) = (if u_periodic { 6 } else { 7 }, 4);
    let poles: Vec<Vec<Point>> = (0..n_u)
        .map(|i| {
            (0..n_v)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    Point::new(
                        x + 0.3 * y.sin(),
                        y + 0.2 * (x * y).cos(),
                        x.sin() * y.cos(),
                    )
                })
                .collect()
        })
        .collect();
    let weights: Vec<Vec<f64>> = (0..n_u)
        .map(|i| {
            (0..n_v)
                .map(|j| 0.5 + ((i * 7 + j * 3) % 5) as f64 * 0.4)
                .collect()
        })
        .collect();
    if rational {
        BSplineSurface::new_rational(
            3, 2, poles, weights, u_knots, u_mults, v_knots, v_mults, u_periodic, v_periodic,
        )
    } else {
        BSplineSurface::new(
            3, 2, poles, u_knots, u_mults, v_knots, v_mults, u_periodic, v_periodic,
        )
    }
    .unwrap()
}

fn insertion_surfaces() -> Vec<BSplineSurface> {
    [
        (false, false, false),
        (true, false, false),
        (false, true, false),
        (true, false, true),
        (true, true, true),
    ]
    .into_iter()
    .map(|(r, u, v)| insertion_surface(r, u, v))
    .collect()
}

/// `refined` is a valid surface evaluating as `surface`, point and
/// derivatives up to second order, within 1e-12 relative to its size,
/// over a grid including every knot of both.
fn assert_same_surface(surface: &BSplineSurface, refined: &BSplineSurface) {
    let rebuilt = match refined.weights() {
        Some(w) => BSplineSurface::new_rational(
            refined.u_degree(),
            refined.v_degree(),
            rows(refined.poles(), refined.v_pole_count()),
            rows(w, refined.v_pole_count()),
            refined.u_knots().to_vec(),
            refined.u_multiplicities().to_vec(),
            refined.v_knots().to_vec(),
            refined.v_multiplicities().to_vec(),
            refined.is_u_periodic(),
            refined.is_v_periodic(),
        ),
        None => BSplineSurface::new(
            refined.u_degree(),
            refined.v_degree(),
            rows(refined.poles(), refined.v_pole_count()),
            refined.u_knots().to_vec(),
            refined.u_multiplicities().to_vec(),
            refined.v_knots().to_vec(),
            refined.v_multiplicities().to_vec(),
            refined.is_u_periodic(),
            refined.is_v_periodic(),
        ),
    };
    assert!(rebuilt.is_ok(), "invalid refinement: {rebuilt:?}");
    assert_eq!(refined.u_param_range(), surface.u_param_range());
    assert_eq!(refined.v_param_range(), surface.v_param_range());
    let size = surface
        .poles()
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
        .fold(1.0, f64::max);
    let params = |(first, last): (f64, f64), a: &[f64], b: &[f64]| {
        let mut out: Vec<f64> = (0..=12)
            .map(|i| first + (last - first) * i as f64 / 12.0)
            .collect();
        out.extend(a.iter().chain(b));
        out
    };
    let us = params(
        surface.u_param_range(),
        surface.u_knots(),
        refined.u_knots(),
    );
    let vs = params(
        surface.v_param_range(),
        surface.v_knots(),
        refined.v_knots(),
    );
    for &u in &us {
        for &v in &vs {
            let (p, [su, sv], [suu, suv, svv]) = surface.eval_d2(u, v);
            let (q, [tu, tv], [tuu, tuv, tvv]) = refined.eval_d2(u, v);
            let pairs = [
                (p - Point::ORIGIN, q - Point::ORIGIN),
                (su, tu),
                (sv, tv),
                (suu, tuu),
                (suv, tuv),
                (svv, tvv),
            ];
            for (k, (a, b)) in pairs.into_iter().enumerate() {
                let tol = 1e-12 * size.max(a.magnitude());
                assert!(
                    (a - b).magnitude() <= tol,
                    "({u}, {v}), output {k}: {a:?} vs {b:?}"
                );
            }
        }
    }
}

#[test]
fn test_insert_u_and_v_knot_keep_shape() {
    for surface in insertion_surfaces() {
        let (u0, u1) = surface.u_param_range();
        let (v0, v1) = surface.v_param_range();
        for (u, times) in [(u0 + 0.37 * (u1 - u0), 3), (surface.u_knots()[1], 1)] {
            let refined = surface.insert_u_knot(u, times).unwrap();
            assert_eq!(
                refined.u_pole_count(),
                surface.u_pole_count() + times as usize
            );
            assert_eq!(refined.v_knots(), surface.v_knots());
            assert_same_surface(&surface, &refined);
        }
        for (v, times) in [(v0 + 0.37 * (v1 - v0), 2), (surface.v_knots()[1], 1)] {
            let refined = surface.insert_v_knot(v, times).unwrap();
            assert_eq!(refined.u_pole_count(), surface.u_pole_count());
            assert_eq!(
                refined.v_pole_count(),
                surface.v_pole_count() + times as usize
            );
            assert_same_surface(&surface, &refined);
        }
        if surface.is_u_periodic() {
            let refined = surface.insert_u_knot(u1, 2).unwrap();
            assert_eq!(refined.u_multiplicities()[0], 3);
            assert_same_surface(&surface, &refined);
        }
        if surface.is_v_periodic() {
            let refined = surface.insert_v_knot(v0, 1).unwrap();
            assert_eq!(refined.v_multiplicities()[0], 2);
            assert_same_surface(&surface, &refined);
            assert_eq!(
                surface.insert_v_knot(v1, 2),
                Err(KnotInsertionError::MultiplicityTooLarge)
            );
        }
    }
}

#[test]
fn test_bezier_patches_cover_every_span_pair_and_evaluate_as_the_surface() {
    let breaks = |knots: &[f64], (first, last): (f64, f64)| {
        let mut breaks = vec![first];
        breaks.extend(knots.iter().filter(|&&k| k > first && k < last));
        breaks.push(last);
        breaks
    };
    for surface in insertion_surfaces() {
        let patches = surface.to_bezier_patches();
        let us = breaks(surface.u_knots(), surface.u_param_range());
        let vs = breaks(surface.v_knots(), surface.v_param_range());
        assert_eq!(patches.len(), us.len() - 1);
        let (p, q) = (surface.u_degree(), surface.v_degree());
        for (i, row) in patches.iter().enumerate() {
            assert_eq!(row.len(), vs.len() - 1);
            for (j, patch) in row.iter().enumerate() {
                assert_eq!(patch.u_bounds(), (us[i], us[i + 1]));
                assert_eq!(patch.v_bounds(), (vs[j], vs[j + 1]));
                let spline = patch.to_bspline();
                assert_eq!(spline.u_multiplicities(), &[p as u32 + 1; 2]);
                assert_eq!(spline.v_multiplicities(), &[q as u32 + 1; 2]);
                assert_eq!(patch.u_pole_count(), p + 1);
                assert_eq!(patch.v_pole_count(), q + 1);
                assert_eq!(patch.weights().is_some(), surface.weights().is_some());
                for a in 0..=4 {
                    for b in 0..=4 {
                        let u = us[i] + (us[i + 1] - us[i]) * a as f64 / 4.0;
                        let v = vs[j] + (vs[j + 1] - vs[j]) * b as f64 / 4.0;
                        let (x, y) = (surface.eval_point(u, v), patch.eval_point(u, v));
                        assert!(x.distance(y) < 1e-12 * 10.0, "({u}, {v}): {x:?} vs {y:?}");
                    }
                }
            }
        }
    }
}

#[test]
fn test_knot_insertion_keeps_untouched_rational_poles_exactly() {
    let surface = insertion_surface(true, false, false);
    let rows_of = |s: &BSplineSurface| {
        let poles = rows(s.poles(), s.v_pole_count());
        (poles, rows(s.weights().unwrap(), s.v_pole_count()))
    };
    let (poles, weights) = rows_of(&surface);
    // In `u`, 2.5 blends rows 3 to 5 only.
    let (new_poles, new_weights) = rows_of(&surface.insert_u_knot(2.5, 1).unwrap());
    assert_eq!(new_poles[..4], poles[..4]);
    assert_eq!(new_weights[..4], weights[..4]);
    assert_eq!(new_poles[7..], poles[6..]);
    assert_eq!(new_weights[7..], weights[6..]);
    // In `v`, 1.5 blends the poles 1 and 2 of every row only.
    let (new_poles, new_weights) = rows_of(&surface.insert_v_knot(1.5, 1).unwrap());
    for i in 0..poles.len() {
        assert_eq!(new_poles[i][..2], poles[i][..2]);
        assert_eq!(new_weights[i][..2], weights[i][..2]);
        assert_eq!(new_poles[i][4..], poles[i][3..]);
        assert_eq!(new_weights[i][4..], weights[i][3..]);
    }
}

#[test]
fn test_refine_u_and_v_knots_match_curve_refinement() {
    for surface in insertion_surfaces() {
        let (u0, u1) = surface.u_param_range();
        let (v0, v1) = surface.v_param_range();
        let u_insert = [
            u0,
            u0 + 0.2 * (u1 - u0),
            u0 + 0.2 * (u1 - u0),
            u0 + 0.7 * (u1 - u0),
        ];
        let v_insert = [v0 + 0.1 * (v1 - v0), v0 + 0.55 * (v1 - v0), v1];
        let start = usize::from(!surface.is_u_periodic());
        let end = v_insert.len() - usize::from(!surface.is_v_periodic());
        let (u_insert, v_insert) = (&u_insert[start..], &v_insert[..end]);
        let by_u = surface.refine_u_knots(u_insert).unwrap();
        assert_same_surface(&surface, &by_u);
        let by_v = surface.refine_v_knots(v_insert).unwrap();
        assert_same_surface(&surface, &by_v);
        let both = by_u.refine_v_knots(v_insert).unwrap();
        assert_same_surface(&surface, &both);
        // Each column (row) is refined as a curve with the `u` (`v`)
        // knots would be.
        let line = |poles: Vec<Point>, weights: Option<Vec<f64>>, u: bool| {
            let (degree, knots, mults, periodic) = if u {
                (
                    3,
                    surface.u_knots(),
                    surface.u_multiplicities(),
                    surface.is_u_periodic(),
                )
            } else {
                (
                    2,
                    surface.v_knots(),
                    surface.v_multiplicities(),
                    surface.is_v_periodic(),
                )
            };
            let (knots, mults) = (knots.to_vec(), mults.to_vec());
            match weights {
                Some(w) => BSplineCurve::new_rational(degree, poles, w, knots, mults, periodic),
                None => BSplineCurve::new(degree, poles, knots, mults, periodic),
            }
            .unwrap()
        };
        let (net, weights) = (
            rows(surface.poles(), surface.v_pole_count()),
            surface.weights().map(|w| rows(w, surface.v_pole_count())),
        );
        let by_u_net = rows(by_u.poles(), by_u.v_pole_count());
        for j in 0..surface.v_pole_count() {
            let column = net.iter().map(|row| row[j]).collect();
            let column_weights = weights
                .as_ref()
                .map(|w| w.iter().map(|row| row[j]).collect());
            let curve = line(column, column_weights, true)
                .refine_knots(u_insert)
                .unwrap();
            for (row, q) in by_u_net.iter().zip(curve.poles()) {
                assert!(row[j].distance(*q) < 1e-12, "{:?} vs {q:?}", row[j]);
            }
        }
        let by_v_net = rows(by_v.poles(), by_v.v_pole_count());
        for (i, row) in net.iter().enumerate() {
            let row_weights = weights.as_ref().map(|w| w[i].clone());
            let curve = line(row.clone(), row_weights, false)
                .refine_knots(v_insert)
                .unwrap();
            for (p, q) in by_v_net[i].iter().zip(curve.poles()) {
                assert!(p.distance(*q) < 1e-12, "{p:?} vs {q:?}");
            }
        }
    }
}

#[test]
fn test_surface_knot_insertion_errors() {
    let surface = insertion_surface(true, false, true);
    assert_eq!(
        surface.insert_u_knot(f64::NAN, 1),
        Err(KnotInsertionError::NonFiniteKnot)
    );
    assert_eq!(
        surface.insert_v_knot(-0.5, 1),
        Err(KnotInsertionError::KnotOutOfRange)
    );
    assert_eq!(
        surface.refine_u_knots(&[1.0, 0.5]),
        Err(KnotInsertionError::KnotsNotSorted)
    );
    assert_eq!(
        surface.refine_v_knots(&[0.5, 0.25]),
        Err(KnotInsertionError::KnotsNotSorted)
    );
    // The clamped `u` ends are full; the `u` double knot takes one more.
    assert_eq!(
        surface.insert_u_knot(0.0, 1),
        Err(KnotInsertionError::MultiplicityTooLarge)
    );
    assert_eq!(
        surface.insert_u_knot(2.0, 2),
        Err(KnotInsertionError::MultiplicityTooLarge)
    );
    assert_eq!(
        surface.insert_v_knot(1.0, 2),
        Err(KnotInsertionError::MultiplicityTooLarge)
    );
    assert_eq!(surface.insert_u_knot(1.0, 0).unwrap(), surface);
    assert_eq!(surface.refine_v_knots(&[]).unwrap(), surface);
}

// ---- knot removal ----

fn net_size(surface: &BSplineSurface) -> f64 {
    surface
        .poles()
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
        .fold(1.0, f64::max)
}

/// `simpler` has the knots and multiplicities of `surface`, and its
/// poles and weights within 1e-12 relative.
fn assert_same_net(surface: &BSplineSurface, simpler: &BSplineSurface) {
    assert_eq!(simpler.u_knots(), surface.u_knots());
    assert_eq!(simpler.u_multiplicities(), surface.u_multiplicities());
    assert_eq!(simpler.v_knots(), surface.v_knots());
    assert_eq!(simpler.v_multiplicities(), surface.v_multiplicities());
    let size = net_size(surface);
    assert_eq!(simpler.poles().len(), surface.poles().len());
    for (p, q) in simpler.poles().iter().zip(surface.poles()) {
        assert!(p.distance(*q) <= 1e-12 * size, "{p:?} vs {q:?}");
    }
    if let (Some(a), Some(b)) = (simpler.weights(), surface.weights()) {
        for (x, y) in a.iter().zip(b) {
            assert!((x - y).abs() <= 1e-12 * y, "{x} vs {y}");
        }
    }
    assert_eq!(simpler.weights().is_some(), surface.weights().is_some());
}

/// Largest distance between the points of `a` and `b` at the same
/// parameters, over a 61x61 grid and every knot of both.
fn sampled_surface_deviation(a: &BSplineSurface, b: &BSplineSurface) -> f64 {
    assert_eq!(a.u_param_range(), b.u_param_range());
    assert_eq!(a.v_param_range(), b.v_param_range());
    let params = |(first, last): (f64, f64), x: &[f64], y: &[f64]| {
        (0..=60)
            .map(|i| first + (last - first) * i as f64 / 60.0)
            .chain(x.iter().chain(y).copied())
            .filter(|t| (first..=last).contains(t))
            .collect::<Vec<f64>>()
    };
    let us = params(a.u_param_range(), a.u_knots(), b.u_knots());
    let vs = params(a.v_param_range(), a.v_knots(), b.v_knots());
    let mut deviation: f64 = 0.0;
    for &u in &us {
        for &v in &vs {
            deviation = deviation.max(a.eval_point(u, v).distance(b.eval_point(u, v)));
        }
    }
    deviation
}

/// `surface` with its poles replaced by `poles`.
fn with_net(surface: &BSplineSurface, poles: Vec<Vec<Point>>) -> BSplineSurface {
    let (uk, um) = (
        surface.u_knots().to_vec(),
        surface.u_multiplicities().to_vec(),
    );
    let (vk, vm) = (
        surface.v_knots().to_vec(),
        surface.v_multiplicities().to_vec(),
    );
    let (up, vp) = (surface.is_u_periodic(), surface.is_v_periodic());
    let (du, dv) = (surface.u_degree(), surface.v_degree());
    match surface.weights() {
        Some(w) => {
            let w = rows(w, surface.v_pole_count());
            BSplineSurface::new_rational(du, dv, poles, w, uk, um, vk, vm, up, vp)
        }
        None => BSplineSurface::new(du, dv, poles, uk, um, vk, vm, up, vp),
    }
    .unwrap()
}

#[test]
fn test_remove_inserted_u_and_v_knots_give_back_the_surface() {
    for surface in insertion_surfaces() {
        let tol = 1e-9 * net_size(&surface);
        let (u0, u1) = surface.u_param_range();
        let (v0, v1) = surface.v_param_range();
        for (u, times) in [(u0 + 0.37 * (u1 - u0), 3), (surface.u_knots()[1], 1)] {
            let refined = surface.insert_u_knot(u, times).unwrap();
            let (simpler, removed) = refined.remove_u_knot(u, times, tol).unwrap();
            assert_eq!(removed, times);
            assert_same_net(&surface, &simpler);
        }
        for (v, times) in [(v0 + 0.37 * (v1 - v0), 2), (surface.v_knots()[1], 1)] {
            let refined = surface.insert_v_knot(v, times).unwrap();
            let (simpler, removed) = refined.remove_v_knot(v, times, tol).unwrap();
            assert_eq!(removed, times);
            assert_same_net(&surface, &simpler);
        }
        if surface.is_u_periodic() {
            let refined = surface.insert_u_knot(u0, 2).unwrap();
            let (simpler, removed) = refined.remove_u_knot(u1, 2, tol).unwrap();
            assert_eq!(removed, 2);
            assert_same_net(&surface, &simpler);
        }
        if surface.is_v_periodic() {
            let refined = surface.insert_v_knot(v1, 1).unwrap();
            let (simpler, removed) = refined.remove_v_knot(v0, 1, tol).unwrap();
            assert_eq!(removed, 1);
            assert_same_net(&surface, &simpler);
        }
    }
}

#[test]
fn test_remove_u_and_v_knots_undo_a_refinement() {
    for surface in insertion_surfaces() {
        let tol = 1e-9 * net_size(&surface);
        let (u0, u1) = surface.u_param_range();
        let (v0, v1) = surface.v_param_range();
        let u_knots = [
            u0 + 0.2 * (u1 - u0),
            u0 + 0.2 * (u1 - u0),
            u0 + 0.7 * (u1 - u0),
        ];
        let v_knots = [v0 + 0.1 * (v1 - v0), v0 + 0.55 * (v1 - v0)];
        let both = surface
            .refine_u_knots(&u_knots)
            .unwrap()
            .refine_v_knots(&v_knots)
            .unwrap();
        let (by_u, removed) = both.remove_u_knots(&u_knots, tol).unwrap();
        assert_eq!(removed, 3);
        let (simpler, removed) = by_u.remove_v_knots(&v_knots, tol).unwrap();
        assert_eq!(removed, 2);
        assert_same_net(&surface, &simpler);
    }
}

#[test]
fn test_surface_removal_keeps_knots_the_shape_needs() {
    for surface in insertion_surfaces() {
        let tol = 1e-6 * net_size(&surface);
        let (u0, u1) = surface.u_param_range();
        let (v0, v1) = surface.v_param_range();
        let interior = |knots: &[f64], (first, last): (f64, f64), periodic: bool| {
            knots
                .iter()
                .copied()
                .filter(|&k| {
                    if periodic {
                        k < last
                    } else {
                        first < k && k < last
                    }
                })
                .collect::<Vec<f64>>()
        };
        let us = interior(surface.u_knots(), (u0, u1), surface.is_u_periodic());
        let vs = interior(surface.v_knots(), (v0, v1), surface.is_v_periodic());
        let (same, removed) = surface.remove_u_knots(&us, tol).unwrap();
        assert_eq!(removed, 0);
        assert_eq!(same, surface);
        let (same, removed) = surface.remove_v_knots(&vs, tol).unwrap();
        assert_eq!(removed, 0);
        assert_eq!(same, surface);
    }
}

#[test]
fn test_surface_removal_needs_every_column_and_row() {
    // One pole moved off an inserted knot: the knot stays in that
    // direction until the tolerance covers the move.
    for surface in insertion_surfaces() {
        let size = net_size(&surface);
        let shift = 1e-4 * size;
        let (u0, u1) = surface.u_param_range();
        let u = u0 + 0.43 * (u1 - u0);
        let refined = surface.insert_u_knot(u, 1).unwrap();
        let mut poles = rows(refined.poles(), refined.v_pole_count());
        let (i, j) = (poles.len() / 2, poles[0].len() - 1);
        poles[i][j] = poles[i][j] + Vector::new(shift, 0.0, -shift);
        let moved = with_net(&refined, poles);
        assert_eq!(moved.remove_u_knot(u, 1, 1e-2 * shift).unwrap().1, 0);
        let tol = 100.0 * shift;
        let (simpler, removed) = moved.remove_u_knot(u, 1, tol).unwrap();
        assert_eq!(removed, 1);
        assert!(sampled_surface_deviation(&moved, &simpler) <= tol);

        let (v0, v1) = surface.v_param_range();
        let v = v0 + 0.61 * (v1 - v0);
        let refined = surface.insert_v_knot(v, 1).unwrap();
        let mut poles = rows(refined.poles(), refined.v_pole_count());
        let (i, j) = (poles.len() - 1, poles[0].len() / 2);
        poles[i][j] = poles[i][j] + Vector::new(-shift, shift, shift);
        let moved = with_net(&refined, poles);
        assert_eq!(moved.remove_v_knot(v, 1, 1e-2 * shift).unwrap().1, 0);
        let (simpler, removed) = moved.remove_v_knot(v, 1, tol).unwrap();
        assert_eq!(removed, 1);
        assert!(sampled_surface_deviation(&moved, &simpler) <= tol);
    }
}

#[test]
fn test_surface_removal_stays_within_max_deviation() {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut removed_all = 0;
    for surface in insertion_surfaces() {
        let (u0, u1) = surface.u_param_range();
        let (v0, v1) = surface.v_param_range();
        let us: Vec<f64> = [0.15, 0.45, 0.8]
            .iter()
            .map(|t| u0 + t * (u1 - u0))
            .collect();
        let vs: Vec<f64> = [0.3, 0.65].iter().map(|t| v0 + t * (v1 - v0)).collect();
        let refined = surface
            .refine_u_knots(&us)
            .unwrap()
            .refine_v_knots(&vs)
            .unwrap();
        let noise = 1e-3 * net_size(&refined);
        let poles = rows(refined.poles(), refined.v_pole_count())
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&p| p + Vector::new(next() - 0.5, next() - 0.5, next() - 0.5) * noise)
                    .collect()
            })
            .collect();
        let noisy = with_net(&refined, poles);
        let mut u_all = us.clone();
        u_all.extend(surface.u_knots().iter().filter(|&&k| u0 < k && k < u1));
        let mut v_all = vs.clone();
        v_all.extend(surface.v_knots().iter().filter(|&&k| v0 < k && k < v1));
        for factor in [0.3, 3.0, 30.0] {
            let tol = factor * noise;
            let (simpler, removed) = noisy.remove_u_knots(&u_all, tol).unwrap();
            assert!(sampled_surface_deviation(&noisy, &simpler) <= tol);
            removed_all += removed;
            let (simpler, removed) = noisy.remove_v_knots(&v_all, tol).unwrap();
            assert!(sampled_surface_deviation(&noisy, &simpler) <= tol);
            removed_all += removed;
        }
    }
    assert!(removed_all > 10, "{removed_all}");
}

#[test]
fn test_surface_knot_removal_errors() {
    let surface = insertion_surface(true, false, true);
    assert_eq!(
        surface.remove_u_knot(0.8, 1, -1.0),
        Err(KnotRemovalError::InvalidDeviation)
    );
    assert_eq!(
        surface.remove_v_knots(&[], f64::NAN),
        Err(KnotRemovalError::InvalidDeviation)
    );
    assert_eq!(
        surface.remove_u_knot(f64::NAN, 1, 1e-3),
        Err(KnotRemovalError::NonFiniteKnot)
    );
    assert_eq!(
        surface.remove_v_knots(&[1.0, f64::INFINITY], 1e-3),
        Err(KnotRemovalError::NonFiniteKnot)
    );
    // The clamped `u` ends stay; the periodic `v` seam may be named.
    assert_eq!(
        surface.remove_u_knot(0.0, 1, 1e-3),
        Err(KnotRemovalError::KnotOutOfRange)
    );
    assert_eq!(
        surface.remove_u_knots(&[0.8, 3.0], 1e-3),
        Err(KnotRemovalError::KnotOutOfRange)
    );
    assert_eq!(
        surface.remove_v_knot(4.5, 1, 1e-3),
        Err(KnotRemovalError::KnotOutOfRange)
    );
    assert_eq!(surface.remove_v_knot(4.0, 1, 1e-3).unwrap().1, 0);
    assert_eq!(
        surface.remove_u_knot(1.5, 1, 1e-3),
        Err(KnotRemovalError::KnotNotFound)
    );
    assert_eq!(
        surface.remove_v_knots(&[1.0, 2.0], 1e-3),
        Err(KnotRemovalError::KnotNotFound)
    );
    assert_eq!(
        surface.remove_u_knot(0.8, 0, 1e-3).unwrap(),
        (surface.clone(), 0)
    );
    assert_eq!(
        surface.remove_v_knots(&[], 1e-3).unwrap(),
        (surface.clone(), 0)
    );
}

#[test]
fn test_rational_surface_removal_bounds_every_column_about_one_centre() {
    // One column of the net sits at the origin with an uneven weight,
    // the other far away: the weight left over by a removal moves no
    // point of the first column, but drags the points in between.
    let column = |x: f64| vec![Point::new(x, 0.0, 0.0); 4];
    let (near, far) = (column(0.0), column(1000.0));
    let poles: Vec<Vec<Point>> = (0..4).map(|i| vec![near[i], far[i]]).collect();
    let weights: Vec<Vec<f64>> = [1.0, 1.0, 1.004, 1.0]
        .iter()
        .map(|&w| vec![w, 1.0])
        .collect();
    let surface = BSplineSurface::new_rational(
        2,
        1,
        poles,
        weights,
        vec![0.0, 0.5, 1.0],
        vec![3, 1, 3],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap();
    let (_, removed) = surface.remove_u_knot(0.5, 1, 0.01).unwrap();
    assert_eq!(removed, 0);
    for tol in [0.01, 1.0, 10.0] {
        let (simpler, _) = surface.remove_u_knot(0.5, 1, tol).unwrap();
        assert!(sampled_surface_deviation(&surface, &simpler) <= tol);
    }
}

// ---- transformed ----

fn transformed_test_patch() -> BSplineSurface {
    let poles = vec![
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
        vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    ];
    BSplineSurface::new(
        1,
        1,
        poles,
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_bspline_surface_transformed_translation_spot_check() {
    let surface = transformed_test_patch();
    let t = Transform::translation(Vector::new(0.0, 0.0, 4.0));
    let moved = surface.transformed(t);
    // Hand-computed patch midpoint, shifted.
    assert_eq!(moved.eval_point(0.5, 0.5), Point::new(0.5, 0.5, 4.25));
    assert_eq!(moved.u_degree(), surface.u_degree());
    assert_eq!(moved.v_degree(), surface.v_degree());
}

#[test]
fn test_bspline_surface_transformed_round_trip_and_mirror() {
    use crate::{Axis, Frame};
    use std::f64::consts::PI;
    let surface = transformed_test_patch();
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let t = Transform::rotation(axis, PI / 2.0).then(Transform::scaling(Point::ORIGIN, 2.0));
    let moved = surface.transformed(t);
    let tol = Tolerance::DEFAULT;
    for (u, v) in [(0.0, 0.0), (0.5, 0.5), (1.0, 0.25)] {
        assert!(moved.contains(t.apply_point(surface.eval_point(u, v)), tol));
    }
    let inv = Transform::scaling(Point::ORIGIN, 0.5).then(Transform::rotation(axis, -PI / 2.0));
    let back = moved.transformed(inv);
    for (u, v) in [(0.0, 0.0), (0.5, 0.5), (1.0, 0.25)] {
        assert!(back.eval_point(u, v).distance(surface.eval_point(u, v)) < 1e-9);
    }
    // Mirror preserves the point set (as a fresh pole net).
    let m = Transform::mirror_plane(Frame::WORLD);
    let mirrored = surface.transformed(m);
    for (u, v) in [(0.0, 0.0), (0.5, 0.5), (1.0, 0.25)] {
        assert!(mirrored.contains(m.apply_point(surface.eval_point(u, v)), tol));
    }
    assert_eq!(mirrored.transformed(m), surface);
}

// ---- exact projection ----

/// A Bernstein net evaluated by de Casteljau, `[rows][cols][h]`.
fn bernstein_point(net: &[f64], rows: usize, cols: usize, h: usize, s: f64, t: f64) -> Point {
    let lerp = |a: &[f64], n: usize, x: f64| {
        let mut a = a.to_vec();
        for r in 1..n {
            for k in 0..(n - r) * h {
                a[k] = (1.0 - x) * a[k] + x * a[k + h];
            }
        }
        a[..h].to_vec()
    };
    let column: Vec<f64> = (0..rows)
        .flat_map(|i| lerp(&net[i * cols * h..(i + 1) * cols * h], cols, t))
        .collect();
    let p = lerp(&column, rows, s);
    let w = if h == 4 { p[3] } else { 1.0 };
    Point::new(p[0] / w, p[1] / w, p[2] / w)
}

/// An extruded degree-6 profile `z(x)`, `x` linear in `u`, with two
/// dips of nearly equal depth, the right one deeper: one knot span.
fn two_basins() -> BSplineSurface {
    let z = [1.0, 1.0, -1.0, 1.2, -1.04, 1.0, 1.0];
    let poles = (0..7)
        .map(|i| {
            (0..2)
                .map(|j| Point::new(i as f64 / 6.0, j as f64, z[i]))
                .collect()
        })
        .collect();
    BSplineSurface::new(
        6,
        1,
        poles,
        vec![0.0, 1.0],
        vec![7, 7],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

/// The distance from `point` to the profile of [`two_basins`] at
/// `v = 0.5`, by dense sampling then golden-section refinement.
fn profile_distance(s: &BSplineSurface, point: Point) -> f64 {
    let d = |u: f64| s.eval_point(u.clamp(0.0, 1.0), 0.5).distance(point);
    let n = 20_000;
    let best = (0..=n)
        .map(|i| i as f64 / n as f64)
        .min_by(|a, b| d(*a).total_cmp(&d(*b)))
        .unwrap();
    let (mut a, mut b) = (best - 1.0 / n as f64, best + 1.0 / n as f64);
    let g = 0.5 * (5f64.sqrt() - 1.0);
    for _ in 0..80 {
        let (x, y) = (b - g * (b - a), a + g * (b - a));
        if d(x) < d(y) {
            b = y;
        } else {
            a = x;
        }
    }
    d(0.5 * (a + b))
}

#[test]
fn test_project_point_two_basins_in_one_patch() {
    let s = two_basins();
    let tol = Tolerance::DEFAULT;
    // Points above and between the dips, some nearer the shallower one.
    for x in [0.38, 0.42, 0.45, 0.48, 0.5, 0.55] {
        let query = Point::new(x, 0.5, 0.6);
        let proj = s.project_point(query, tol);
        let reference = profile_distance(&s, query);
        assert!(
            (proj.distance - reference).abs() < 1e-12,
            "x = {x}: {} vs {reference}",
            proj.distance
        );
        // The foot is a critical point: the residual is normal there.
        let (p, su, sv) = s.eval_surface_d1(proj.u, proj.v);
        let r = p - query;
        assert!(r.dot(su).abs() < 1e-9 * su.magnitude() * r.magnitude());
        assert!(r.dot(sv).abs() < 1e-9 * sv.magnitude() * r.magnitude());
        let extrema = s.extrema(query, tol);
        assert_eq!(extrema[0].distance, proj.distance);
        // Both dips are local minima.
        assert!(extrema.len() >= 2, "x = {x}: {extrema:?}");
    }
}

/// A folded net of degree 2 by 5, three `u` spans and two `v` spans,
/// whose sheets pass close to each other: many points have a second
/// local minimum a short distance away in the same patch.
fn folded_net() -> BSplineSurface {
    let row = |p: &[(f64, f64, f64)]| -> Vec<Point> {
        p.iter().map(|&(x, y, z)| Point::new(x, y, z)).collect()
    };
    #[rustfmt::skip]
    let poles = vec![
        row(&[(4.52, -4.72, -0.97), (0.94, -2.11, -2.98), (0.99, -4.89, 0.49), (4.03, 2.02, 2.08), (2.35, -4.38, -1.07),
              (-2.22, -2.8, -3.92), (-4.6, 2.81, -3.13), (3.21, 2.17, 2.1), (0.74, 0.06, -3.3), (3.03, 3.51, 1.31)]),
        row(&[(2.04, -4.83, 3.15), (-4.97, 1.53, 2.88), (0.8, 4.4, 3.79), (-3.33, -3.15, -4.61), (3.81, -4.51, 1.64),
              (0.84, -1.12, -4.38), (-1.41, -3.42, -4.87), (-1.57, -3.53, 0.75), (-4.0, 4.32, -1.94), (-2.98, 4.31, -3.09)]),
        row(&[(-4.47, 0.45, 0.78), (0.81, -0.98, 2.47), (4.44, 2.43, 2.29), (-1.33, -4.42, 0.08), (1.79, -2.52, -1.4),
              (-0.55, -4.36, -4.45), (-4.85, 3.25, -2.06), (0.55, 3.93, -4.22), (2.3, 0.75, 0.34), (2.87, 3.87, -3.4)]),
        row(&[(-1.7, 2.4, -2.3), (-0.07, -0.04, 1.97), (-4.41, -0.31, -0.26), (0.5, 1.27, 2.44), (4.69, -2.87, -4.96),
              (3.99, -1.46, 4.44), (-4.07, 2.78, -0.86), (3.23, -2.36, -2.98), (1.69, -2.35, 4.7), (-3.49, -0.13, -0.89)]),
        row(&[(4.62, -2.18, -3.87), (4.1, -0.82, 0.15), (3.88, -2.59, 4.99), (0.8, 3.99, -0.12), (-3.78, 1.26, -0.58),
              (3.27, -1.56, 0.59), (-0.55, 0.41, -2.35), (-2.7, 3.86, 4.63), (-1.19, -3.46, 2.63), (-1.62, -1.04, 1.81)]),
    ];
    BSplineSurface::new(
        2,
        5,
        poles,
        vec![0.0, 1.0, 2.0, 3.0],
        vec![3, 1, 1, 3],
        vec![0.0, 0.56, 1.0],
        vec![6, 4, 6],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_project_point_close_basins_on_a_folded_net() {
    let s = folded_net();
    let tol = Tolerance::DEFAULT;
    // In patch `[0, 1] x [0, 0.56]`: the foot is the point itself, and a
    // second minimum lies on the edge `u = 0` about 0.08 away.
    let on = s.eval_point(0.05, 0.1);
    let proj = s.project_point(on, tol);
    assert_eq!(proj.distance, 0.0);
    let extrema = s.extrema(on, tol);
    assert_eq!(extrema[0].distance, 0.0);
    assert!(
        extrema[1..]
            .iter()
            .any(|e| e.u == 0.0 && e.distance > 0.05 && e.distance < 0.1),
        "{extrema:?}"
    );
    // Every point of the surface projects onto itself.
    let n = 40;
    for i in 0..=n {
        for j in 0..=n {
            let (u, v) = (3.0 * i as f64 / n as f64, j as f64 / n as f64);
            let on = s.eval_point(u, v);
            let proj = s.project_point(on, tol);
            assert!(proj.distance < 1e-12, "({u}, {v}): {proj:?}");
        }
    }
}

#[test]
fn test_project_point_on_surface_is_zero() {
    let tube = periodic_tube();
    let (u0, u1) = tube.u_bounds();
    let tol = Tolerance::DEFAULT;
    for s in [two_basins(), tube, rational_quarter_cylinder(2.0)] {
        let (u0, u1) = if s.is_u_periodic() {
            (u0, u1)
        } else {
            s.u_bounds()
        };
        let (v0, v1) = s.v_bounds();
        for i in 0..=7 {
            for j in 0..=5 {
                let (u, v) = (
                    u0 + (u1 - u0) * i as f64 / 7.0,
                    v0 + (v1 - v0) * j as f64 / 5.0,
                );
                let on = s.eval_point(u, v);
                let proj = s.project_point(on, tol);
                assert_eq!(proj.distance, 0.0);
                assert!(s.eval_point(proj.u, proj.v).distance(on) < 1e-12);
            }
        }
    }
}

/// An exact quarter cylinder of radius `r` about the z axis: a rational
/// quadratic arc in `u`, extruded over `z` in `[0, 1]` in `v`.
fn rational_quarter_cylinder(r: f64) -> BSplineSurface {
    let arc = [(r, 0.0), (r, r), (0.0, r)];
    let poles = arc
        .iter()
        .map(|&(x, y)| vec![Point::new(x, y, 0.0), Point::new(x, y, 1.0)])
        .collect();
    let w = std::f64::consts::FRAC_1_SQRT_2;
    BSplineSurface::new_rational(
        2,
        1,
        poles,
        vec![vec![1.0, 1.0], vec![w, w], vec![1.0, 1.0]],
        vec![0.0, 1.0],
        vec![3, 3],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_project_point_exact_foot_on_rational_cylinder() {
    let radius = 2.0;
    let s = rational_quarter_cylinder(radius);
    let tol = Tolerance::DEFAULT;
    for (rho, theta, z) in [(0.5, 0.3, 0.2), (3.0, 1.1, 0.7), (7.5, 0.05, 0.5)] {
        let query = Point::new(rho * f64::cos(theta), rho * f64::sin(theta), z);
        let proj = s.project_point(query, tol);
        assert!((proj.distance - (rho - radius).abs()).abs() < 1e-12);
        let foot = s.eval_point(proj.u, proj.v);
        let expected = Point::new(radius * f64::cos(theta), radius * f64::sin(theta), z);
        assert!(foot.distance(expected) < 1e-12);
    }
}

#[test]
fn test_project_point_exact_foot_on_edge_corner_and_crease() {
    let tol = Tolerance::DEFAULT;
    // Bilinear square [0, 2]^2 in z = 0 over [0, 1]^2: an edge foot
    // and a corner foot.
    let s = bilinear_patch();
    let proj = s.project_point(Point::new(3.0, 1.0, 1.0), tol);
    assert_eq!((proj.u, proj.v), (1.0, 0.5));
    assert!((proj.distance - 2f64.sqrt()).abs() < 1e-12);
    let proj = s.project_point(Point::new(3.0, 3.0, 1.0), tol);
    assert_eq!((proj.u, proj.v), (1.0, 1.0));
    assert!((proj.distance - 3f64.sqrt()).abs() < 1e-12);
    // A ridge: degree 1 across a crease knot, the foot on the crease.
    let ridge = BSplineSurface::new(
        1,
        1,
        vec![
            vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
            vec![Point::new(1.0, 0.0, 1.0), Point::new(1.0, 1.0, 1.0)],
            vec![Point::new(2.0, 0.0, 0.0), Point::new(2.0, 1.0, 0.0)],
        ],
        vec![0.0, 1.0, 2.0],
        vec![2, 1, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap();
    let proj = ridge.project_point(Point::new(1.0, 0.25, 3.0), tol);
    assert_eq!(proj.u, 1.0);
    assert!((proj.v - 0.25).abs() < 1e-12);
    assert!((proj.distance - 2.0).abs() < 1e-12);
}

#[test]
fn test_patch_nets_match_the_surface_and_halve() {
    for s in [
        two_basins(),
        periodic_tube(),
        rational_quarter_cylinder(1.5),
    ] {
        let patches = s.patch_tree();
        let mut search = PatchSearch::new(&s, Point::ORIGIN, 0.0, true);
        let (rows, cols, h) = (s.u_degree() + 1, s.v_degree() + 1, s.net_stride());
        for k in 0..patches.len() {
            let (rect, _) = patches.patch(k);
            let root = search.root(patches, k).unwrap();
            let net = search.nets[root.slot..root.slot + search.size].to_vec();
            let mut halves = vec![vec![0.0; net.len()]; 2];
            let [left, right] = &mut halves[..] else {
                unreachable!()
            };
            math::halve_net(&net, rows, cols, h, true, left, right);
            for (a, b) in [(0.0, 0.0), (0.3, 0.8), (1.0, 0.5), (0.7, 1.0)] {
                let (u, v) = (
                    rect[0] + a * (rect[1] - rect[0]),
                    rect[2] + b * (rect[3] - rect[2]),
                );
                let expected = s.eval_point(u, v);
                let got = bernstein_point(&net, rows, cols, h, a, b);
                assert!(got.distance(expected) < 1e-12);
                // The halves along `u` cover `[0, 1/2]` and `[1/2, 1]`.
                let half = if a < 0.5 {
                    bernstein_point(left, rows, cols, h, 2.0 * a, b)
                } else {
                    bernstein_point(right, rows, cols, h, 2.0 * a - 1.0, b)
                };
                assert!(half.distance(expected) < 1e-12);
            }
        }
    }
}

// ---- split, extract, reverse ----

/// The shortest nonzero knot span of `knots`.
fn shortest_span(knots: &[f64]) -> f64 {
    knots
        .windows(2)
        .map(|w| w[1] - w[0])
        .fold(f64::INFINITY, f64::min)
}

/// The largest pole coordinate of `surface`, at least 1.
fn size(surface: &BSplineSurface) -> f64 {
    surface
        .poles()
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
        .fold(1.0, f64::max)
}

/// `piece` is a valid surface evaluating as `surface` over its bounds,
/// points and derivatives up to second order on a dense grid, points
/// on a grid of its knots. A direction given a range is clamped,
/// non-periodic and bounded exactly by it; the others are unchanged.
/// Derivative tolerances grow as in the curve tests where the piece's
/// spans are shorter than the surface's.
fn assert_patch(
    surface: &BSplineSurface,
    piece: &BSplineSurface,
    u_range: Option<(f64, f64)>,
    v_range: Option<(f64, f64)>,
) {
    assert!(
        BSplineSurface::new(
            piece.u_degree(),
            piece.v_degree(),
            rows(piece.poles(), piece.v_pole_count()),
            piece.u_knots().to_vec(),
            piece.u_multiplicities().to_vec(),
            piece.v_knots().to_vec(),
            piece.v_multiplicities().to_vec(),
            piece.is_u_periodic(),
            piece.is_v_periodic(),
        )
        .is_ok()
    );
    assert_eq!(piece.is_rational(), surface.is_rational());
    let check = |range: Option<(f64, f64)>,
                 bounds: (f64, f64),
                 periodic: bool,
                 degree: usize,
                 mults: &[u32],
                 original: ((f64, f64), bool)| {
        match range {
            Some(range) => {
                assert_eq!(bounds, range);
                assert!(!periodic);
                let end = degree as u32 + 1;
                assert_eq!((mults[0], mults[mults.len() - 1]), (end, end));
            }
            None => assert_eq!((bounds, periodic), original),
        }
    };
    check(
        u_range,
        piece.u_param_range(),
        piece.is_u_periodic(),
        piece.u_degree(),
        piece.u_multiplicities(),
        (surface.u_param_range(), surface.is_u_periodic()),
    );
    check(
        v_range,
        piece.v_param_range(),
        piece.is_v_periodic(),
        piece.v_degree(),
        piece.v_multiplicities(),
        (surface.v_param_range(), surface.is_v_periodic()),
    );
    let size = size(surface);
    let stretch = (shortest_span(surface.u_knots()) / shortest_span(piece.u_knots()))
        .max(shortest_span(surface.v_knots()) / shortest_span(piece.v_knots()))
        .max(1.0);
    let grid = |(first, last): (f64, f64)| -> Vec<f64> {
        (0..16)
            .map(|i| first + (last - first) * (i as f64 + 0.5) / 16.0)
            .collect()
    };
    let (us, vs) = (grid(piece.u_param_range()), grid(piece.v_param_range()));
    for &u in &us {
        for &v in &vs {
            let (p, [su, sv], [suu, suv, svv]) = surface.eval_d2(u, v);
            let (q, [tu, tv], [tuu, tuv, tvv]) = piece.eval_d2(u, v);
            let pairs = [
                (p - Point::ORIGIN, q - Point::ORIGIN, 1.0),
                (su, tu, stretch),
                (sv, tv, stretch),
                (suu, tuu, stretch * stretch),
                (suv, tuv, stretch * stretch),
                (svv, tvv, stretch * stretch),
            ];
            for (k, (a, b, scale)) in pairs.into_iter().enumerate() {
                let tol = 1e-12 * (size * scale).max(a.magnitude());
                assert!(
                    (a - b).magnitude() <= tol,
                    "({u}, {v}), output {k}: {a:?} vs {b:?}"
                );
            }
        }
    }
    for &u in piece.u_knots() {
        for &v in piece.v_knots() {
            let (p, q) = (surface.eval_point(u, v), piece.eval_point(u, v));
            assert!(p.distance(q) <= 1e-12 * size, "({u}, {v}): {p:?} vs {q:?}");
        }
    }
}

/// Ranges of a direction with the given bounds and knots: a middle
/// one, the whole, from and to an interior knot and, if periodic,
/// across the seam and over a whole period from a knot.
fn ranges((first, last): (f64, f64), knots: &[f64], periodic: bool) -> Vec<(f64, f64)> {
    let at = |t: f64| first + t * (last - first);
    let mut out = vec![(at(0.13), at(0.71)), (first, last)];
    if let Some(&k) = knots.iter().find(|&&k| k > first && k < last) {
        out.extend([(k, at(0.9)), (at(0.05), k)]);
        if periodic {
            out.push((k, k + (last - first)));
        }
    }
    if periodic {
        out.extend([(at(0.6), at(1.3)), (at(-0.4), at(0.2))]);
    }
    out
}

#[test]
fn test_extract_sub_patches_evaluate_as_the_surface() {
    for surface in insertion_surfaces() {
        let u_ranges = ranges(
            surface.u_param_range(),
            surface.u_knots(),
            surface.is_u_periodic(),
        );
        let v_ranges = ranges(
            surface.v_param_range(),
            surface.v_knots(),
            surface.is_v_periodic(),
        );
        for &u in &u_ranges {
            let piece = surface.extract_u(u.0, u.1).unwrap();
            assert_patch(&surface, &piece, Some(u), None);
            for &v in &v_ranges {
                let piece = surface.extract(u, v).unwrap();
                assert_patch(&surface, &piece, Some(u), Some(v));
            }
        }
        for &v in &v_ranges {
            let piece = surface.extract_v(v.0, v.1).unwrap();
            assert_patch(&surface, &piece, None, Some(v));
        }
    }
}

#[test]
fn test_extract_whole_clamped_surface_is_the_surface() {
    let surface = insertion_surface(true, false, false);
    let whole = surface
        .extract(surface.u_param_range(), surface.v_param_range())
        .unwrap();
    assert_eq!(whole, surface);
    // The `u` double knot at 2 of the cubic takes one more copy.
    let (low, high) = surface.split_u_at(2.0).unwrap();
    assert_eq!(
        low.u_pole_count() + high.u_pole_count(),
        surface.u_pole_count() + 2
    );
}

#[test]
fn test_extract_matches_curve_extraction() {
    let surface = insertion_surface(true, true, false);
    let (u, v) = ((4.5, 7.25), (0.25, 1.5));
    let patch = surface.extract(u, v).unwrap();
    // Row `i` of a `u` extraction is the `u` extraction of column `i`.
    let by_u = surface.extract_u(u.0, u.1).unwrap();
    let n_v = surface.v_pole_count();
    let (net, weight_net) = (surface.poles(), surface.weights().unwrap());
    let by_u_net = rows(by_u.poles(), n_v);
    for j in 0..n_v {
        let column: Vec<Point> = net.iter().skip(j).step_by(n_v).copied().collect();
        let weights = weight_net.iter().skip(j).step_by(n_v).copied().collect();
        let curve = BSplineCurve::new_rational(
            3,
            column,
            weights,
            surface.u_knots().to_vec(),
            surface.u_multiplicities().to_vec(),
            true,
        )
        .unwrap()
        .extract(u.0, u.1)
        .unwrap();
        assert_eq!(curve.knots(), by_u.u_knots());
        for (row, q) in by_u_net.iter().zip(curve.poles()) {
            assert!(row[j].distance(*q) < 1e-12, "{:?} vs {q:?}", row[j]);
        }
    }
    assert_patch(&surface, &patch, Some(u), Some(v));
}

#[test]
fn test_split_u_and_v_meet_along_the_split() {
    for surface in insertion_surfaces() {
        let (u0, u1) = surface.u_param_range();
        let (v0, v1) = surface.v_param_range();
        for u in [u0 + 0.37 * (u1 - u0), surface.u_knots()[1]] {
            let (low, high) = surface.split_u_at(u).unwrap();
            assert_patch(&surface, &low, Some((u0, u)), None);
            assert_patch(&surface, &high, Some((u, u1)), None);
            let (low_rows, high_rows) = (
                rows(low.poles(), low.v_pole_count()),
                rows(high.poles(), high.v_pole_count()),
            );
            assert_eq!(low_rows.last(), high_rows.first());
        }
        for v in [v0 + 0.37 * (v1 - v0), surface.v_knots()[1]] {
            let (low, high) = surface.split_v_at(v).unwrap();
            assert_patch(&surface, &low, None, Some((v0, v)));
            assert_patch(&surface, &high, None, Some((v, v1)));
            let low_rows = rows(low.poles(), low.v_pole_count());
            for (a, b) in low_rows.iter().zip(rows(high.poles(), high.v_pole_count())) {
                assert_eq!(a.last(), b.first());
            }
        }
    }
}

#[test]
fn test_reversed_u_and_v_mirror_the_parameters() {
    for surface in insertion_surfaces() {
        let (u0, u1) = surface.u_param_range();
        let (v0, v1) = surface.v_param_range();
        let size = size(&surface);
        for along_u in [true, false] {
            let reversed = if along_u {
                surface.reversed_u()
            } else {
                surface.reversed_v()
            };
            assert_eq!(reversed.u_param_range(), (u0, u1));
            assert_eq!(reversed.v_param_range(), (v0, v1));
            assert_eq!(reversed.is_u_periodic(), surface.is_u_periodic());
            assert_eq!(reversed.is_v_periodic(), surface.is_v_periodic());
            for i in 0..=12 {
                for j in 0..=12 {
                    let u = u0 + (u1 - u0) * (i as f64 + 0.3) / 13.0;
                    let v = v0 + (v1 - v0) * (j as f64 + 0.3) / 13.0;
                    let (mu, mv) = if along_u {
                        (u0 + u1 - u, v)
                    } else {
                        (u, v0 + v1 - v)
                    };
                    let (p, [su, sv], [suu, suv, svv]) = surface.eval_d2(mu, mv);
                    let (q, [tu, tv], [tuu, tuv, tvv]) = reversed.eval_d2(u, v);
                    let (su, sv) = if along_u { (-su, sv) } else { (su, -sv) };
                    let pairs = [
                        (p - Point::ORIGIN, q - Point::ORIGIN),
                        (su, tu),
                        (sv, tv),
                        (suu, tuu),
                        (-suv, tuv),
                        (svv, tvv),
                    ];
                    for (k, (a, b)) in pairs.into_iter().enumerate() {
                        let tol = 1e-12 * size.max(a.magnitude());
                        assert!(
                            (a - b).magnitude() <= tol,
                            "({u}, {v}), output {k}: {a:?} vs {b:?}"
                        );
                    }
                }
            }
            let back = if along_u {
                reversed.reversed_u()
            } else {
                reversed.reversed_v()
            };
            assert_eq!(back.poles(), surface.poles());
            assert_eq!(back.weights(), surface.weights());
            assert_eq!(back.u_multiplicities(), surface.u_multiplicities());
            assert_eq!(back.v_multiplicities(), surface.v_multiplicities());
            for (knots, original, (first, last)) in [
                (back.u_knots(), surface.u_knots(), (u0, u1)),
                (back.v_knots(), surface.v_knots(), (v0, v1)),
            ] {
                let scale = first.abs().max(last.abs());
                for (k, l) in knots.iter().zip(original) {
                    assert!((k - l).abs() <= f64::EPSILON * scale, "{k} vs {l}");
                }
            }
        }
    }
}

#[test]
fn test_surface_split_and_extract_errors() {
    let surface = insertion_surface(true, false, true);
    let (u0, u1) = surface.u_param_range();
    let (v0, v1) = surface.v_param_range();
    assert_eq!(
        surface.split_u_at(f64::NAN),
        Err(BSplineSplitError::NonFinite)
    );
    assert_eq!(surface.split_u_at(u0), Err(BSplineSplitError::EmptyRange));
    assert_eq!(surface.split_v_at(v1), Err(BSplineSplitError::EmptyRange));
    assert_eq!(
        surface.split_v_at(v1 + 1.0),
        Err(BSplineSplitError::OutOfRange)
    );
    assert_eq!(
        surface.extract((u0 - 0.1, u1), (v0, v1)),
        Err(BSplineSplitError::OutOfRange)
    );
    assert_eq!(
        surface.extract((u0, u1), (v0, v0 + 2.0 * (v1 - v0))),
        Err(BSplineSplitError::ExceedsPeriod)
    );
    assert_eq!(
        surface.extract_v(v1, v0),
        Err(BSplineSplitError::EmptyRange)
    );
    assert_eq!(
        surface.extract_u(u0, f64::INFINITY),
        Err(BSplineSplitError::NonFinite)
    );
    // The periodic `v` direction takes a range across its seam.
    assert!(surface.extract_v(v1 - 0.5, v1 + 0.5).is_ok());
}

#[test]
fn test_extract_keeps_untouched_rational_poles_exactly() {
    // Periodic both ways: a cubic `u` with knots 0..6 cut at 0.5 and
    // 5.5 keeps three middle rows; a quadratic `v` with knots 0, 1,
    // 2.5, 3, 4 cut at 0.5 and 3.5 keeps two middle columns.
    let surface = insertion_surface(true, true, true);
    // Each pole with its weight, row by row.
    let net = |s: &BSplineSurface| {
        let poles = s
            .poles()
            .iter()
            .copied()
            .zip(s.weights().unwrap().iter().copied());
        rows(&poles.collect::<Vec<_>>(), s.v_pole_count())
    };
    let column = |net: &[Vec<(Point, f64)>], j: usize| -> Vec<(Point, f64)> {
        net.iter().map(|row| row[j]).collect()
    };
    let original = net(&surface);
    let piece = surface.extract_u(0.5, 5.5).unwrap();
    assert_patch(&surface, &piece, Some((0.5, 5.5)), None);
    let piece_net = net(&piece);
    assert_eq!(piece_net.len(), 9);
    for row in &piece_net[3..6] {
        assert!(original.contains(row), "{row:?} not kept exactly");
    }
    let piece = surface.extract_v(0.5, 3.5).unwrap();
    assert_patch(&surface, &piece, None, Some((0.5, 3.5)));
    let piece_net = net(&piece);
    assert_eq!(piece_net[0].len(), 6);
    for j in 2..4 {
        let kept = (0..original[0].len()).any(|k| column(&original, k) == column(&piece_net, j));
        assert!(kept, "column {j} not kept exactly");
    }
}

// ---- one-sided evaluation ----

fn creased_u_patch() -> BSplineSurface {
    // Degree 2 in `u` with a double interior knot at 0.5 (C0 across
    // it), degree 1 in `v`.
    let row = |x: f64, z: f64| vec![Point::new(x, 0.0, z), Point::new(x, 2.0, z)];
    BSplineSurface::new(
        2,
        1,
        vec![
            row(0.0, 0.0),
            row(1.0, 0.0),
            row(2.0, 1.0),
            row(3.0, 1.0),
            row(4.0, 0.0),
        ],
        vec![0.0, 0.5, 1.0],
        vec![3, 2, 3],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_sided_derivatives_match_split_patches() {
    let surface = creased_u_patch();
    let (low, high) = surface.split_u_at(0.5).unwrap();
    let v = 0.7;
    let below = surface.eval_derivative_sided(0.5, v, 1, 0, Side::Below, Side::Above);
    let above = surface.eval_derivative_sided(0.5, v, 1, 0, Side::Above, Side::Above);
    // C0 in `u`: positions agree, first derivatives do not.
    let (pb, _) = surface.eval_d1_sided(0.5, v, Side::Below, Side::Above);
    let (pa, _) = surface.eval_d1_sided(0.5, v, Side::Above, Side::Above);
    assert!(pb.distance(pa) < 1e-12);
    assert!((below - above).magnitude() > 1e-6);
    // Each side matches its split patch evaluated at the cut.
    let (_, low_hi) = low.u_bounds();
    let (high_lo, _) = high.u_bounds();
    assert!((low.eval_derivative(low_hi, v, 1, 0) - below).magnitude() < 1e-9);
    assert!((high.eval_derivative(high_lo, v, 1, 0) - above).magnitude() < 1e-9);
    let below2 = surface.eval_derivative_sided(0.5, v, 2, 0, Side::Below, Side::Above);
    let above2 = surface.eval_derivative_sided(0.5, v, 2, 0, Side::Above, Side::Above);
    assert!((low.eval_derivative(low_hi, v, 2, 0) - below2).magnitude() < 1e-9);
    assert!((high.eval_derivative(high_lo, v, 2, 0) - above2).magnitude() < 1e-9);
}

#[test]
fn test_sided_evaluation_agrees_off_knots() {
    let surface = creased_u_patch();
    for (u, v) in [(0.1, 0.2), (0.49, 0.7), (0.51, 0.3), (0.9, 0.9)] {
        let (p, [su, sv]) = surface.eval_d1(u, v);
        assert_eq!(
            surface.eval_d1_sided(u, v, Side::Below, Side::Below),
            (p, [su, sv])
        );
        assert_eq!(
            surface.eval_d1_sided(u, v, Side::Above, Side::Above),
            (p, [su, sv])
        );
    }
}

#[test]
fn test_sided_derivatives_at_creased_periodic_seam() {
    // Periodic quadratic in `u`, C0 at its seam (multiplicity 2 at
    // both ends), linear in `v`.
    let rows = (0..6)
        .map(|i| {
            let angle = i as f64;
            let (x, y) = (2.0 * angle.cos(), 2.0 * angle.sin());
            let z = (i % 2) as f64;
            vec![Point::new(x, y, z), Point::new(x, y, z + 3.0)]
        })
        .collect();
    let surface = BSplineSurface::new(
        2,
        1,
        rows,
        vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        vec![2, 1, 1, 1, 1, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        true,
        false,
    )
    .unwrap();
    let (v, step) = (0.4, 1e-9);
    for seam in [0.0, 5.0] {
        let below = surface.eval_derivative_sided(seam, v, 1, 0, Side::Below, Side::Above);
        let above = surface.eval_derivative_sided(seam, v, 1, 0, Side::Above, Side::Above);
        assert!((below - above).magnitude() > 0.1);
        assert!((surface.eval_derivative(5.0 - step, v, 1, 0) - below).magnitude() < 1e-6);
        assert!((surface.eval_derivative(step, v, 1, 0) - above).magnitude() < 1e-6);
    }
    // Plain evaluation at the lower bound takes the above side.
    assert_eq!(
        surface.eval_d2(0.0, v),
        surface.eval_d2_sided(0.0, v, Side::Above, Side::Above)
    );
}

#[test]
fn test_sided_evaluation_on_periodic_matches_plain() {
    let surface = periodic_tube();
    for (u, v) in [(0.0, 0.5), (3.0, 0.2)] {
        let (p, [su, sv]) = surface.eval_d1(u, v);
        assert_eq!(
            surface.eval_d1_sided(u, v, Side::Below, Side::Above),
            (p, [su, sv])
        );
    }
}

#[test]
fn test_sided_matches_neighbor_spans_and_point_is_exact() {
    // Bicubic patch with a double u-knot at 0.6: sided derivatives
    // match epsilon-step neighbors, positions are bitwise identical.
    let net = |x: f64, z: f64| {
        vec![
            Point::new(x, 0.0, z),
            Point::new(x, 1.0, z + 1.0),
            Point::new(x, 2.0, z),
        ]
    };
    let surface = BSplineSurface::new(
        3,
        2,
        vec![
            net(0.0, 0.0),
            net(1.0, 1.0),
            net(2.0, 0.0),
            net(3.0, 2.0),
            net(4.0, 1.0),
            net(5.0, 0.0),
        ],
        vec![0.0, 0.6, 1.0],
        vec![4, 2, 4],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
        false,
    )
    .unwrap();
    let v = 0.35;
    for du in [1u32, 2] {
        let below = surface.eval_derivative_sided(0.6, v, du, 0, Side::Below, Side::Above);
        let above = surface.eval_derivative_sided(0.6, v, du, 0, Side::Above, Side::Above);
        let left = surface.eval_derivative(0.6 - 1e-13, v, du, 0);
        let right = surface.eval_derivative(0.6 + 1e-13, v, du, 0);
        assert!((below - left).magnitude() < 1e-8, "du = {du}");
        assert!((above - right).magnitude() < 1e-8, "du = {du}");
        assert_eq!(above, surface.eval_derivative(0.6, v, du, 0));
    }
    // The double knot leaves the patch C1 in u: only the second
    // derivative tells the sides apart, so it must differ.
    let jump = surface.eval_derivative_sided(0.6, v, 2, 0, Side::Below, Side::Above)
        - surface.eval_derivative_sided(0.6, v, 2, 0, Side::Above, Side::Above);
    assert!(jump.magnitude() > 1.0, "{jump:?}");
    let (pb, _) = surface.eval_d1_sided(0.6, v, Side::Below, Side::Above);
    let (pa, _) = surface.eval_d1_sided(0.6, v, Side::Above, Side::Above);
    assert_eq!(pb, pa);
}
