//! `BezierSurface` through the public API: construction, evaluation
//! against hand-computed values, editing, and the round trip with
//! `BSplineSurface`.

use geomcore::surfaces::ParametricSurface;
use geomcore::{
    BSplineSurface, BezierConstructionError, BezierSurface, BoundingBox, DegreeElevationError,
    HasBoundingBox, Point, Side, Tolerance, Vector,
};

/// The grid `(i, j, i * j)`: the surface `(2u, 2v, 4uv)` for degrees 2 x 2
/// (the Bernstein form reproduces a linear function of the index).
fn saddle_poles() -> Vec<Vec<Point>> {
    (0..3)
        .map(|i| {
            (0..3)
                .map(|j| Point::new(i as f64, j as f64, (i * j) as f64))
                .collect()
        })
        .collect()
}

fn saddle() -> BezierSurface {
    BezierSurface::new(2, 2, saddle_poles()).unwrap()
}

/// A rational quarter cylinder of radius 1 and height 1 about `z`: a
/// quadratic arc in `u`, a line in `v`.
fn quarter_cylinder() -> BezierSurface {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let arc = |z: f64| {
        vec![
            Point::new(1.0, 0.0, z),
            Point::new(1.0, 1.0, z),
            Point::new(0.0, 1.0, z),
        ]
    };
    // poles[i_u][i_v]: three `u` rows of two `v` poles.
    let poles = (0..3).map(|i| vec![arc(0.0)[i], arc(1.0)[i]]).collect();
    let weights = [1.0, w, 1.0].iter().map(|&w| vec![w, w]).collect();
    BezierSurface::new_rational(2, 1, poles, weights).unwrap()
}

fn assert_close(a: Point, b: Point, tol: f64) {
    assert!(a.distance(b) <= tol, "{a:?} vs {b:?}");
}

fn assert_same_shape(a: &BezierSurface, b: &BezierSurface) {
    for i in 0..=6 {
        for j in 0..=6 {
            let u = a.u_bounds().0 + (a.u_bounds().1 - a.u_bounds().0) * i as f64 / 6.0;
            let v = a.v_bounds().0 + (a.v_bounds().1 - a.v_bounds().0) * j as f64 / 6.0;
            assert_close(a.eval_point(u, v), b.eval_point(u, v), 1e-12);
        }
    }
}

// ---- construction / accessors ----

#[test]
fn test_new_accessors() {
    let surface = saddle();
    assert_eq!((surface.u_degree(), surface.v_degree()), (2, 2));
    assert_eq!((surface.u_pole_count(), surface.v_pole_count()), (3, 3));
    assert_eq!(surface.poles().len(), 9);
    // Row-major: a row per `u` pole, as `BSplineSurface`.
    assert_eq!(surface.poles()[1], Point::new(0.0, 1.0, 0.0));
    assert_eq!(surface.poles()[3], Point::new(1.0, 0.0, 0.0));
    assert_eq!(surface.weights(), None);
    assert!(!surface.is_rational());
    assert_eq!(surface.u_bounds(), (0.0, 1.0));
    assert_eq!(surface.v_bounds(), (0.0, 1.0));
}

#[test]
fn test_new_mixed_degrees() {
    let poles = (0..2)
        .map(|i| {
            (0..4)
                .map(|j| Point::new(i as f64, j as f64, 0.0))
                .collect()
        })
        .collect();
    let surface = BezierSurface::new(1, 3, poles).unwrap();
    assert_eq!((surface.u_degree(), surface.v_degree()), (1, 3));
    assert_eq!((surface.u_pole_count(), surface.v_pole_count()), (2, 4));
}

#[test]
fn test_new_rational_accessors() {
    let surface = quarter_cylinder();
    let w = std::f64::consts::FRAC_1_SQRT_2;
    assert!(surface.is_rational());
    assert_eq!(
        surface.weights(),
        Some([1.0, 1.0, w, w, 1.0, 1.0].as_slice())
    );
}

#[test]
fn test_new_errors() {
    let grid = |rows: usize, cols: usize| vec![vec![Point::ORIGIN; cols]; rows];
    let err = |u, v, poles| BezierSurface::new(u, v, poles).unwrap_err();
    assert_eq!(
        err(0, 1, grid(1, 2)),
        BezierConstructionError::InvalidDegree
    );
    assert_eq!(
        err(1, 26, grid(2, 27)),
        BezierConstructionError::InvalidDegree
    );
    assert!(BezierSurface::new(25, 1, grid(26, 2)).is_ok());
    assert_eq!(
        err(2, 2, grid(2, 3)),
        BezierConstructionError::PoleCountMismatch
    );
    assert_eq!(
        err(2, 2, grid(3, 4)),
        BezierConstructionError::PoleCountMismatch
    );
    assert_eq!(
        err(1, 1, Vec::new()),
        BezierConstructionError::PoleCountMismatch
    );
    let mut ragged = grid(2, 2);
    ragged[1].pop();
    assert_eq!(
        err(1, 1, ragged),
        BezierConstructionError::PoleCountMismatch
    );
}

#[test]
fn test_new_rational_errors() {
    let poles = vec![vec![Point::ORIGIN; 2]; 2];
    let build = |weights| BezierSurface::new_rational(1, 1, poles.clone(), weights);
    assert_eq!(
        build(vec![vec![1.0; 2]]).unwrap_err(),
        BezierConstructionError::WeightCountMismatch
    );
    assert_eq!(
        build(vec![vec![1.0; 2], vec![1.0; 3]]).unwrap_err(),
        BezierConstructionError::WeightCountMismatch
    );
    assert_eq!(
        build(vec![vec![1.0, 1.0], vec![1.0, 0.0]]).unwrap_err(),
        BezierConstructionError::NonPositiveWeight
    );
}

// ---- evaluation against known values ----

#[test]
fn test_eval_point_by_hand() {
    let surface = saddle();
    assert_eq!(surface.eval_point(0.0, 0.0), Point::new(0.0, 0.0, 0.0));
    assert_eq!(surface.eval_point(1.0, 1.0), Point::new(2.0, 2.0, 4.0));
    assert_eq!(surface.eval_point(1.0, 0.0), Point::new(2.0, 0.0, 0.0));
    // (2u, 2v, 4uv)
    assert_close(
        surface.eval_point(0.25, 0.5),
        Point::new(0.5, 1.0, 0.5),
        1e-14,
    );
    assert_eq!(
        surface.eval_points(&[(0.25, 0.5), (1.0, 1.0)]),
        vec![surface.eval_point(0.25, 0.5), Point::new(2.0, 2.0, 4.0)]
    );
}

#[test]
fn test_eval_derivatives_by_hand() {
    let surface = saddle();
    let (u, v) = (0.25, 0.5);
    let (point, [su, sv], [suu, suv, svv]) = surface.eval_d2(u, v);
    assert_close(point, Point::new(0.5, 1.0, 0.5), 1e-14);
    // Su = (2, 0, 4v), Sv = (0, 2, 4u), Suv = (0, 0, 4).
    assert!((su - Vector::new(2.0, 0.0, 2.0)).magnitude() < 1e-12);
    assert!((sv - Vector::new(0.0, 2.0, 1.0)).magnitude() < 1e-12);
    assert!(suu.magnitude() < 1e-12);
    assert!((suv - Vector::new(0.0, 0.0, 4.0)).magnitude() < 1e-12);
    assert!(svv.magnitude() < 1e-12);
    assert_eq!(surface.eval_d1(u, v).1, [su, sv]);
    assert_eq!(surface.eval_derivative(u, v, 1, 1), suv);
}

#[test]
fn test_sided_evaluation_agrees_with_plain() {
    let surface = saddle();
    let (u, v) = (0.3, 0.7);
    for (us, vs) in [
        (Side::Below, Side::Below),
        (Side::Above, Side::Above),
        (Side::Below, Side::Above),
    ] {
        assert_eq!(surface.eval_d1_sided(u, v, us, vs), surface.eval_d1(u, v));
        assert_eq!(surface.eval_d2_sided(u, v, us, vs), surface.eval_d2(u, v));
        assert_eq!(
            surface.eval_derivative_sided(u, v, 1, 1, us, vs),
            surface.eval_derivative(u, v, 1, 1)
        );
    }
}

#[test]
fn test_rational_quarter_cylinder_stays_on_the_cylinder() {
    let surface = quarter_cylinder();
    for i in 0..=10 {
        for j in 0..=4 {
            let (u, v) = (i as f64 / 10.0, j as f64 / 4.0);
            let p = surface.eval_point(u, v);
            assert!((p.x.hypot(p.y) - 1.0).abs() < 1e-12, "{p:?}");
            assert!((p.z - v).abs() < 1e-12, "{p:?}");
        }
    }
    let end = surface.eval_point(1.0, 0.5);
    assert_close(end, Point::new(0.0, 1.0, 0.5), 1e-12);
}

#[test]
fn test_parametric_surface_trait_delegates() {
    let surface = saddle();
    assert_eq!(ParametricSurface::u_bounds(&surface), (0.0, 1.0));
    assert_eq!(ParametricSurface::v_bounds(&surface), (0.0, 1.0));
    assert_eq!(ParametricSurface::u_period(&surface), None);
    assert_eq!(ParametricSurface::v_period(&surface), None);
    assert!(!ParametricSurface::is_u_closed(
        &surface,
        Tolerance::DEFAULT
    ));
    assert_eq!(
        ParametricSurface::eval_point(&surface, 0.25, 0.5),
        surface.eval_point(0.25, 0.5)
    );
    assert_eq!(
        ParametricSurface::eval_d2(&surface, 0.25, 0.5),
        surface.eval_d2(0.25, 0.5)
    );
    let normal = ParametricSurface::normal(&surface, 0.0, 0.0, Tolerance::DEFAULT).unwrap();
    assert!((normal - Vector::Z).magnitude() < 1e-12);
}

#[test]
fn test_project_point_and_contains() {
    let surface = saddle();
    let on = surface.eval_point(0.4, 0.6);
    assert!(surface.contains(on, Tolerance::DEFAULT));
    let projection = surface.project_point(on, Tolerance::DEFAULT);
    assert_eq!(projection.distance, 0.0);
    assert!((projection.u - 0.4).abs() < 1e-7 && (projection.v - 0.6).abs() < 1e-7);
    assert!(!surface.contains(Point::new(0.0, 0.0, 5.0), Tolerance::DEFAULT));
}

// ---- B-spline round trip ----

#[test]
fn test_to_bspline_is_one_clamped_span_each_way() {
    let spline = saddle().to_bspline();
    assert_eq!(spline.u_knots(), &[0.0, 1.0]);
    assert_eq!(spline.v_knots(), &[0.0, 1.0]);
    assert_eq!(spline.u_multiplicities(), &[3, 3]);
    assert_eq!(spline.v_multiplicities(), &[3, 3]);
    assert!(!spline.is_u_periodic() && !spline.is_v_periodic());
    assert_eq!(spline.poles(), saddle().poles());
    assert_close(
        spline.eval_point(0.25, 0.5),
        saddle().eval_point(0.25, 0.5),
        1e-15,
    );
}

#[test]
fn test_bspline_round_trip() {
    for surface in [saddle(), quarter_cylinder()] {
        let round = BezierSurface::from_bspline(&surface.to_bspline()).unwrap();
        assert_eq!(round, surface);
    }
}

#[test]
fn test_from_bspline_rejects_what_is_not_one_span() {
    let spline = saddle().to_bspline();
    let refined_u = spline.insert_u_knot(0.5, 1).unwrap();
    let refined_v = spline.insert_v_knot(0.5, 1).unwrap();
    assert_eq!(
        BezierSurface::from_bspline(&refined_u),
        Err(BezierConstructionError::NotABezier)
    );
    assert_eq!(
        BezierSurface::from_bspline(&refined_v),
        Err(BezierConstructionError::NotABezier)
    );
    // Unclamped: one span but end multiplicities below degree + 1.
    let open = BSplineSurface::new(
        1,
        1,
        vec![vec![Point::ORIGIN; 2]; 2],
        vec![0.0, 1.0, 2.0, 3.0],
        vec![1, 1, 1, 1],
        vec![0.0, 1.0, 2.0, 3.0],
        vec![1, 1, 1, 1],
        false,
        false,
    )
    .unwrap();
    assert_eq!(
        BezierSurface::from_bspline(&open),
        Err(BezierConstructionError::NotABezier)
    );
}

#[test]
fn test_from_bspline_keeps_the_span_as_bounds() {
    let spline = BSplineSurface::new(
        1,
        1,
        vec![
            vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
            vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
        ],
        vec![2.0, 5.0],
        vec![2, 2],
        vec![-1.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap();
    let surface = BezierSurface::from_bspline(&spline).unwrap();
    assert_eq!(surface.u_bounds(), (2.0, 5.0));
    assert_eq!(surface.v_bounds(), (-1.0, 1.0));
    assert_eq!(surface.to_bspline(), spline);
}

#[test]
fn test_to_bezier_patches_returns_bezier_surfaces() {
    let poles: Vec<Vec<Point>> = (0..4)
        .map(|i| {
            (0..3)
                .map(|j| Point::new(i as f64, j as f64, (i * j) as f64))
                .collect()
        })
        .collect();
    let spline = BSplineSurface::new(
        2,
        2,
        poles,
        vec![0.0, 1.0, 2.0],
        vec![3, 1, 3],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
        false,
    )
    .unwrap();
    let patches: Vec<Vec<BezierSurface>> = spline.to_bezier_patches();
    assert_eq!((patches.len(), patches[0].len()), (2, 1));
    let patch = &patches[1][0];
    assert_eq!(patch.u_bounds(), (1.0, 2.0));
    assert_eq!(patch.v_bounds(), (0.0, 1.0));
    assert_eq!((patch.u_pole_count(), patch.v_pole_count()), (3, 3));
    for (u, v) in [(1.0, 0.0), (1.5, 0.5), (2.0, 1.0), (1.2, 0.9)] {
        assert_close(patch.eval_point(u, v), spline.eval_point(u, v), 1e-12);
    }
    assert_eq!(
        BezierSurface::from_bspline(&patch.to_bspline()).as_ref(),
        Ok(patch)
    );
}

// ---- editing ----

#[test]
fn test_split_u_at() {
    let surface = saddle();
    let (low, high) = surface.split_u_at(0.25).unwrap();
    assert_eq!(
        (low.u_bounds(), high.u_bounds()),
        ((0.0, 0.25), (0.25, 1.0))
    );
    assert_eq!(low.v_bounds(), (0.0, 1.0));
    assert_eq!((low.u_degree(), low.u_pole_count()), (2, 3));
    // The last row of `low` is the first of `high`: the iso-curve at `u`.
    assert_eq!(&low.poles()[6..9], &high.poles()[0..3]);
    for (u, v) in [(0.1, 0.3), (0.2, 0.9)] {
        assert_close(low.eval_point(u, v), surface.eval_point(u, v), 1e-12);
    }
    for (u, v) in [(0.5, 0.3), (1.0, 1.0)] {
        assert_close(high.eval_point(u, v), surface.eval_point(u, v), 1e-12);
    }
}

#[test]
fn test_split_v_at() {
    let surface = quarter_cylinder();
    let (low, high) = surface.split_v_at(0.4).unwrap();
    assert_eq!((low.v_bounds(), high.v_bounds()), ((0.0, 0.4), (0.4, 1.0)));
    assert_eq!(low.u_bounds(), (0.0, 1.0));
    assert!(low.is_rational() && high.is_rational());
    for (u, v) in [(0.3, 0.2), (0.9, 0.4)] {
        assert_close(low.eval_point(u, v), surface.eval_point(u, v), 1e-12);
    }
    for (u, v) in [(0.3, 0.7), (0.9, 1.0)] {
        assert_close(high.eval_point(u, v), surface.eval_point(u, v), 1e-12);
    }
}

#[test]
fn test_split_errors() {
    use geomcore::BSplineSplitError;
    let surface = saddle();
    assert_eq!(surface.split_u_at(0.0), Err(BSplineSplitError::EmptyRange));
    assert_eq!(surface.split_v_at(1.0), Err(BSplineSplitError::EmptyRange));
    assert_eq!(surface.split_u_at(2.0), Err(BSplineSplitError::OutOfRange));
    assert_eq!(
        surface.split_v_at(f64::NAN),
        Err(BSplineSplitError::NonFinite)
    );
}

#[test]
fn test_elevate_degree_keeps_the_shape() {
    for surface in [saddle(), quarter_cylinder()] {
        let up = surface.elevate_u_degree(2).unwrap();
        assert_eq!(up.u_degree(), surface.u_degree() + 2);
        assert_eq!(up.u_pole_count(), surface.u_pole_count() + 2);
        assert_eq!(up.v_degree(), surface.v_degree());
        assert_eq!(up.is_rational(), surface.is_rational());
        assert_same_shape(&surface, &up);

        let vp = surface.elevate_v_degree(1).unwrap();
        assert_eq!(vp.v_degree(), surface.v_degree() + 1);
        assert_eq!(vp.u_degree(), surface.u_degree());
        assert_same_shape(&surface, &vp);
        // One span still: it converts back from a B-spline.
        assert!(BezierSurface::from_bspline(&vp.to_bspline()).is_ok());
    }
    assert_eq!(saddle().elevate_u_degree(0).unwrap(), saddle());
}

#[test]
fn test_elevate_degree_too_high() {
    let surface = BezierSurface::new(25, 1, vec![vec![Point::ORIGIN; 2]; 26]).unwrap();
    assert_eq!(
        surface.elevate_u_degree(1),
        Err(DegreeElevationError::DegreeTooHigh)
    );
}

#[test]
fn test_reversed_mirrors_the_parameters() {
    for surface in [saddle(), quarter_cylinder()] {
        let ru = surface.reversed_u();
        let rv = surface.reversed_v();
        for (u, v) in [(0.0, 0.0), (0.2, 0.7), (1.0, 0.4)] {
            assert_close(ru.eval_point(u, v), surface.eval_point(1.0 - u, v), 1e-12);
            assert_close(rv.eval_point(u, v), surface.eval_point(u, 1.0 - v), 1e-12);
        }
        assert_eq!(ru.reversed_u(), surface);
        assert_eq!(rv.reversed_v(), surface);
    }
}

// ---- bounding box ----

#[test]
fn test_bounding_box_is_the_pole_hull() {
    let b: BoundingBox = saddle().bounding_box().unwrap();
    assert_eq!(b.min(), Point::new(0.0, 0.0, 0.0));
    assert_eq!(b.max(), Point::new(2.0, 2.0, 4.0));
}
