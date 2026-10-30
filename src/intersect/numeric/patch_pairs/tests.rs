use super::*;
use crate::Point;
use crate::curves::BSplineCurve;
use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::intersect::numeric::marching::grid_marching;
use crate::intersect::numeric::marching::marched_curves;
use crate::surfaces::ParametricSurface;

const TOL: Tolerance = Tolerance::DEFAULT;

/// `surface` as its own patch.
fn whole(surface: &BSplineSurface) -> SurfacePatches<'_> {
    SurfacePatches::new(surface, || None).unwrap()
}

/// A bicubic B-spline over `[0, 5]^2`, clamped, with the interior
/// knots `knots` in both directions and its poles above the Greville
/// abscissae, so that `S(u, v) = (u, v, z)` when polynomial: heights
/// `height(x, y)`, weights `weight(i, j)` if given.
fn graph(
    knots: &[f64],
    height: impl Fn(f64, f64) -> f64,
    weight: Option<&dyn Fn(usize, usize) -> f64>,
) -> BSplineSurface {
    let flat: Vec<f64> = [0.0; 4]
        .iter()
        .chain(knots)
        .chain(&[5.0; 4])
        .copied()
        .collect();
    let n = knots.len() + 4;
    let greville: Vec<f64> = (0..n)
        .map(|i| flat[i + 1..i + 4].iter().sum::<f64>() / 3.0)
        .collect();
    let poles = greville
        .iter()
        .map(|&x| {
            greville
                .iter()
                .map(|&y| Point::new(x, y, height(x, y)))
                .collect()
        })
        .collect();
    let breaks: Vec<f64> = [0.0].iter().chain(knots).chain(&[5.0]).copied().collect();
    let mut mults = vec![1; breaks.len()];
    (mults[0], mults[breaks.len() - 1]) = (4, 4);
    match weight {
        None => BSplineSurface::new(
            3,
            3,
            poles,
            breaks.clone(),
            mults.clone(),
            breaks,
            mults,
            false,
            false,
        ),
        Some(w) => BSplineSurface::new_rational(
            3,
            3,
            poles,
            (0..n).map(|i| (0..n).map(|j| w(i, j)).collect()).collect(),
            breaks.clone(),
            mults.clone(),
            breaks,
            mults,
            false,
            false,
        ),
    }
    .unwrap()
}

const THIRDS: [f64; 2] = [5.0 / 3.0, 10.0 / 3.0];

fn wavy() -> BSplineSurface {
    graph(&THIRDS, |x, y| (0.25 * x * y).sin(), None)
}

/// Surfaces crossing [`wavy`] in several curves, open and closed.
fn crossing_wavy() -> Vec<BSplineSurface> {
    vec![
        graph(
            &THIRDS,
            |x, y| 0.5 * (0.8 * x).cos() * (0.7 * y).cos(),
            None,
        ),
        graph(&THIRDS, |_, _| 0.3, None),
        graph(&THIRDS, |x, y| 0.1 * (x - y), None),
        graph(&THIRDS, |x, y| 0.6 * x.sin() * (0.6 * y).cos() - 0.1, None),
    ]
}

fn samples(curve: &BSplineCurve, n: usize) -> impl Iterator<Item = Point> + '_ {
    let (t0, t1) = curve.bounds();
    (0..=n).map(move |k| curve.eval_point(t0 + (t1 - t0) * k as f64 / n as f64))
}

fn is_closed(curve: &BSplineCurve) -> bool {
    let (t0, t1) = curve.bounds();
    curve.eval_point(t0).distance(curve.eval_point(t1)) < 1e-9
}

fn distance_to(curves: &[BSplineCurve], p: Point) -> f64 {
    curves
        .iter()
        .map(|c| c.project_point(p, TOL).distance)
        .fold(f64::INFINITY, f64::min)
}

/// Asserts every curve lies on both surfaces, and that no curve runs
/// along another (each loop is traced once).
#[track_caller]
fn assert_on_both_once(curves: &[BSplineCurve], a: &BSplineSurface, b: &BSplineSurface) {
    for (i, curve) in curves.iter().enumerate() {
        for p in samples(curve, 40) {
            assert!(a.project_point(p, TOL).distance < 1e-3, "{p:?}");
            assert!(b.project_point(p, TOL).distance < 1e-3, "{p:?}");
        }
        for (j, other) in curves.iter().enumerate() {
            let near = samples(curve, 20)
                .filter(|&p| other.project_point(p, TOL).distance < 1e-3)
                .count();
            assert!(i == j || near <= 2, "curves {i} and {j} overlap");
        }
    }
}

/// Asserts that `curves` pass near every sign change of the height
/// difference of `a` and `b` (graphs over the same parameters) along
/// the lines of a fine grid.
#[track_caller]
fn assert_every_crossing(curves: &[BSplineCurve], a: &BSplineSurface, b: &BSplineSurface) {
    let m = 120;
    let at = |k: usize| 5.0 * k as f64 / m as f64;
    let d = |u: f64, v: f64| a.eval_point(u, v).z - b.eval_point(u, v).z;
    let mut crossings = 0;
    for i in 0..m {
        for j in 0..=m {
            for (u0, v0, u1, v1) in [
                (at(i), at(j), at(i + 1), at(j)),
                (at(j), at(i), at(j), at(i + 1)),
            ] {
                let (d0, d1) = (d(u0, v0), d(u1, v1));
                if d0 * d1 < 0.0 {
                    crossings += 1;
                    let f = d0 / (d0 - d1);
                    let p = a.eval_point(u0 + f * (u1 - u0), v0 + f * (v1 - v0));
                    assert!(distance_to(curves, p) < 0.01, "{p:?} missed");
                }
            }
        }
    }
    assert!(crossings > 0);
}

#[test]
fn test_bspline_pairs_trace_every_curve_once() {
    let a = wavy();
    for b in crossing_wavy() {
        let curves = marched_curves(&a, &b, TOL);
        assert!(patch_pairs(&whole(&a), &whole(&b), TOL).complete);
        assert_on_both_once(&curves, &a, &b);
        assert_every_crossing(&curves, &a, &b);
        let grid = grid_marching(&a, &b, TOL);
        assert_eq!(curves.len(), grid.len());
        for curve in &grid {
            assert!(samples(curve, 20).all(|p| distance_to(&curves, p) < 1e-3));
        }
    }
}

#[test]
fn test_small_loop_through_a_nearly_flat_surface() {
    // A bump of width 0.1 about (2.6, 2.6), within one cell of the
    // seed grid (lines spaced 5 / 24), pokes 0.04 through a gently
    // tilted sheet, which also crosses the bump's flat surroundings
    // along x = 2.1.
    let knots = [1.0, 2.0, 2.4, 2.55, 2.575, 2.6, 2.625, 2.65, 2.8, 3.5, 4.25];
    let at_bump = |x: f64, y: f64| (x - 2.6).abs() < 1e-9 && (y - 2.6).abs() < 1e-9;
    let bump = graph(&knots, |x, y| if at_bump(x, y) { 0.2 } else { 0.0 }, None);
    let sheet = graph(&THIRDS, |x, y| 0.1 * (x - 2.1) + 0.005 * y.sin(), None);
    let spacing = 5.0 / 24.0;
    for (a, b) in [(&bump, &sheet), (&sheet, &bump)] {
        let curves = marched_curves(a, b, TOL);
        assert_on_both_once(&curves, a, b);
        assert_eq!(curves.len(), 2, "{curves:?}");
        let loops: Vec<_> = curves.iter().filter(|c| is_closed(c)).collect();
        assert_eq!(loops.len(), 1);
        let points: Vec<Point> = samples(loops[0], 40).collect();
        let diameter = points
            .iter()
            .flat_map(|p| points.iter().map(move |q| p.distance(*q)))
            .fold(0.0, f64::max);
        assert!(diameter < 0.5 * spacing, "{diameter}");
        assert!(
            points
                .iter()
                .all(|p| p.distance(Point::new(2.6, 2.6, p.z)) < 0.05)
        );
        // The grid seeds trace the line but miss the loop.
        let grid = grid_marching(a, b, TOL);
        assert_eq!(grid.len(), 1);
        assert!(!is_closed(&grid[0]));
        assert!(points.iter().all(|&p| distance_to(&grid, p) > 0.4));
    }
}

#[test]
fn test_open_curves_end_on_the_domain_bounds() {
    // A sheet crossing a flat surface along x = 2.5 - 0.5 sin(0.4 y),
    // from the edge y = 0 to the edge y = 5.
    let flat = graph(&THIRDS, |_, _| 0.0, None);
    let sheet = graph(
        &THIRDS,
        |x, y| 0.4 * (x - 2.5) + 0.2 * (0.4 * y).sin(),
        None,
    );
    for curves in [
        marched_curves(&sheet, &flat, TOL),
        grid_marching(&sheet, &flat, TOL),
    ] {
        assert_on_both_once(&curves, &sheet, &flat);
        assert_eq!(curves.len(), 1);
        let (t0, t1) = curves[0].bounds();
        let mut ends = [curves[0].eval_point(t0).y, curves[0].eval_point(t1).y];
        ends.sort_by(f64::total_cmp);
        assert!(
            ends[0].abs() < 1e-7 && (ends[1] - 5.0).abs() < 1e-7,
            "{ends:?}"
        );
    }
}

#[test]
fn test_rational_pairs_match_the_grid_seeds() {
    let weight = |i: usize, j: usize| 0.6 + 0.2 * ((3 * i + 5 * j) % 4) as f64;
    let a = graph(&THIRDS, |x, y| (0.25 * x * y).sin(), Some(&weight));
    for b in crossing_wavy().into_iter().take(2) {
        let curves = marched_curves(&a, &b, TOL);
        assert_on_both_once(&curves, &a, &b);
        let grid = grid_marching(&a, &b, TOL);
        assert!(!grid.is_empty());
        assert_eq!(curves.len(), grid.len());
        for curve in &grid {
            assert!(samples(curve, 20).all(|p| distance_to(&curves, p) < 1e-3));
        }
    }
}

/// A periodic quadratic tube about the z axis, radius about 2, over
/// `z` in `[-1, 3]`.
fn tube() -> BSplineSurface {
    let poles = (0..6)
        .map(|k| {
            let angle = k as f64 * std::f64::consts::TAU / 6.0;
            let (c, s) = (angle.cos(), angle.sin());
            vec![
                Point::new(2.0 * c, 2.0 * s, -1.0),
                Point::new(2.2 * c, 2.2 * s, 3.0),
            ]
        })
        .collect();
    BSplineSurface::new(
        2,
        1,
        poles,
        (0..=6).map(f64::from).collect(),
        vec![1; 7],
        vec![0.0, 1.0],
        vec![2, 2],
        true,
        false,
    )
    .unwrap()
}

#[test]
fn test_periodic_surface_gives_one_loop_across_its_seam() {
    let tube = tube();
    // A tilted sheet over `[-2.5, 2.5]^2`, through the tube.
    let shift = crate::Transform::translation(crate::Vector::new(-2.5, -2.5, 0.0));
    let sheet = graph(&THIRDS, |x, y| 0.2 * x + 0.1 * (y - 2.5), None);
    let sheet = crate::Transformable::transformed(&sheet, shift).unwrap();
    let curves = marched_curves(&tube, &sheet, TOL);
    assert_on_both_once(&curves, &tube, &sheet);
    assert_eq!(curves.len(), 1);
    assert!(is_closed(&curves[0]));
    let grid = grid_marching(&tube, &sheet, TOL);
    assert_eq!(grid.len(), 1);
    assert!(samples(&grid[0], 20).all(|p| distance_to(&curves, p) < 1e-3));
}

#[test]
fn test_overlapping_surfaces_exhaust_the_search() {
    let a = wavy();
    assert!(!patch_pairs(&whole(&a), &whole(&a), TOL).complete);
    // The grid seeds take over; the surfaces are tangent everywhere,
    // so nothing traces.
    assert!(marched_curves(&a, &a, TOL).is_empty());
}
