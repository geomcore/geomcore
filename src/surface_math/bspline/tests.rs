use super::*;

// ---- fixed-degree point path ----

#[test]
fn test_point_of_low_degrees_matches_general_path() {
    // Every pair of degrees 1 to 3 and a degree-4 direction, which takes
    // the general path, over uneven knots, polynomial and rational.
    for (p, q) in [(1, 3), (2, 2), (3, 1), (3, 3), (2, 4)] {
        let direction = |degree: usize| {
            KnotVector::new(
                degree,
                vec![0.0, 0.4, 0.55, 1.0],
                vec![degree as u32 + 1, 1, 1, degree as u32 + 1],
                false,
                degree + 3,
            )
            .unwrap()
        };
        let (u_knots, v_knots) = (direction(p), direction(q));
        let count = (p + 3) * (q + 3);
        let poles: Vec<Point> = (0..count)
            .map(|k| {
                let (i, j) = ((k / (q + 3)) as f64, (k % (q + 3)) as f64);
                Point::new(i, j, (0.7 * i + 1.3 * j).sin())
            })
            .collect();
        let weights: Vec<f64> = (0..count).map(|k| 0.6 + 0.2 * (k % 4) as f64).collect();
        for weights in [None, Some(weights.as_slice())] {
            let view = Poles {
                points: &poles,
                weights,
            };
            for a in -1..=11 {
                for b in -1..=11 {
                    let (u, v) = (a as f64 / 10.0, b as f64 / 10.0);
                    let point = surface_point(&u_knots, &v_knots, view, u, v);
                    let [s, ..] = surface_derivatives::<0>(&u_knots, &v_knots, view, u, v);
                    let close = (point.x - s.x)
                        .abs()
                        .max((point.y - s.y).abs())
                        .max((point.z - s.z).abs());
                    assert!(
                        close < 1e-12,
                        "degrees {p} {q} at ({u}, {v}): {point:?} vs {s:?}"
                    );
                }
            }
        }
    }
}

// ---- bilinear patch by hand: the smallest tensor product ----

// Degree-1 x degree-1 patch over the unit square with corner poles
//   P00=(0,0,0) P01=(0,1,0)
//   P10=(1,0,0) P11=(1,1,1)
// (row-major [i_u][i_v]). Bilinear interpolation gives
//   S(u,v) = (u, v, u*v), Su = (1, 0, v), Sv = (0, 1, u).
fn bilinear_poles() -> Vec<Point> {
    vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 1.0),
    ]
}

fn unit_linear() -> KnotVector {
    KnotVector::new(1, vec![0.0, 1.0], vec![2, 2], false, 2).unwrap()
}

#[test]
fn test_bilinear_patch_point_by_hand() {
    let poles = bilinear_poles();
    let knots = unit_linear();
    let view = Poles {
        points: &poles,
        weights: None,
    };
    let s = surface_point(&knots, &knots, view, 0.5, 0.5);
    assert_close3(&[s.x, s.y, s.z], &[0.5, 0.5, 0.25]);
}

#[test]
fn test_bilinear_patch_derivatives_by_hand() {
    let poles = bilinear_poles();
    let knots = unit_linear();
    let view = Poles {
        points: &poles,
        weights: None,
    };
    // At (u, v) = (0.3, 0.7): S=(0.3,0.7,0.21), Su=(1,0,0.7), Sv=(0,1,0.3).
    let [s, su, sv, ..] = surface_derivatives::<1>(&knots, &knots, view, 0.3, 0.7);
    assert_close3(&[s.x, s.y, s.z], &[0.3, 0.7, 0.21]);
    assert_close3(&[su.x, su.y, su.z], &[1.0, 0.0, 0.7]);
    assert_close3(&[sv.x, sv.y, sv.z], &[0.0, 1.0, 0.3]);
}

// ---- fixture replay: the correctness gate ----

#[test]
fn test_fixture_replay_all_cases() {
    replay_fixture();
}

// ---- helpers ----

fn assert_close3(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (i, (&a, &e)) in actual.iter().zip(expected.iter()).enumerate() {
        let tol = 1e-9 * f64::max(1.0, e.abs());
        assert!((a - e).abs() <= tol, "index {i}: {a} vs {e}");
    }
}

use serde_json::Value;

fn replay_fixture() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/surfaces_bspline.json"
    );
    let text = std::fs::read_to_string(path).expect("read fixture");
    let json: Value = serde_json::from_str(&text).expect("parse fixture");
    for case in json["cases"].as_array().expect("cases array") {
        replay_case(case);
    }
}

fn arr_f64(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn arr_u32(v: &Value) -> Vec<u32> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_u64().unwrap() as u32)
        .collect()
}

fn replay_case(case: &Value) {
    let name = case["name"].as_str().unwrap();
    let u_degree = case["u_degree"].as_u64().unwrap() as usize;
    let v_degree = case["v_degree"].as_u64().unwrap() as usize;
    let u_periodic = case["u_periodic"].as_bool().unwrap();
    let v_periodic = case["v_periodic"].as_bool().unwrap();
    let u_knots = arr_f64(&case["u_knots"]);
    let v_knots = arr_f64(&case["v_knots"]);
    let u_mults = arr_u32(&case["u_mults"]);
    let v_mults = arr_u32(&case["v_mults"]);

    // poles[i_u][i_v] = [x, y, z]
    let pole_grid: Vec<Vec<Vec<f64>>> = case["poles"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| row.as_array().unwrap().iter().map(arr_f64).collect())
        .collect();
    let n_u = pole_grid.len();
    let n_v = pole_grid[0].len();

    let weights: Option<Vec<Vec<f64>>> = case.get("weights").and_then(|w| {
        if w.is_null() {
            None
        } else {
            Some(w.as_array().unwrap().iter().map(arr_f64).collect())
        }
    });
    let poles: Vec<Point> = pole_grid
        .iter()
        .flatten()
        .map(|p| Point::new(p[0], p[1], p[2]))
        .collect();
    let weights: Option<Vec<f64>> = weights.map(|w| w.concat());
    let view = Poles {
        points: &poles,
        weights: weights.as_deref(),
    };
    let u_vector = KnotVector::new(u_degree, u_knots, u_mults, u_periodic, n_u).unwrap();
    let v_vector = KnotVector::new(v_degree, v_knots, v_mults, v_periodic, n_v).unwrap();

    for sample in case["samples"].as_array().unwrap() {
        let u = sample["u"].as_f64().unwrap();
        let v = sample["v"].as_f64().unwrap();

        let s0 = surface_point(&u_vector, &v_vector, view, u, v);
        let [s, su, sv, ..] = surface_derivatives::<1>(&u_vector, &v_vector, view, u, v);
        let [_, _, _, suu, suv, svv] = surface_derivatives::<2>(&u_vector, &v_vector, view, u, v);
        for (key, d) in [("d2u", suu), ("d2uv", suv), ("d2v", svv)] {
            check(name, u, v, key, &[d.x, d.y, d.z], &arr_f64(&sample[key]));
        }
        // d0 and d1 agree on the point.
        check(name, u, v, "d0", &[s0.x, s0.y, s0.z], &[s.x, s.y, s.z]);
        check(
            name,
            u,
            v,
            "point",
            &[s.x, s.y, s.z],
            &arr_f64(&sample["point"]),
        );
        check(
            name,
            u,
            v,
            "d1u",
            &[su.x, su.y, su.z],
            &arr_f64(&sample["d1u"]),
        );
        check(
            name,
            u,
            v,
            "d1v",
            &[sv.x, sv.y, sv.z],
            &arr_f64(&sample["d1v"]),
        );
    }
}

fn check(name: &str, u: f64, v: f64, key: &str, got: &[f64], expected: &[f64]) {
    for (c, (&g, &e)) in got.iter().zip(expected.iter()).enumerate() {
        let tol = 1e-7 * f64::max(1.0, e.abs());
        assert!(
            (g - e).abs() <= tol,
            "{name} (u={u}, v={v}) {key}[{c}]: got {g}, expected {e} (tol {tol})"
        );
    }
}
