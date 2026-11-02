use super::*;
use crate::{IntersectExact, Line, LineCircleIntersection, Point, Vector};

#[test]
fn test_intersect_curve_curve_line_circle() {
    use crate::{Circle, Line};
    let tol = Tolerance::DEFAULT;
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let hits = search_curve_curve(&line, &circle, tol);
    assert_eq!(hits.len(), 2);
    assert!((hits[0].first_parameter.unwrap() + 2.0).abs() < 1e-9);
    assert!((hits[1].first_parameter.unwrap() - 2.0).abs() < 1e-9);
    // Cross-check against the exact intersection.
    match line.intersect_exact(&circle, tol) {
        LineCircleIntersection::Points((t1, _), (t2, _)) => {
            assert!((hits[0].first_parameter.unwrap() - t1).abs() < 1e-9);
            assert!((hits[1].first_parameter.unwrap() - t2).abs() < 1e-9);
        }
        _ => panic!("exact intersection expected two points"),
    }
}

/// A deterministic stream of numbers in `[0, 1)` (xorshift).
struct Random(u64);

impl Random {
    fn unit(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A number in `[lo, hi)` rounded to two decimals.
    fn rounded(&mut self, lo: f64, hi: f64) -> f64 {
        ((lo + (hi - lo) * self.unit()) * 100.0).round() / 100.0
    }

    fn unit_vector(&mut self) -> Vector {
        loop {
            let v = Vector::new(
                2.0 * self.unit() - 1.0,
                2.0 * self.unit() - 1.0,
                2.0 * self.unit() - 1.0,
            );
            if let Some(unit) = v.normalized()
                && v.magnitude() > 0.1
            {
                return unit;
            }
        }
    }
}

/// A plane for the random conic pairs: the XY plane, or a random
/// frame centred within 5 of the origin.
struct Frame {
    origin: Point,
    normal: Vector,
    x: Vector,
    y: Vector,
}

impl Frame {
    fn new(random: &mut Random, tilted: bool) -> Frame {
        if !tilted {
            return Frame {
                origin: Point::ORIGIN,
                normal: Vector::Z,
                x: Vector::X,
                y: Vector::Y,
            };
        }
        let origin = Point::ORIGIN + random.unit_vector() * (5.0 * random.unit());
        let normal = random.unit_vector();
        let x = loop {
            let w = random.unit_vector();
            if let Some(x) = (w - normal * w.dot(normal)).normalized() {
                break x;
            }
        };
        Frame {
            origin,
            normal,
            x,
            y: normal.cross(x),
        }
    }

    fn point(&self, u: f64, v: f64) -> Point {
        self.origin + self.x * u + self.y * v
    }

    fn direction(&self, angle: f64) -> Vector {
        self.x * angle.cos() + self.y * angle.sin()
    }

    fn circle(&self, random: &mut Random) -> crate::geometry::curves::ElementaryCurve {
        let center = self.point(random.rounded(-3.0, 3.0), random.rounded(-3.0, 3.0));
        let radius = random.rounded(0.2, 4.0);
        crate::Circle::new(center, self.normal, radius)
            .unwrap()
            .into()
    }

    /// An ellipse, at the origin along the frame's axes when
    /// `centred`, else anywhere at any angle.
    fn ellipse(
        &self,
        random: &mut Random,
        centred: bool,
    ) -> crate::geometry::curves::ElementaryCurve {
        let major = random.rounded(0.2, 4.0);
        let minor = random.rounded(0.1, major);
        let (center, axis) = if centred {
            (self.origin, self.x)
        } else {
            let center = self.point(random.rounded(-2.0, 2.0), random.rounded(-2.0, 2.0));
            (
                center,
                self.direction(random.rounded(0.0, std::f64::consts::PI)),
            )
        };
        crate::Ellipse::new(center, self.normal, axis, major, minor)
            .unwrap()
            .into()
    }

    fn line(&self, random: &mut Random) -> Line {
        let through = self.point(random.rounded(-3.0, 3.0), random.rounded(-3.0, 3.0));
        Line::new(
            through,
            self.direction(random.rounded(0.0, std::f64::consts::PI)),
        )
        .unwrap()
    }
}

/// Whether `intersect_curve_curve` finds, within `tol.confusion`, every
/// transversal hit of the exact solver on `first` and `second` whose
/// parameter on `first` lies in `window`, and only points of the exact
/// solution; the number of transversal hits compared.
fn finds_exact_hits<A: ParametricCurve>(
    first: &A,
    exact_first: &crate::geometry::curves::ElementaryCurve,
    second: &crate::geometry::curves::ElementaryCurve,
    window: (f64, f64),
    tol: Tolerance,
) -> Result<usize, String> {
    let exact = crate::Intersect::intersect(exact_first, second, tol);
    let exact = exact.hits();
    if exact.is_empty() {
        return Ok(0);
    }
    let numeric = intersect_curve_curve(first, second, tol);
    let inside = |u: f64| u > window.0 + 1e-6 && u < window.1 - 1e-6;
    let transversal: Vec<_> = exact
        .iter()
        .filter(|h| h.multiplicity == 1 && h.first_parameter.is_some_and(inside))
        .collect();
    for hit in &transversal {
        if !numeric
            .hits()
            .iter()
            .any(|n| n.point.distance(hit.point) <= tol.confusion)
        {
            return Err(format!(
                "missed {hit:?}: exact {exact:?}, numeric {numeric:?}"
            ));
        }
    }
    for n in numeric.hits() {
        if !exact.iter().any(|h| h.point.distance(n.point) <= 1e-6) {
            return Err(format!("{n:?} is not exact: exact {exact:?}"));
        }
    }
    Ok(transversal.len())
}

/// Compares the numeric and exact hits of `pairs` random coplanar
/// conic pairs of each kind (circle and centred ellipse, two ellipses,
/// a line or a segment of it and a circle or an ellipse), in the XY
/// plane or in random frames; the failures and the number of
/// transversal hits compared.
fn compare_conic_pairs(seed: u64, pairs: usize, tilted: bool) -> (Vec<String>, usize) {
    use crate::geometry::curves::ElementaryCurve;
    let tol = Tolerance::DEFAULT;
    let mut random = Random(seed);
    let (mut failures, mut compared) = (Vec::new(), 0);
    let mut check = |outcome: Result<usize, String>| match outcome {
        Ok(hits) => compared += hits,
        Err(failure) => failures.push(failure),
    };
    let unbounded = (f64::NEG_INFINITY, f64::INFINITY);
    for k in 0..pairs {
        let frame = Frame::new(&mut random, tilted);
        let circle = frame.circle(&mut random);
        let centred = frame.ellipse(&mut random, true);
        check(finds_exact_hits(&circle, &circle, &centred, unbounded, tol));
        let (first, second) = (
            frame.ellipse(&mut random, false),
            frame.ellipse(&mut random, false),
        );
        check(finds_exact_hits(&first, &first, &second, unbounded, tol));
        let conic = if k % 2 == 0 { circle } else { second };
        let line = frame.line(&mut random);
        let exact_line = ElementaryCurve::from(line);
        check(finds_exact_hits(&line, &exact_line, &conic, unbounded, tol));
        let range = (random.rounded(-6.0, 0.0), random.rounded(0.01, 6.0));
        let segment = crate::TrimmedCurve::new(line, range.0, range.1).unwrap();
        check(finds_exact_hits(&segment, &exact_line, &conic, range, tol));
    }
    (failures, compared)
}

#[test]
fn test_intersect_curve_curve_finds_close_circle_ellipse_crossings() {
    // Issue #13: each pair has a crossing the grid alone missed, close
    // to another crossing or to a near touch.
    let tol = Tolerance::DEFAULT;
    let cases = [
        ((3.44, 1.03), (-1.77, 0.84), 1.9, (-3.42685, -0.08997)),
        ((1.04, 0.22), (0.73, -1.52), 2.31, (-1.03236, -0.02662)),
        ((2.56, 1.21), (-0.57, 0.97), 2.21, (-0.18799, -1.20673)),
        ((2.72, 0.66), (0.86, 1.86), 2.49, (1.45751, -0.55725)),
    ];
    for ((major, minor), (x, y), radius, missed) in cases {
        let ellipse =
            crate::Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, major, minor).unwrap();
        let circle = crate::Circle::new(Point::new(x, y, 0.0), Vector::Z, radius).unwrap();
        let exact = crate::Intersect::intersect(&circle, &ellipse, tol);
        let numeric = intersect_curve_curve(&circle, &ellipse, tol);
        assert_eq!(numeric.hits().len(), exact.hits().len(), "{numeric:?}");
        for hit in exact.hits() {
            assert!(
                numeric
                    .hits()
                    .iter()
                    .any(|n| n.point.distance(hit.point) <= tol.confusion),
                "missed {hit:?}: {numeric:?}"
            );
        }
        let missed = Point::new(missed.0, missed.1, 0.0);
        assert!(
            numeric
                .hits()
                .iter()
                .any(|n| n.point.distance(missed) < 1e-5)
        );
    }
}

#[test]
fn test_intersect_curve_curve_finds_every_exact_conic_hit() {
    for (seed, tilted) in [
        (0x9E37_79B9_7F4A_7C15, false),
        (0x2545_F491_4F6C_DD1D, true),
    ] {
        let (failures, compared) = compare_conic_pairs(seed, 1000, tilted);
        assert!(compared >= 3000, "only {compared} hits compared");
        assert!(
            failures.is_empty(),
            "{} failures, first: {}",
            failures.len(),
            failures[0]
        );
    }
}
