use super::*;
use crate::{
    BSplineCurve, Completeness, CurveCurveIntersection, IntersectExact, Line,
    LineCircleIntersection, Point, Vector,
};

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

/// The exact nine-pole rational B-spline circle of radius `radius` about
/// the origin, in the XY plane.
fn nine_pole_circle(radius: f64) -> BSplineCurve {
    let corners = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
    ];
    let poles = corners
        .into_iter()
        .map(|(x, y)| Point::new(radius * x, radius * y, 0.0))
        .collect();
    let weights = (0..corners.len())
        .map(|k| {
            if k % 2 == 1 {
                std::f64::consts::FRAC_1_SQRT_2
            } else {
                1.0
            }
        })
        .collect();
    BSplineCurve::new_rational(
        2,
        poles,
        weights,
        vec![0.0, 1.0, 2.0, 3.0, 4.0],
        vec![3, 2, 2, 2, 3],
        false,
    )
    .unwrap()
}

/// The B-spline segment from `from` to `to`.
fn segment(from: Point, to: Point) -> BSplineCurve {
    BSplineCurve::new(1, vec![from, to], vec![0.0, 1.0], vec![2, 2], false).unwrap()
}

/// Whether `result` shows `point` as a hit, or as lying in a region where
/// hits may be missing.
fn shows_point(result: &CurveCurveIntersection, point: Point, tol: Tolerance) -> bool {
    result
        .hits()
        .iter()
        .any(|hit| hit.point.distance(point) <= tol.confusion)
        || match result.completeness() {
            Completeness::Incomplete { regions, .. } => regions.iter().any(|region| {
                let (lo, hi) = (region.min_corner, region.max_corner);
                (lo.x..=hi.x).contains(&point.x)
                    && (lo.y..=hi.y).contains(&point.y)
                    && (lo.z..=hi.z).contains(&point.z)
            }),
            _ => false,
        }
}

#[test]
fn test_numeric_keeps_crossings_closer_than_sqrt_confusion() {
    // The two crossings of a circle and the chord y = s / 2 are 1.73 s
    // apart: from s = 1e-4 on they lie within sqrt(confusion) of each
    // other, the merge radius that used to drop one of them.
    let tol = Tolerance::DEFAULT;
    for s in [1e-3, 1e-4, 1e-5, 1e-6] {
        let circle = nine_pole_circle(s);
        let chord = segment(
            Point::new(-2.0 * s, 0.5 * s, 0.0),
            Point::new(2.0 * s, 0.5 * s, 0.0),
        );
        let result = intersect_curve_curve(&circle, &chord, tol);
        assert!(result.is_complete(), "s = {s}: {result:?}");
        assert_eq!(result.hits().len(), 2, "s = {s}: {result:?}");
        let half_chord = 3f64.sqrt() / 2.0 * s;
        for x in [-half_chord, half_chord] {
            let crossing = Point::new(x, 0.5 * s, 0.0);
            assert!(
                result
                    .hits()
                    .iter()
                    .any(|hit| hit.point.distance(crossing) <= tol.confusion),
                "s = {s}: no hit at {crossing:?} in {result:?}"
            );
        }
        assert!(result.hits().iter().all(|hit| hit.multiplicity == 1));
    }
}

#[test]
fn test_numeric_keeps_crossings_closer_than_sqrt_confusion_at_small_confusion() {
    // At confusion 1e-12 the merge radius is 1e-6: the crossings at scales
    // 3e-7 and 1e-9, 5e-7 and 2e-9 apart, fell inside it.
    let tol = Tolerance {
        confusion: 1e-12,
        ..Tolerance::DEFAULT
    };
    for s in [3e-7, 1e-9] {
        let circle = nine_pole_circle(s);
        let chord = segment(
            Point::new(-2.0 * s, 0.5 * s, 0.0),
            Point::new(2.0 * s, 0.5 * s, 0.0),
        );
        let result = intersect_curve_curve(&circle, &chord, tol);
        let half_chord = 3f64.sqrt() / 2.0 * s;
        for x in [-half_chord, half_chord] {
            assert!(
                shows_point(&result, Point::new(x, 0.5 * s, 0.0), tol),
                "s = {s}: {result:?}"
            );
        }
        if result.is_complete() {
            assert_eq!(result.hits().len(), 2, "s = {s}: {result:?}");
        }
    }
}

#[test]
fn test_numeric_near_tangent_line_keeps_its_two_hits() {
    // The line y = 1 - 1e-6 meets the unit circle at x = +-1.414e-3, a
    // near-tangent pair whose hits the merge rule must not change.
    let tol = Tolerance::DEFAULT;
    let circle = crate::Circle::new(Point::ORIGIN, Vector::Z, 1.0)
        .unwrap()
        .to_bspline();
    let chord = segment(
        Point::new(-2.0, 1.0 - 1e-6, 0.0),
        Point::new(2.0, 1.0 - 1e-6, 0.0),
    );
    let result = intersect_curve_curve(&circle, &chord, tol);
    assert!(result.is_complete(), "{result:?}");
    assert_eq!(result.hits().len(), 2, "{result:?}");
    for x in [-1.414_213_2e-3, 1.414_213_2e-3] {
        assert!(
            result
                .hits()
                .iter()
                .any(|hit| (hit.point.x - x).abs() < 1e-8
                    && (hit.point.y - (1.0 - 1e-6)).abs() < 1e-8
                    && hit.multiplicity == 1),
            "no hit at x = {x}: {result:?}"
        );
    }
}

#[test]
fn test_near_tangent_line_matches_closed_form_with_two_crossings() {
    // The line y = 1 - 1e-6 crosses the unit circle twice, 2.8e-3 apart:
    // two crossings of multiplicity 1, as the closed form gives.
    let tol = Tolerance::DEFAULT;
    let circle = crate::Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let line = Line::new(Point::new(0.0, 1.0 - 1e-6, 0.0), Vector::X).unwrap();
    let chord = segment(
        Point::new(-2.0, 1.0 - 1e-6, 0.0),
        Point::new(2.0, 1.0 - 1e-6, 0.0),
    );
    let closed = crate::Intersect::intersect(&circle, &line, tol);
    assert_eq!(closed.hits().len(), 2, "{closed:?}");
    assert!(closed.hits().iter().all(|hit| hit.multiplicity == 1));
    let numeric = intersect_curve_curve(&circle.to_bspline(), &chord, tol);
    assert_agrees_with_closed_form(&closed, &numeric);
}

#[test]
fn test_shallow_crossing_of_two_lines_stays_one_crossing() {
    // Two B-spline lines cross at the origin at 1e-4 rad. Moving either by
    // confusion never makes them touch, so the hit is a crossing.
    let tol = Tolerance::DEFAULT;
    let first = segment(Point::new(-1.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0));
    let second = segment(Point::new(-1.0, -1e-4, 0.0), Point::new(1.0, 1e-4, 0.0));
    let result = intersect_curve_curve(&first, &second, tol);
    assert!(result.is_complete(), "{result:?}");
    assert_eq!(result.hits().len(), 1, "{result:?}");
    assert_eq!(result.hits()[0].multiplicity, 1, "{result:?}");
}

#[test]
fn test_circle_tangent_to_line_at_three_scales_is_one_contact() {
    // A circle of radius r about the origin touches y = r at (0, r). The
    // curvature 1 / r sets how far the two can part within confusion, so the
    // contact is a touch at every scale: multiplicity 2, complete.
    let tol = Tolerance::DEFAULT;
    for radius in [1.0, 1e-3, 1e3] {
        let circle = crate::Circle::new(Point::ORIGIN, Vector::Z, radius).unwrap();
        let line = Line::new(Point::new(0.0, radius, 0.0), Vector::X).unwrap();
        let chord = segment(
            Point::new(-2.0 * radius, radius, 0.0),
            Point::new(2.0 * radius, radius, 0.0),
        );
        let closed = crate::Intersect::intersect(&circle, &line, tol);
        assert_eq!(closed.hits().len(), 1, "radius {radius}: {closed:?}");
        assert_eq!(closed.hits()[0].multiplicity, 2, "radius {radius}");
        let numeric = intersect_curve_curve(&nine_pole_circle(radius), &chord, tol);
        assert_eq!(numeric.hits().len(), 1, "radius {radius}: {numeric:?}");
        assert_eq!(
            numeric.hits()[0].multiplicity,
            2,
            "radius {radius}: {numeric:?}"
        );
        assert!(numeric.is_complete(), "radius {radius}: {numeric:?}");
        assert!(
            numeric.hits()[0]
                .point
                .distance(Point::new(0.0, radius, 0.0))
                <= 1e-6 * radius,
            "radius {radius}: {numeric:?}"
        );
    }
}

/// The quadratic B-spline `y = x^2` over `x` in `[-1, 1]`, the Bézier
/// through the control point `(0, -1)`.
fn parabola_bspline() -> BSplineCurve {
    BSplineCurve::new(
        2,
        vec![
            Point::new(-1.0, 1.0, 0.0),
            Point::new(0.0, -1.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
        ],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
    )
    .unwrap()
}

/// Asserts that the numeric result on B-splines gives the closed form's
/// hits: the same points and multiplicities, and complete.
#[track_caller]
fn assert_agrees_with_closed_form(
    closed: &CurveCurveIntersection,
    numeric: &CurveCurveIntersection,
) {
    assert_eq!(
        closed.hits().len(),
        numeric.hits().len(),
        "{closed:?} vs {numeric:?}"
    );
    for hit in closed.hits() {
        assert!(
            numeric.hits().iter().any(|n| {
                n.point.distance(hit.point) < 1e-7 && n.multiplicity == hit.multiplicity
            }),
            "closed {hit:?} missing from {numeric:?}"
        );
    }
    assert!(numeric.is_complete(), "{numeric:?}");
}

#[test]
fn test_parabola_contact_with_flat_line_matches_closed_form() {
    // The parabola y = x^2 touches y = -1e-9 at the origin: its distance
    // there is within the default confusion (1e-7), so one contact.
    let tol = Tolerance::DEFAULT;
    let parabola = crate::Parabola::new(Point::ORIGIN, Vector::Z, Vector::Y, 0.25).unwrap();
    let line = Line::new(Point::new(-2.0, -1e-9, 0.0), Vector::X).unwrap();
    let flat = segment(Point::new(-2.0, -1e-9, 0.0), Point::new(2.0, -1e-9, 0.0));
    let closed = crate::Intersect::intersect(&parabola, &line, tol);
    assert_eq!(closed.hits().len(), 1, "{closed:?}");
    assert_eq!(closed.hits()[0].multiplicity, 2);
    let numeric = intersect_curve_curve(&parabola_bspline(), &flat, tol);
    assert_eq!(numeric.hits().len(), 1, "{numeric:?}");
    assert_eq!(numeric.hits()[0].multiplicity, 2, "{numeric:?}");
    assert!(
        numeric.hits()[0].point.distance(Point::ORIGIN) < 1e-7,
        "{numeric:?}"
    );
    assert!(numeric.is_complete(), "{numeric:?}");
}

#[test]
fn test_parabola_contact_vanishes_at_confusion_1e_12() {
    // At confusion 1e-12 the gap 1e-9 is a separation: no hit, complete.
    let tol = Tolerance {
        confusion: 1e-12,
        ..Tolerance::DEFAULT
    };
    let flat = segment(Point::new(-2.0, -1e-9, 0.0), Point::new(2.0, -1e-9, 0.0));
    let numeric = intersect_curve_curve(&parabola_bspline(), &flat, tol);
    assert!(numeric.hits().is_empty(), "{numeric:?}");
    assert!(numeric.is_complete(), "{numeric:?}");
}

#[test]
fn test_circle_tangent_to_line_matches_closed_form() {
    let tol = Tolerance::DEFAULT;
    let circle = crate::Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let line = Line::new(Point::new(0.0, 1.0, 0.0), Vector::X).unwrap();
    let chord = segment(Point::new(-2.0, 1.0, 0.0), Point::new(2.0, 1.0, 0.0));
    let closed = crate::Intersect::intersect(&circle, &line, tol);
    let numeric = intersect_curve_curve(&circle.to_bspline(), &chord, tol);
    assert_agrees_with_closed_form(&closed, &numeric);
}

#[test]
fn test_circles_tangent_externally_and_internally_match_closed_form() {
    let tol = Tolerance::DEFAULT;
    let external = (
        crate::Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap(),
        crate::Circle::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 1.0).unwrap(),
    );
    let internal = (
        crate::Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap(),
        crate::Circle::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 1.0).unwrap(),
    );
    for (first, second) in [external, internal] {
        let closed = crate::Intersect::intersect(&first, &second, tol);
        let numeric = intersect_curve_curve(&first.to_bspline(), &second.to_bspline(), tol);
        assert_agrees_with_closed_form(&closed, &numeric);
    }
}

#[test]
fn test_ellipse_tangent_to_line_matches_closed_form() {
    let tol = Tolerance::DEFAULT;
    // The ellipse with semi-axes 2 and 1 touches y = 1 at (0, 1).
    let ellipse = crate::Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let line = Line::new(Point::new(0.0, 1.0, 0.0), Vector::X).unwrap();
    let chord = segment(Point::new(-3.0, 1.0, 0.0), Point::new(3.0, 1.0, 0.0));
    let closed = crate::Intersect::intersect(&ellipse, &line, tol);
    let numeric = intersect_curve_curve(&ellipse.to_bspline(), &chord, tol);
    assert_agrees_with_closed_form(&closed, &numeric);
}

#[test]
fn test_ellipse_and_circle_touching_at_two_points_give_both_points() {
    // The ellipse with semi-axes 2 and 1 and the unit circle touch at
    // (0, 1) and (0, -1).
    let tol = Tolerance::DEFAULT;
    let ellipse = crate::Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let circle = crate::Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let closed = crate::Intersect::intersect(&circle, &ellipse, tol);
    let numeric = intersect_curve_curve(&circle.to_bspline(), &ellipse.to_bspline(), tol);
    for y in [1.0, -1.0] {
        let contact = Point::new(0.0, y, 0.0);
        assert!(
            numeric
                .hits()
                .iter()
                .any(|hit| hit.point.distance(contact) < 1e-7 && hit.multiplicity == 2),
            "no contact at {contact:?}: {numeric:?}"
        );
    }
    assert_eq!(numeric.hits().len(), 2, "{numeric:?}");
    assert!(
        numeric.is_complete(),
        "closed {closed:?}, numeric {numeric:?}"
    );
}
