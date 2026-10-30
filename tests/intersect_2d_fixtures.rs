//! Golden fixtures for 2D curve-curve intersections: one case per pair
//! family, with the hits and overlaps worked out by hand (or by an
//! independent derivation, noted in each case) and replayed in both
//! argument orders. `tests/fixtures/intersect_2d.json` holds the curves
//! and the expected result; the test checks the method, the completeness,
//! every hit's point and multiplicity, and every overlap's ranges.

mod common;

use common::load;
use geomcore::curves::{Curve2D, TrimmedCurve2D};
use geomcore::prelude::*;
use geomcore::{
    BSplineCurve2D, Circle2D, CurveCurveIntersection2D, Ellipse2D, Frame2D, Hyperbola2D,
    IntersectionMethod, Line2D, Parabola2D, Point2D, Tolerance, Vector2D,
};
use serde::Deserialize;

const TOL: Tolerance = Tolerance::DEFAULT;
const EXACT_EPS: f64 = 1e-9;
const NUMERIC_EPS: f64 = 1e-7;

#[derive(Deserialize)]
struct Fixture {
    cases: Vec<Case>,
}

#[derive(Deserialize)]
struct Case {
    name: String,
    first: CurveSpec,
    second: CurveSpec,
    expect: Expect,
}

#[derive(Deserialize)]
struct CurveSpec {
    shape: Shape,
    trim: Option<[f64; 2]>,
}

#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
enum Shape {
    Line {
        origin: [f64; 2],
        direction: [f64; 2],
    },
    Circle {
        center: [f64; 2],
        radius: f64,
        #[serde(default)]
        clockwise: bool,
    },
    Ellipse {
        center: [f64; 2],
        x_dir: [f64; 2],
        major: f64,
        minor: f64,
        #[serde(default)]
        clockwise: bool,
    },
    Parabola {
        apex: [f64; 2],
        x_dir: [f64; 2],
        focal: f64,
        #[serde(default)]
        clockwise: bool,
    },
    Hyperbola {
        center: [f64; 2],
        x_dir: [f64; 2],
        major: f64,
        minor: f64,
        #[serde(default)]
        clockwise: bool,
    },
    /// A degree-1 B-spline through the points, parametrized 0, 1, 2, ...
    Polyline { points: Vec<[f64; 2]> },
    /// The exact B-spline of a circle or an ellipse.
    BsplineOf { of: Box<Shape> },
}

#[derive(Deserialize)]
struct Expect {
    method: String,
    complete: bool,
    hits: Vec<ExpectedHit>,
    #[serde(default)]
    overlaps: Vec<ExpectedOverlap>,
}

#[derive(Deserialize)]
struct ExpectedHit {
    point: [f64; 2],
    multiplicity: u32,
}

#[derive(Deserialize)]
struct ExpectedOverlap {
    first: [f64; 2],
    second: [f64; 2],
    same_direction: bool,
}

fn p(a: [f64; 2]) -> Point2D {
    Point2D::new(a[0], a[1])
}

fn v(a: [f64; 2]) -> Vector2D {
    Vector2D::new(a[0], a[1])
}

/// The frame at `origin` with its x axis along `x_dir`, and its y axis
/// turned a quarter turn left, or right when `clockwise`.
fn frame(origin: [f64; 2], x_dir: [f64; 2], clockwise: bool) -> Frame2D {
    let x = v(x_dir).normalized().unwrap();
    let y = if clockwise {
        Vector2D::new(x.y, -x.x)
    } else {
        Vector2D::new(-x.y, x.x)
    };
    Frame2D::new(p(origin), x, y).unwrap()
}

fn build_shape(shape: &Shape) -> Curve2D {
    match shape {
        Shape::Line { origin, direction } => Line2D::new(p(*origin), v(*direction)).unwrap().into(),
        Shape::Circle {
            center,
            radius,
            clockwise,
        } => Circle2D::from_frame(frame(*center, [1.0, 0.0], *clockwise), *radius)
            .unwrap()
            .into(),
        Shape::Ellipse {
            center,
            x_dir,
            major,
            minor,
            clockwise,
        } => Ellipse2D::from_frame(frame(*center, *x_dir, *clockwise), *major, *minor)
            .unwrap()
            .into(),
        Shape::Parabola {
            apex,
            x_dir,
            focal,
            clockwise,
        } => Parabola2D::from_frame(frame(*apex, *x_dir, *clockwise), *focal)
            .unwrap()
            .into(),
        Shape::Hyperbola {
            center,
            x_dir,
            major,
            minor,
            clockwise,
        } => Hyperbola2D::from_frame(frame(*center, *x_dir, *clockwise), *major, *minor)
            .unwrap()
            .into(),
        Shape::Polyline { points } => {
            let last = points.len() - 1;
            let mut multiplicities = vec![1; points.len()];
            multiplicities[0] = 2;
            multiplicities[last] = 2;
            BSplineCurve2D::new(
                1,
                points.iter().copied().map(p).collect(),
                (0..points.len()).map(|i| i as f64).collect(),
                multiplicities,
                false,
            )
            .unwrap()
            .into()
        }
        Shape::BsplineOf { of } => match build_shape(of) {
            Curve2D::Circle(c) => c.to_bspline().into(),
            Curve2D::Ellipse(e) => e.to_bspline().into(),
            _ => panic!("only a circle or an ellipse converts whole"),
        },
    }
}

/// A fixture curve: whole, or restricted to a parameter range.
enum Operand {
    Whole(Curve2D),
    Trimmed(TrimmedCurve2D<Curve2D>),
}

fn build(spec: &CurveSpec) -> Operand {
    let curve = build_shape(&spec.shape);
    match spec.trim {
        Some([first, last]) => Operand::Trimmed(TrimmedCurve2D::new(curve, first, last).unwrap()),
        None => Operand::Whole(curve),
    }
}

impl Operand {
    fn point_at(&self, parameter: f64) -> Point2D {
        match self {
            Operand::Whole(c) => c.eval_point(parameter),
            Operand::Trimmed(c) => c.eval_point(parameter),
        }
    }

    fn intersect(&self, other: &Operand) -> CurveCurveIntersection2D {
        match (self, other) {
            (Operand::Whole(a), Operand::Whole(b)) => a.intersect(b, TOL),
            (Operand::Whole(a), Operand::Trimmed(b)) => a.intersect(b, TOL),
            (Operand::Trimmed(a), Operand::Whole(b)) => a.intersect(b, TOL),
            (Operand::Trimmed(a), Operand::Trimmed(b)) => a.intersect(b, TOL),
        }
    }
}

fn near(a: [f64; 2], b: [f64; 2], eps: f64) -> bool {
    (a[0] - b[0]).abs() <= eps && (a[1] - b[1]).abs() <= eps
}

/// Checks `result` against `expect`; `swapped` says the result reads the
/// pair in the reverse of the fixture's order.
#[track_caller]
fn check(
    name: &str,
    result: &CurveCurveIntersection2D,
    expect: &Expect,
    operands: (&Operand, &Operand),
    swapped: bool,
) {
    let eps = if expect.method == "numeric" {
        NUMERIC_EPS
    } else {
        EXACT_EPS
    };
    let label = format!("{name} (swapped: {swapped})");
    let method = match result.method() {
        IntersectionMethod::Exact => "exact",
        IntersectionMethod::Numeric => "numeric",
    };
    assert_eq!(method, expect.method, "{label}: method");
    assert_eq!(
        result.is_complete(),
        expect.complete,
        "{label}: completeness"
    );

    let hits = result.hits();
    assert_eq!(hits.len(), expect.hits.len(), "{label}: hits {hits:?}");
    let mut taken = vec![false; hits.len()];
    for want in &expect.hits {
        let found = hits.iter().enumerate().position(|(i, hit)| {
            !taken[i]
                && near([hit.point.x, hit.point.y], want.point, eps)
                && hit.multiplicity == want.multiplicity
        });
        let i = found.unwrap_or_else(|| panic!("{label}: no hit at {:?} in {hits:?}", want.point));
        taken[i] = true;
    }
    // A stored parameter lands on the hit's point, on its own curve.
    let (first, second) = operands;
    for hit in hits {
        for (parameter, curve) in [(hit.first_parameter, first), (hit.second_parameter, second)] {
            if let Some(t) = parameter {
                let at = curve.point_at(t);
                assert!(
                    at.distance(hit.point) <= 1e-7,
                    "{label}: parameter {t} is at {at:?}, not {:?}",
                    hit.point
                );
            }
        }
    }

    let overlaps = result.overlaps();
    assert_eq!(
        overlaps.len(),
        expect.overlaps.len(),
        "{label}: overlaps {overlaps:?}"
    );
    for want in &expect.overlaps {
        let (own, other) = if swapped {
            (want.second, want.first)
        } else {
            (want.first, want.second)
        };
        assert!(
            overlaps.iter().any(|o| {
                near([o.first_range.0, o.first_range.1], own, eps)
                    && near([o.second_range.0, o.second_range.1], other, eps)
                    && o.same_direction == want.same_direction
            }),
            "{label}: no overlap {own:?} / {other:?} in {overlaps:?}"
        );
    }
}

#[test]
fn test_2d_intersection_fixtures_in_both_orders() {
    let fixture: Fixture = serde_json::from_value(load("intersect_2d.json")).unwrap();
    assert!(fixture.cases.len() >= 40);
    for case in &fixture.cases {
        let (first, second) = (build(&case.first), build(&case.second));
        let forward = first.intersect(&second);
        check(&case.name, &forward, &case.expect, (&first, &second), false);
        let backward = second.intersect(&first);
        check(&case.name, &backward, &case.expect, (&second, &first), true);
    }
}
