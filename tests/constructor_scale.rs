//! The constructors that decide whether their points are degenerate read
//! the caller's tolerance (#473): scaling the points and `tol.confusion`
//! by one factor accepts or refuses the same configurations (ADR 0051).

use geomcore::{
    ArcConstructionError, Circle, Circle2D, CircleConstructionError, Cone, ConeConstructionError,
    Ellipse, Ellipse2D, EllipseConstructionError, Hyperbola, Hyperbola2D,
    HyperbolaConstructionError, Line, Line2D, Parabola, Parabola2D, ParabolaConstructionError,
    Plane, PlaneConstructionError, Point, Point2D, Tolerance, TrimmedCurve, TrimmedCurve2D, Vector,
    Vector2D,
};
use std::fmt::Debug;

/// Half the default confusion distance: a gap the default tolerance
/// treats as no gap.
const HAIR: f64 = 5e-8;

/// Checks that `build` gives `expected` at the unit scale and at the
/// scales 1e-3 and 1e3, with the confusion distance scaled alike.
fn at_every_scale<E: PartialEq + Debug>(
    expected: Result<(), E>,
    build: impl Fn(f64, Tolerance) -> Result<(), E>,
) {
    let default = Tolerance::DEFAULT;
    for scale in [1.0, 1e-3, 1e3] {
        let tol = Tolerance::new(
            default.confusion * scale,
            default.angular,
            default.parametric,
        )
        .unwrap();
        assert_eq!(build(scale, tol), expected, "at scale {scale}");
    }
}

fn point(scale: f64, x: f64, y: f64, z: f64) -> Point {
    Point::new(x * scale, y * scale, z * scale)
}

fn point_2d(scale: f64, x: f64, y: f64) -> Point2D {
    Point2D::new(x * scale, y * scale)
}

#[test]
fn test_circle_through_three_points_at_every_scale() {
    let circle = |p: [[f64; 3]; 3]| {
        move |s: f64, tol| {
            let [a, b, c] = p.map(|[x, y, z]| point(s, x, y, z));
            Circle::from_three_points(a, b, c, tol).map(drop)
        }
    };
    let circle_2d = |p: [[f64; 2]; 3]| {
        move |s: f64, tol| {
            let [a, b, c] = p.map(|[x, y]| point_2d(s, x, y));
            Circle2D::from_three_points(a, b, c, tol).map(drop)
        }
    };
    let turn = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]];
    let confused = [[1.0, 0.0, 0.0], [1.0 + HAIR, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    let collinear = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, HAIR, 0.0]];
    at_every_scale(Ok(()), circle(turn));
    at_every_scale(
        Err(CircleConstructionError::ConfusedPoints),
        circle(confused),
    );
    at_every_scale(
        Err(CircleConstructionError::CollinearPoints),
        circle(collinear),
    );
    let flat = |p: [[f64; 3]; 3]| p.map(|[x, y, _]| [x, y]);
    at_every_scale(Ok(()), circle_2d(flat(turn)));
    at_every_scale(
        Err(CircleConstructionError::ConfusedPoints),
        circle_2d(flat(confused)),
    );
    at_every_scale(
        Err(CircleConstructionError::CollinearPoints),
        circle_2d(flat(collinear)),
    );
}

#[test]
fn test_arc_through_three_points_at_every_scale() {
    let arc = |p: [[f64; 2]; 3]| {
        move |s: f64, tol| {
            let [a, b, c] = p.map(|[x, y]| point(s, x, y, 0.0));
            TrimmedCurve::arc_through(a, b, c, tol).map(drop)
        }
    };
    let arc_2d = |p: [[f64; 2]; 3]| {
        move |s: f64, tol| {
            let [a, b, c] = p.map(|[x, y]| point_2d(s, x, y));
            TrimmedCurve2D::arc_through(a, b, c, tol).map(drop)
        }
    };
    let half = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]];
    let confused = [[1.0, 0.0], [1.0 + HAIR, 0.0], [-1.0, 0.0]];
    let collinear = [[0.0, 0.0], [1.0, 0.0], [2.0, HAIR]];
    for (points, expected) in [
        (half, Ok(())),
        (confused, Err(ArcConstructionError::ConfusedPoints)),
        (collinear, Err(ArcConstructionError::CollinearPoints)),
    ] {
        at_every_scale(expected, arc(points));
        at_every_scale(expected, arc_2d(points));
    }
}

#[test]
fn test_plane_through_three_points_at_every_scale() {
    let plane = |p: [[f64; 3]; 3]| {
        move |s: f64, tol| {
            let [a, b, c] = p.map(|[x, y, z]| point(s, x, y, z));
            Plane::from_three_points(a, b, c, tol).map(drop)
        }
    };
    at_every_scale(
        Ok(()),
        plane([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    );
    at_every_scale(
        Err(PlaneConstructionError::ConfusedPoints),
        plane([[0.0, 0.0, 0.0], [HAIR, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    );
    at_every_scale(
        Err(PlaneConstructionError::CollinearPoints),
        plane([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, HAIR, 0.0]]),
    );
}

#[test]
fn test_ellipse_from_center_and_points_at_every_scale() {
    let ellipse = |s1: [f64; 2], s2: [f64; 2]| {
        move |s: f64, tol| {
            let at = |[x, y]: [f64; 2]| point(s, x, y, 0.0);
            Ellipse::from_center_and_points(at([0.0, 0.0]), at(s1), at(s2), tol).map(drop)
        }
    };
    let ellipse_2d = |s1: [f64; 2], s2: [f64; 2]| {
        move |s: f64, tol| {
            let at = |[x, y]: [f64; 2]| point_2d(s, x, y);
            Ellipse2D::from_center_and_points(at([0.0, 0.0]), at(s1), at(s2), tol).map(drop)
        }
    };
    for (s1, s2, expected) in [
        ([3.0, 0.0], [0.0, 1.5], Ok(())),
        (
            [HAIR, 0.0],
            [0.0, 1.5],
            Err(EllipseConstructionError::NullAxis),
        ),
        (
            [3.0, 0.0],
            [1.0, HAIR],
            Err(EllipseConstructionError::InvertedAxis),
        ),
        (
            [3.0, 0.0],
            [0.0, 4.0],
            Err(EllipseConstructionError::InvertedAxis),
        ),
    ] {
        at_every_scale(expected, ellipse(s1, s2));
        at_every_scale(expected, ellipse_2d(s1, s2));
    }
}

#[test]
fn test_hyperbola_from_center_and_points_at_every_scale() {
    let hyperbola = |s1: [f64; 2], s2: [f64; 2]| {
        move |s: f64, tol| {
            let at = |[x, y]: [f64; 2]| point(s, x, y, 0.0);
            Hyperbola::from_center_and_points(at([0.0, 0.0]), at(s1), at(s2), tol).map(drop)
        }
    };
    let hyperbola_2d = |s1: [f64; 2], s2: [f64; 2]| {
        move |s: f64, tol| {
            let at = |[x, y]: [f64; 2]| point_2d(s, x, y);
            Hyperbola2D::from_center_and_points(at([0.0, 0.0]), at(s1), at(s2), tol).map(drop)
        }
    };
    for (s1, s2, expected) in [
        ([2.0, 0.0], [0.0, 1.0], Ok(())),
        (
            [HAIR, 0.0],
            [0.0, 1.0],
            Err(HyperbolaConstructionError::ConfusedPoints),
        ),
        (
            [2.0, 0.0],
            [3.0, HAIR],
            Err(HyperbolaConstructionError::CollinearPoints),
        ),
    ] {
        at_every_scale(expected, hyperbola(s1, s2));
        at_every_scale(expected, hyperbola_2d(s1, s2));
    }
}

#[test]
fn test_cone_through_two_sections_at_every_scale() {
    let cone = |height: f64, r1: f64, r2: f64| {
        move |s: f64, tol| {
            let (p1, p2) = (point(s, 0.0, 0.0, 0.0), point(s, 0.0, 0.0, height));
            Cone::from_two_points_and_radii(p1, p2, r1 * s, r2 * s, tol).map(drop)
        }
    };
    at_every_scale(Ok(()), cone(4.0, 2.0, 1.0));
    at_every_scale(Err(ConeConstructionError::NullAxis), cone(HAIR, 2.0, 1.0));
    at_every_scale(Err(ConeConstructionError::NullAngle), cone(4.0, 2.0, 2.0));
}

#[test]
fn test_parabola_from_focus_and_directrix_at_every_scale() {
    let parabola = |focus_height: f64| {
        move |s: f64, tol| {
            let directrix = Line::new(point(s, 0.0, 0.0, -1.0), Vector::X).unwrap();
            Parabola::from_focus_and_directrix(point(s, 0.0, 0.0, focus_height), &directrix, tol)
                .map(drop)
        }
    };
    let parabola_2d = |focus_height: f64| {
        move |s: f64, tol| {
            let directrix = Line2D::new(point_2d(s, 0.0, -1.0), Vector2D::X).unwrap();
            Parabola2D::from_focus_and_directrix(point_2d(s, 0.0, focus_height), &directrix, tol)
                .map(drop)
        }
    };
    at_every_scale(Ok(()), parabola(1.0));
    at_every_scale(Ok(()), parabola_2d(1.0));
    let on_directrix = Err(ParabolaConstructionError::NonPositiveFocal);
    at_every_scale(on_directrix, parabola(-1.0 + HAIR));
    at_every_scale(on_directrix, parabola_2d(-1.0 + HAIR));
}
