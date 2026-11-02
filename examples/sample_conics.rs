//! Samples a circle and an ellipse over a full turn, and demonstrates
//! analytic curve-on-surface parametrization by placing a circle on a
//! cylinder and reading off its `(u, v)` image.
//!
//! Run with `cargo run --example sample_conics`.

use geomcore::prelude::*;
use geomcore::{Circle, Curve2D, Cylinder, Ellipse, Point, Tolerance, Vector};
use std::f64::consts::TAU;

const SAMPLE_COUNT: usize = 100;

fn main() {
    let params: Vec<f64> = (0..SAMPLE_COUNT)
        .map(|i| TAU * i as f64 / SAMPLE_COUNT as f64)
        .collect();

    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let circle_points = circle.eval_points(&params);

    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let ellipse_points = ellipse.eval_points(&params);

    println!("Circle (radius 2.0), first/last 3 of {SAMPLE_COUNT} points over [0, TAU):");
    print_head_and_tail(&params, &circle_points);

    println!(
        "\nEllipse (major 3.0, minor 1.5), first/last 3 of {SAMPLE_COUNT} points over [0, TAU):"
    );
    print_head_and_tail(&params, &ellipse_points);

    // Curve-on-surface parametrization: a circle coaxial with a cylinder
    // has an exact UV curve, a horizontal line at v = height in (u, v).
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let section = Circle::new(Point::new(0.0, 0.0, 3.0), Vector::Z, 2.0).unwrap();
    let uv_curve = cylinder
        .curve_3d_to_uv(&section.into(), Tolerance::DEFAULT)
        .unwrap();
    match &uv_curve {
        Curve2D::Line(line) => {
            println!(
                "\nCircle-on-cylinder UV curve: (u, v) = ({:.3}, {:.3}) + t * ({:.3}, {:.3})",
                line.origin().x,
                line.origin().y,
                line.direction().x,
                line.direction().y
            );
        }
        _ => unreachable!("a coaxial circle on a cylinder always maps to a line"),
    }
    let q = uv_curve.eval_point(1.0);
    println!("uv_curve(1.0) = ({:.6}, {:.6})", q.x, q.y);
}

fn print_head_and_tail(params: &[f64], points: &[Point]) {
    let n = points.len();
    for i in 0..3 {
        println!("  u={:.4}  ->  {:?}", params[i], points[i]);
    }
    println!("  ...");
    for i in n - 3..n {
        println!("  u={:.4}  ->  {:?}", params[i], points[i]);
    }
}
