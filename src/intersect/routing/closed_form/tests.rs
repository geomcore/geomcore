use super::*;

use crate::curves::{BSplineCurve2D, BezierCurve, BezierCurve2D};
use crate::intersect::routing::{CurveOperand, CurveOperand2D, SurfaceOperand};
use crate::intersect::{Intersect, IntersectionMethod};
use crate::surfaces::{
    BezierSurface, ExtrusionSurface, ParametricSurface, RevolutionSurface, Surface,
};
use crate::{Axis, Point, Point2D, Vector, Vector2D};

const TOL: Tolerance = Tolerance::DEFAULT;

/// One curve of every [`Curve`] variant, placed where every pair the table
/// lists with the samples below is solved in closed form.
fn curves() -> Vec<Curve> {
    let lifted = |x: f64, z: f64| Point::new(x, 0.1, z);
    let bezier = BezierCurve::new(
        2,
        vec![lifted(-3.0, -1.0), lifted(0.0, 3.0), lifted(3.0, -1.0)],
    )
    .unwrap();
    vec![
        Line::new(Point::new(0.3, -0.2, -5.0), Vector::new(0.1, 0.2, 1.0))
            .unwrap()
            .into(),
        Circle::new(Point::new(0.0, 0.0, 0.3), Vector::Z, 2.0)
            .unwrap()
            .into(),
        Ellipse::new(Point::new(0.0, 0.0, -0.4), Vector::Z, Vector::X, 3.0, 1.5)
            .unwrap()
            .into(),
        Parabola::new(Point::new(0.0, 0.0, -3.0), Vector::Y, Vector::Z, 0.5)
            .unwrap()
            .into(),
        Hyperbola::new(Point::new(0.0, 0.0, 0.2), Vector::Y, Vector::X, 1.0, 2.0)
            .unwrap()
            .into(),
        bezier.to_bspline().into(),
        bezier.into(),
    ]
}

/// One surface of every [`Surface`] variant, coaxial about z, so every
/// pair the table lists among them is solved in closed form.
fn surfaces() -> Vec<Surface> {
    let z = Vector::Z;
    let axis = Axis::new(Point::ORIGIN, z).unwrap();
    let profile: Curve = Line::new(Point::new(2.0, 0.0, -1.0), Vector::new(0.5, 0.0, 1.0))
        .unwrap()
        .into();
    let circle: Curve = Circle::new(Point::ORIGIN, z, 1.5).unwrap().into();
    let row = |x: f64, height: f64| {
        vec![
            Point::new(x, -1.0, 0.0),
            Point::new(x, 0.0, height),
            Point::new(x, 1.0, 0.0),
        ]
    };
    let dome =
        BezierSurface::new(2, 2, vec![row(-1.0, 0.0), row(0.0, 2.0), row(1.0, 0.0)]).unwrap();
    vec![
        Plane::new(Point::new(0.0, 0.0, 0.5), z).unwrap().into(),
        Cylinder::new(Point::ORIGIN, z, 1.0).unwrap().into(),
        Cone::new(Point::new(0.0, 0.0, 1.0), z, 0.4, 0.5)
            .unwrap()
            .into(),
        Sphere::new(Point::ORIGIN, 2.0).unwrap().into(),
        Torus::new(Point::ORIGIN, z, 4.0, 1.0).unwrap().into(),
        ExtrusionSurface::new(circle, z).unwrap().into(),
        RevolutionSurface::new(profile, axis).unwrap().into(),
        dome.to_bspline().into(),
        dome.into(),
    ]
}

/// One curve of every [`Curve2D`] variant.
fn curves_2d() -> Vec<Curve2D> {
    let bezier = BezierCurve2D::new(
        2,
        vec![
            Point2D::new(-3.0, -1.0),
            Point2D::new(0.0, 3.0),
            Point2D::new(3.0, -1.0),
        ],
    )
    .unwrap();
    let spline: BSplineCurve2D = bezier.to_bspline();
    vec![
        Line2D::new(Point2D::new(0.3, -5.0), Vector2D::new(0.1, 1.0))
            .unwrap()
            .into(),
        Circle2D::new(Point2D::new(0.2, 0.0), 2.0).unwrap().into(),
        Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5)
            .unwrap()
            .into(),
        Parabola2D::new(Point2D::new(0.0, -3.0), Vector2D::Y, 0.5)
            .unwrap()
            .into(),
        Hyperbola2D::new(Point2D::new(0.0, 0.2), Vector2D::X, 1.0, 2.0)
            .unwrap()
            .into(),
        spline.into(),
        bezier.into(),
    ]
}

// The expected entries, written from the shapes rather than read from the
// table.

fn is_elementary_curve(curve: &Curve) -> bool {
    !matches!(curve, Curve::BSpline(_) | Curve::Bezier(_))
}

fn is_elementary_surface(surface: &Surface) -> bool {
    !matches!(
        surface,
        Surface::Extrusion(_) | Surface::Revolution(_) | Surface::BSpline(_) | Surface::Bezier(_)
    )
}

/// A line meets every elementary surface in closed form, a conic every
/// quadric.
fn curve_surface_listed(curve: &Curve, surface: &Surface) -> bool {
    match (curve, surface) {
        _ if !is_elementary_curve(curve) || !is_elementary_surface(surface) => false,
        (Curve::Line(_), _) => true,
        (_, Surface::Torus(_)) => false,
        _ => true,
    }
}

/// Every pair of elementary surfaces, and a plane with a swept surface.
fn surface_surface_listed(first: &Surface, second: &Surface) -> bool {
    let swept =
        |surface: &Surface| matches!(surface, Surface::Extrusion(_) | Surface::Revolution(_));
    let plane = |surface: &Surface| matches!(surface, Surface::Plane(_));
    (is_elementary_surface(first) && is_elementary_surface(second))
        || (plane(first) && swept(second))
        || (swept(first) && plane(second))
}

fn is_elementary_curve_2d(curve: &Curve2D) -> bool {
    !matches!(curve, Curve2D::BSpline(_) | Curve2D::Bezier(_))
}

#[test]
fn test_curve_surface_pairs_route_to_exact_exactly_when_listed() {
    for curve in &curves() {
        for surface in &surfaces() {
            let listed = curve_surface_listed(curve, surface);
            let (basis, _) = curve.basis_and_range();
            assert_eq!(
                basis.has_closed_form(&surface.surface_basis()),
                listed,
                "{curve:?} and {surface:?}"
            );
            let result = curve.intersect(surface, TOL);
            assert_eq!(
                result.method() == IntersectionMethod::Exact,
                listed,
                "{curve:?} and {surface:?}: {result:?}"
            );
            assert_eq!(surface.intersect(curve, TOL), result);
        }
    }
}

#[test]
fn test_surface_pairs_route_to_exact_exactly_when_listed() {
    for first in &surfaces() {
        for second in &surfaces() {
            let listed = surface_surface_listed(first, second);
            let (first_kind, second_kind) = (first.surface_basis(), second.surface_basis());
            assert_eq!(
                first_kind.has_closed_form(&second_kind),
                listed,
                "{first:?} and {second:?}"
            );
            // The march stands in for nothing here: the router's one call
            // is checked for whether it reaches it.
            let mut marched = false;
            first_kind.exact_family(&second_kind, (None, None), TOL, || {
                marched = true;
                SurfaceSurfaceIntersection::exact(Vec::new(), Vec::new(), false)
            });
            assert_eq!(!marched, listed, "{first:?} and {second:?}");
        }
    }
}

#[test]
fn test_curve_pairs_route_to_exact_exactly_when_listed() {
    for first in &curves() {
        for second in &curves() {
            let listed = is_elementary_curve(first) && is_elementary_curve(second);
            let (first_basis, _) = first.basis_and_range();
            let (second_basis, _) = second.basis_and_range();
            assert_eq!(
                first_basis.has_closed_form(&second_basis),
                listed,
                "{first:?} and {second:?}"
            );
            let result = first.intersect(second, TOL);
            assert_eq!(
                result.method() == IntersectionMethod::Exact,
                listed,
                "{first:?} and {second:?}: {result:?}"
            );
        }
    }
}

#[test]
fn test_2d_curve_pairs_route_to_exact_exactly_when_listed() {
    for first in &curves_2d() {
        for second in &curves_2d() {
            let listed = is_elementary_curve_2d(first) && is_elementary_curve_2d(second);
            let (first_basis, _) = first.basis_and_range();
            let (second_basis, _) = second.basis_and_range();
            assert_eq!(
                first_basis.has_closed_form(&second_basis),
                listed,
                "{first:?} and {second:?}"
            );
            let result = first.intersect(second, TOL);
            assert_eq!(
                result.method() == IntersectionMethod::Exact,
                listed,
                "{first:?} and {second:?}: {result:?}"
            );
        }
    }
}

#[test]
fn test_typed_and_enum_bases_agree_on_the_table() {
    // A concrete shape is its own basis: the typed entry and the one read
    // through the enums give the same answer.
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    for surface in surfaces() {
        let Some(elementary) = surface.as_elementary() else {
            continue;
        };
        let kind = surface.surface_basis();
        assert_eq!(
            circle.has_closed_form(&elementary),
            CurveKind::from(&Curve::from(circle)).has_closed_form(&kind)
        );
        assert!(line.has_closed_form(&elementary));
    }
}
