use super::*;

const TOL: Tolerance = Tolerance::DEFAULT;

#[test]
fn test_trimmed_operand_reports_basis_and_range() {
    let circle = crate::Circle::new(crate::Point::ORIGIN, crate::Vector::Z, 1.0).unwrap();
    let arc = crate::geometry::curves::TrimmedCurve::new(circle, -1.0, 1.0).unwrap();
    let (basis, range) = arc.basis_and_range();
    // A trimmed circle routes through the typed circle solver.
    assert_eq!(basis, circle);
    assert_eq!(range, Some((-1.0, 1.0)));
}

#[test]
fn test_trimmed_of_trimmed_keeps_the_outer_range() {
    let line = crate::Line::new(crate::Point::ORIGIN, crate::Vector::X).unwrap();
    let outer = TrimmedCurve::new(TrimmedCurve::new(line, 0.0, 10.0).unwrap(), 2.0, 3.0).unwrap();
    assert_eq!(outer.basis_and_range().1, Some((2.0, 3.0)));
    let spline_free: Curve = line.into();
    assert_eq!(spline_free.basis_and_range().1, None);
}

#[test]
fn test_exact_surface_curves_carry_their_closed_form_curves_in_parameters() {
    use crate::{Curve2D, Cylinder, ParametricCurve, ParametricCurve2D, Plane, Vector};
    let tol = Tolerance::DEFAULT;
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    // An oblique plane cuts an ellipse: an ellipse in the plane's
    // parameters, and no closed form in the cylinder's.
    let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0)).unwrap();
    let result = plane.intersect(&cylinder, tol);
    let piece = &result.curves()[0];
    assert!(matches!(piece.first_uv_curve, Some(Curve2D::Ellipse(_))));
    assert_eq!(piece.second_uv_curve, None);
    // In the other order, the other way round.
    let result = cylinder.intersect(&plane, tol);
    let piece = &result.curves()[0];
    assert_eq!(piece.first_uv_curve, None);
    assert!(matches!(piece.second_uv_curve, Some(Curve2D::Ellipse(_))));
    // A plane holding the axis cuts two rulings: a line in both.
    let plane = Plane::new(Point::ORIGIN, Vector::Y).unwrap();
    let result = cylinder.intersect(&plane, tol);
    assert_eq!(result.curves().len(), 2);
    for piece in result.curves() {
        let (Some(Curve2D::Line(on_cylinder)), Some(Curve2D::Line(_))) =
            (&piece.first_uv_curve, &piece.second_uv_curve)
        else {
            panic!("a ruling is a line in both parameters: {piece:?}");
        };
        for t in [-1.0, 0.0, 2.5] {
            let uv = ParametricCurve2D::point_at(on_cylinder, t);
            let on_surface = cylinder.point_at(uv.x, uv.y);
            let expected = ParametricCurve::point_at(&piece.curve, t);
            assert!(on_surface.distance(expected) < 1e-12);
        }
    }
}

#[test]
fn test_the_tee_of_two_equal_cylinders_is_two_exact_ellipses() {
    use crate::{Curve, Cylinder, IntersectionMethod, Vector};
    let along_x = Cylinder::new(Point::ORIGIN, Vector::X, 5.0).unwrap();
    let along_z = Cylinder::new(Point::ORIGIN, Vector::Z, 5.0).unwrap();
    let result = along_x.intersect(&along_z, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    assert_eq!(result.curves().len(), 2);
    for piece in result.curves() {
        assert!(matches!(piece.curve, Curve::Ellipse(_)) && !piece.tangent);
    }
    // Both points where the cylinders touch, where the ellipses cross.
    let touching = [Point::new(0.0, -5.0, 0.0), Point::new(0.0, 5.0, 0.0)];
    assert_eq!(result.contact_points().len(), 2);
    for (found, expected) in result.contact_points().iter().zip(touching) {
        assert!(found.distance(expected) < 1e-12, "{found:?}");
    }
}
