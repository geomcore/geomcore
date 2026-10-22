//! Constructors reject NaN and infinite inputs, and the conics reject the
//! zero sizes that collapse them onto lines or segments.

use geomcore::{
    Axis, Axis2D, Circle, Circle2D, CircleConstructionError, Cone, ConeConstructionError, Cylinder,
    CylinderConstructionError, Ellipse, Ellipse2D, EllipseConstructionError, Frame, Frame2D,
    FrameConstructionError, Hyperbola, HyperbolaConstructionError, Line, Line2D,
    LineConstructionError, Parabola, Parabola2D, ParabolaConstructionError, Plane,
    PlaneConstructionError, Point, Point2D, Sphere, SphereConstructionError, Torus,
    TorusConstructionError, Vector, Vector2D,
};

const BAD: [f64; 3] = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
const O: Point = Point::ORIGIN;
const X: Vector = Vector::X;
const Z: Vector = Vector::Z;

fn bad_points() -> impl Iterator<Item = Point> {
    BAD.into_iter().map(|b| Point::new(1.0, b, 0.0))
}

fn bad_vectors() -> impl Iterator<Item = Vector> {
    BAD.into_iter().map(|b| Vector::new(b, 0.0, 1.0))
}

#[test]
fn test_frames_and_axes_reject_non_finite() {
    let nf = Some(FrameConstructionError::NonFinite);
    for p in bad_points() {
        assert_eq!(Axis::new(p, Z).err(), nf);
        assert_eq!(Frame::new(p, Z, X).err(), nf);
        assert_eq!(Frame::from_z(p, Z).err(), nf);
    }
    for v in bad_vectors() {
        assert_eq!(Axis::new(O, v).err(), nf);
        assert_eq!(Frame::new(O, v, X).err(), nf);
        assert_eq!(Frame::new(O, Z, v).err(), nf);
        assert_eq!(Frame::from_z(O, v).err(), nf);
    }
    for b in BAD {
        let (p, v) = (Point2D::new(b, 0.0), Vector2D::new(1.0, b));
        assert_eq!(Axis2D::new(p, Vector2D::X).err(), nf);
        assert_eq!(Axis2D::new(Point2D::ORIGIN, v).err(), nf);
        assert_eq!(Frame2D::new(p, Vector2D::X, Vector2D::Y).err(), nf);
        assert_eq!(Frame2D::new(Point2D::ORIGIN, v, Vector2D::Y).err(), nf);
        assert_eq!(Frame2D::from_x(p, Vector2D::X).err(), nf);
        assert_eq!(Frame2D::from_x(Point2D::ORIGIN, v).err(), nf);
    }
}

#[test]
fn test_lines_and_planes_reject_non_finite() {
    for p in bad_points() {
        assert_eq!(Line::new(p, X), Err(LineConstructionError::NonFinite));
        assert_eq!(
            Line::from_two_points(O, p),
            Err(LineConstructionError::NonFinite)
        );
        assert_eq!(Plane::new(p, Z), Err(PlaneConstructionError::NonFinite));
        assert_eq!(
            Plane::from_three_points(O, Point::new(1.0, 0.0, 0.0), p),
            Err(PlaneConstructionError::NonFinite)
        );
    }
    for v in bad_vectors() {
        assert_eq!(Line::new(O, v), Err(LineConstructionError::NonFinite));
        // Used to panic while picking the frame's x direction.
        assert_eq!(Plane::new(O, v), Err(PlaneConstructionError::NonFinite));
    }
    for b in BAD {
        let p = Point2D::new(b, 0.0);
        assert_eq!(
            Line2D::new(p, Vector2D::X),
            Err(LineConstructionError::NonFinite)
        );
        assert_eq!(
            Line2D::from_two_points(Point2D::ORIGIN, p),
            Err(LineConstructionError::NonFinite)
        );
        assert_eq!(
            Plane::from_coefficients(b, 0.0, 1.0, 0.0),
            Err(PlaneConstructionError::NonFinite)
        );
        assert_eq!(
            Plane::from_coefficients(0.0, 0.0, 1.0, b),
            Err(PlaneConstructionError::NonFinite)
        );
    }
}

#[test]
fn test_circles_reject_non_finite() {
    let nf = Some(CircleConstructionError::NonFinite);
    for b in BAD {
        assert_eq!(Circle::new(O, Z, b).err(), nf);
        assert_eq!(Circle::from_frame(Frame::WORLD, b).err(), nf);
        assert_eq!(Circle2D::new(Point2D::ORIGIN, b).err(), nf);
        assert_eq!(Circle2D::new(Point2D::new(b, 0.0), 1.0).err(), nf);
    }
    for p in bad_points() {
        assert_eq!(Circle::new(p, Z, 1.0).err(), nf);
        assert_eq!(
            Circle::from_three_points(O, Point::new(1.0, 0.0, 0.0), p).err(),
            nf
        );
    }
    // A zero radius is a valid (point) circle.
    assert!(Circle::new(O, Z, 0.0).is_ok());
}

#[test]
fn test_conics_reject_non_finite_and_zero_sizes() {
    for b in BAD {
        assert_eq!(
            Ellipse::new(O, Z, X, b, 1.0),
            Err(EllipseConstructionError::NonFinite)
        );
        assert_eq!(
            Ellipse::new(O, Z, X, 3.0, b),
            Err(EllipseConstructionError::NonFinite)
        );
        assert_eq!(
            Hyperbola::new(O, Z, X, b, 1.0),
            Err(HyperbolaConstructionError::NonFinite)
        );
        assert_eq!(
            Hyperbola::new(O, Z, X, 2.0, b),
            Err(HyperbolaConstructionError::NonFinite)
        );
        assert_eq!(
            Parabola::new(O, Z, X, b),
            Err(ParabolaConstructionError::NonFinite)
        );
    }
    for p in bad_points() {
        assert_eq!(
            Ellipse::new(p, Z, X, 3.0, 1.0),
            Err(EllipseConstructionError::NonFinite)
        );
        assert_eq!(
            Ellipse::from_center_and_points(O, Point::new(3.0, 0.0, 0.0), p),
            Err(EllipseConstructionError::NonFinite)
        );
        assert_eq!(
            Hyperbola::from_center_and_points(O, Point::new(2.0, 0.0, 0.0), p),
            Err(HyperbolaConstructionError::NonFinite)
        );
    }
    assert_eq!(
        Ellipse::new(O, Z, X, 3.0, 0.0),
        Err(EllipseConstructionError::NonPositiveRadius)
    );
    assert_eq!(
        Hyperbola::new(O, Z, X, 0.0, 1.0),
        Err(HyperbolaConstructionError::NonPositiveRadius)
    );
    assert_eq!(
        Parabola::new(O, Z, X, 0.0),
        Err(ParabolaConstructionError::NonPositiveFocal)
    );

    for b in BAD {
        assert_eq!(
            Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, b, 1.0),
            Err(EllipseConstructionError::NonFinite)
        );
        assert_eq!(
            Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, b),
            Err(EllipseConstructionError::NonFinite)
        );
        assert_eq!(
            Ellipse2D::new(Point2D::new(b, 0.0), Vector2D::X, 3.0, 1.0),
            Err(EllipseConstructionError::NonFinite)
        );
        assert_eq!(
            Ellipse2D::new(Point2D::ORIGIN, Vector2D::new(1.0, b), 3.0, 1.0),
            Err(EllipseConstructionError::NonFinite)
        );
    }
    assert_eq!(
        Ellipse2D::new(Point2D::ORIGIN, Vector2D::new(0.0, 0.0), 3.0, 1.0),
        Err(EllipseConstructionError::NullNormal)
    );
    assert_eq!(
        Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 0.0),
        Err(EllipseConstructionError::NonPositiveRadius)
    );
    assert_eq!(
        Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 1.0, 3.0),
        Err(EllipseConstructionError::InvertedRadii)
    );
    assert_eq!(
        Ellipse2D::from_center_and_points(
            Point2D::ORIGIN,
            Point2D::new(3.0, 0.0),
            Point2D::new(6.0, 0.0)
        ),
        Err(EllipseConstructionError::InvertedAxis)
    );
    assert_eq!(
        Ellipse2D::from_center_and_points(Point2D::ORIGIN, Point2D::ORIGIN, Point2D::new(0.0, 1.0)),
        Err(EllipseConstructionError::NullAxis)
    );

    for b in BAD {
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::X, b),
            Err(ParabolaConstructionError::NonFinite)
        );
        assert_eq!(
            Parabola2D::new(Point2D::new(b, 0.0), Vector2D::X, 1.0),
            Err(ParabolaConstructionError::NonFinite)
        );
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::new(1.0, b), 1.0),
            Err(ParabolaConstructionError::NonFinite)
        );
        assert_eq!(
            Parabola2D::from_frame(Frame2D::WORLD, b),
            Err(ParabolaConstructionError::NonFinite)
        );
    }
    assert_eq!(
        Parabola2D::new(Point2D::ORIGIN, Vector2D::new(0.0, 0.0), 1.0),
        Err(ParabolaConstructionError::NullNormal)
    );
    for focal in [0.0, -1.0] {
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::X, focal),
            Err(ParabolaConstructionError::NonPositiveFocal)
        );
    }
}

#[test]
fn test_elementary_surfaces_reject_non_finite() {
    for b in BAD {
        assert_eq!(
            Cylinder::new(O, Z, b),
            Err(CylinderConstructionError::NonFinite)
        );
        assert_eq!(Sphere::new(O, b), Err(SphereConstructionError::NonFinite));
        assert_eq!(
            Torus::new(O, Z, b, 1.0),
            Err(TorusConstructionError::NonFinite)
        );
        assert_eq!(
            Torus::new(O, Z, 4.0, b),
            Err(TorusConstructionError::NonFinite)
        );
        assert_eq!(
            Cone::new(O, Z, b, 1.0),
            Err(ConeConstructionError::NonFinite)
        );
        assert_eq!(
            Cone::new(O, Z, 0.4, b),
            Err(ConeConstructionError::NonFinite)
        );
        assert_eq!(
            Cone::from_two_points_and_radii(O, Point::new(0.0, 0.0, 1.0), 1.0, b),
            Err(ConeConstructionError::NonFinite)
        );
    }
    for p in bad_points() {
        assert_eq!(
            Cylinder::new(p, Z, 1.0),
            Err(CylinderConstructionError::NonFinite)
        );
        assert_eq!(Sphere::new(p, 1.0), Err(SphereConstructionError::NonFinite));
        assert_eq!(
            Torus::new(p, Z, 4.0, 1.0),
            Err(TorusConstructionError::NonFinite)
        );
        assert_eq!(
            Cone::new(p, Z, 0.4, 1.0),
            Err(ConeConstructionError::NonFinite)
        );
        assert_eq!(
            Cone::from_two_points_and_radii(O, p, 1.0, 2.0),
            Err(ConeConstructionError::NonFinite)
        );
    }
    // Zero radii stay valid: point sphere, line cylinder, circle torus,
    // cone with its apex at the origin.
    assert!(Sphere::new(O, 0.0).is_ok());
    assert!(Cylinder::new(O, Z, 0.0).is_ok());
    assert!(Torus::new(O, Z, 4.0, 0.0).is_ok());
    assert!(Cone::new(O, Z, 0.4, 0.0).is_ok());
}
