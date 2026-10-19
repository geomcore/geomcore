//! The 0.36 names still compile in 0.37, with a deprecation warning
//! (ADR 0013): constructors, enum variants, generic bounds and trait
//! methods all go through the old names.
#![allow(deprecated)]

use geomcore::{
    Bound3D, Circle3D, Curve3D, DiscretizationCriteria, Discretize3D, Frame3D, Line3D,
    ParametricCurve3D, Point3D, TrimmedCurve3D, Vector3D,
};

fn start_of<C: ParametricCurve3D>(curve: &C) -> Point3D {
    curve.eval_point(curve.bounds().0)
}

#[test]
fn old_names_still_build_and_evaluate() {
    let circle = Circle3D::new(Point3D::ORIGIN, Vector3D::Z, 2.0).unwrap();
    assert_eq!(circle.eval_point(0.0), Point3D::new(2.0, 0.0, 0.0));
    assert_eq!(Frame3D::WORLD.z_direction(), Vector3D::Z);

    let curve = Curve3D::Circle(circle);
    let Curve3D::Circle(back) = curve else {
        panic!("expected the circle variant");
    };
    assert_eq!(back, circle);

    let segment: TrimmedCurve3D<Line3D> =
        TrimmedCurve3D::new(Line3D::new(Point3D::ORIGIN, Vector3D::X).unwrap(), 0.0, 1.0).unwrap();
    assert_eq!(start_of(&segment), Point3D::ORIGIN);
    assert!(segment.bounding_box().is_some());

    let criteria = DiscretizationCriteria::new(0.1).unwrap();
    let discretization = segment.discretize(criteria);
    assert!(discretization.points.len() >= 2);
}
