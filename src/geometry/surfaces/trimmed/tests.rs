use super::*;
use crate::{Cylinder, Plane, Sphere};
use std::f64::consts::{PI, TAU};

const TOL: Tolerance = Tolerance::DEFAULT;

fn cylinder() -> Cylinder {
    Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap()
}

#[test]
fn test_new_refuses_empty_non_finite_and_outside_ranges() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    assert_eq!(
        TrimmedSurface::new(plane, (0.0, f64::INFINITY), (0.0, 1.0)),
        Err(TrimmedSurfaceConstructionError::NonFinite)
    );
    assert_eq!(
        TrimmedSurface::new(plane, (0.0, 1.0), (2.0, 1.0)),
        Err(TrimmedSurfaceConstructionError::EmptyRange)
    );
    // The sphere's `v` runs from pole to pole only.
    assert_eq!(
        TrimmedSurface::new(sphere, (0.0, 1.0), (-2.0, 0.0)),
        Err(TrimmedSurfaceConstructionError::OutsideBasis)
    );
    assert_eq!(
        TrimmedSurface::new(cylinder(), (0.0, 7.0), (0.0, 1.0)),
        Err(TrimmedSurfaceConstructionError::ExceedsPeriod)
    );
    assert!(TrimmedSurface::new(cylinder(), (-1.0, 1.0), (0.0, 1.0)).is_ok());
    assert!(TrimmedSurface::new(cylinder(), (0.0, TAU), (0.0, 1.0)).is_ok());
}

#[test]
fn test_bounds_are_the_ranges_and_evaluation_is_the_basis() {
    let patch = cylinder()
        .trimmed((1.5 * PI, 2.5 * PI), (0.0, 2.0))
        .unwrap();
    assert_eq!(patch.u_bounds(), (1.5 * PI, 2.5 * PI));
    assert_eq!(patch.v_bounds(), (0.0, 2.0));
    assert_eq!(patch.u_period(), None);
    assert_eq!(patch.point_at(2.0, 1.0), cylinder().point_at(2.0, 1.0));
    assert!(
        patch
            .point_at(2.0 * PI, 1.0)
            .distance(cylinder().point_at(0.0, 1.0))
            <= 1e-12
    );
}

#[test]
fn test_parameters_of_shift_into_a_range_across_the_seam() {
    let patch = cylinder()
        .trimmed((1.5 * PI, 2.5 * PI), (0.0, 2.0))
        .unwrap();
    let point = cylinder().point_at(0.25 * PI, 1.0);
    let (u, v) = patch.parameters_of(point, TOL);
    assert!((u - 2.25 * PI).abs() <= 1e-12);
    assert!((v - 1.0).abs() <= 1e-12);
}

#[test]
fn test_project_point_lands_on_the_boundary_outside_the_ranges() {
    let patch = cylinder().trimmed((0.0, PI), (0.0, 2.0)).unwrap();
    // Above the patch: the closest point is on the top boundary curve.
    let above = patch.project_point(Point::new(0.0, 2.0, 5.0), TOL);
    assert!((above.v - 2.0).abs() <= 1e-9);
    assert!((above.distance - (1.0f64 + 9.0).sqrt()).abs() <= 1e-9);
    // Inside: the basis's projection.
    let inside = patch.project_point(Point::new(0.0, 2.0, 1.0), TOL);
    assert!((inside.u - 0.5 * PI).abs() <= 1e-12);
    assert!((inside.distance - 1.0).abs() <= 1e-12);
}

#[test]
fn test_holds_accepts_a_point_within_confusion_past_an_end() {
    let patch = cylinder().trimmed((0.0, PI), (0.0, 2.0)).unwrap();
    let restriction = patch.restriction();
    assert!(restriction.holds((PI + 0.5 * TOL.confusion, 1.0), TOL));
    assert!(!restriction.holds((PI + 1e-3, 1.0), TOL));
    assert!(restriction.holds((0.5 * PI + TAU, 1.0), TOL));
}

#[test]
fn test_the_bounding_box_holds_the_patch_only() {
    let patch = cylinder().trimmed((0.0, PI), (0.0, 2.0)).unwrap();
    let bounds = patch.bounding_box().unwrap();
    assert!(
        bounds.min_corner().y >= -1e-12,
        "the patch is the half y >= 0"
    );
    for i in 0..=8 {
        for j in 0..=8 {
            let point = patch.point_at(PI * f64::from(i) / 8.0, 2.0 * f64::from(j) / 8.0);
            assert!(bounds.enlarged(1e-12).encloses(point));
        }
    }
}
