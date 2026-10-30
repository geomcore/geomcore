//! Golden replay for `ExtrusionSurface` and `RevolutionSurface`
//! (`fixtures/swept_surfaces.json`), to 1e-7: evaluation, derivatives and
//! point projections against the reference kernel. Each case is also
//! checked against its exact B-spline patch, and the circle-swept cases
//! against the `Cylinder` and `Torus` they sweep. The plane cuts
//! (`fixtures/swept_plane_cuts.json`) are replayed through the router.

mod common;

use std::f64::consts::TAU;

use common::shapes::{
    AxisJson, CurveJson, Projection, SurfaceCase, assert_surface_projections,
    assert_surface_samples, tagged_cases,
};
use common::{assert_close, load, point3, vector3};
use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineSurface, Curve, Cylinder, ExtrusionSurface, IntersectionMethod, Plane, Point,
    RevolutionSurface, Surface, SurfaceSurfaceIntersection, Tolerance, Torus, Vector,
};
use serde::Deserialize;

#[derive(Deserialize)]
struct PlaneJson {
    origin: [f64; 3],
    normal: [f64; 3],
}

#[derive(Deserialize)]
struct ExtrusionCut {
    name: String,
    profile: CurveJson,
    direction: [f64; 3],
    plane: PlaneJson,
    /// Points on the cut; they need not cover it.
    cut: Vec<[f64; 3]>,
}

#[derive(Deserialize)]
struct RevolutionCut {
    name: String,
    profile: CurveJson,
    axis: AxisJson,
    plane: PlaneJson,
    cut: Vec<[f64; 3]>,
}

#[derive(Deserialize)]
struct CutFixture {
    extrusions: Vec<ExtrusionCut>,
    revolutions: Vec<RevolutionCut>,
}

/// The extrusion cases of the fixture, then its revolution cases.
fn fixture() -> (Vec<SurfaceCase>, Vec<SurfaceCase>) {
    let fixture = load("swept_surfaces.json");
    (
        tagged_cases(&fixture, "extrusions", "extrusion"),
        tagged_cases(&fixture, "revolutions", "revolution"),
    )
}

/// The parameter window of `profile` the patches cover: its bounds when
/// finite, otherwise the window the fixture's projections were searched
/// in, `line` for a line.
fn profile_window(profile: &Curve, line: (f64, f64)) -> (f64, f64) {
    match profile {
        Curve::Line(_) => line,
        Curve::Parabola(_) => (-6.0, 6.0),
        Curve::Hyperbola(_) => (-3.0, 3.0),
        _ => profile.bounds(),
    }
}

/// The patch lies on the surface, and projects each fixture query at the
/// fixture's distance.
#[track_caller]
fn assert_patch(
    name: &str,
    surface: &impl ParametricSurface,
    patch: &BSplineSurface,
    projections: &[Projection<[f64; 3]>],
) {
    let ((u0, u1), (v0, v1)) = (patch.u_bounds(), patch.v_bounds());
    for i in 0..=6 {
        for j in 0..=6 {
            let (u, v) = (
                u0 + (u1 - u0) * f64::from(i) / 6.0,
                v0 + (v1 - v0) * f64::from(j) / 6.0,
            );
            let on_patch = patch.eval_point(u, v);
            let scale = (on_patch - Point::ORIGIN).magnitude().max(1.0);
            let distance = surface.project_point(on_patch, Tolerance::DEFAULT).distance;
            assert!(
                distance <= 1e-7 * scale,
                "{name}: patch point at ({u}, {v}) is {distance} off the surface"
            );
        }
    }
    for p in projections {
        let distance = patch
            .project_point(point3(p.point), Tolerance::DEFAULT)
            .distance;
        assert_close(distance, p.distance);
    }
}

fn cut_fixture() -> CutFixture {
    serde_json::from_value(load("swept_plane_cuts.json")).unwrap()
}

fn plane(json: &PlaneJson) -> Plane {
    Plane::new(point3(json.origin), vector3(json.normal)).unwrap()
}

/// Checks the routed cut of `swept` by `plane`: complete, exact unless the
/// profile is a spline, through every reference point, on both surfaces,
/// and the same in either order and through the `Surface` enum.
fn assert_cut<S>(name: &str, swept: &S, plane: &Plane, spline_profile: bool, reference: &[[f64; 3]])
where
    S: ParametricSurface
        + Clone
        + Into<Surface>
        + Intersect<Plane, Output = SurfaceSurfaceIntersection>,
    Plane: Intersect<S, Output = SurfaceSurfaceIntersection>,
{
    let tol = Tolerance::DEFAULT;
    let result = swept.intersect(plane, tol);
    assert!(result.is_complete(), "{name}: incomplete");
    if !spline_profile {
        assert_eq!(result.method(), IntersectionMethod::Exact, "{name}");
    }
    for &point in reference {
        let point = point3(point);
        let nearest = result
            .curves()
            .iter()
            .map(|cut| cut.curve.project_point(point, tol).distance)
            .fold(f64::INFINITY, f64::min);
        assert!(
            nearest < 1e-7,
            "{name}: reference point {point:?} is {nearest} off the cut"
        );
    }
    for cut in result.curves() {
        let (first, last) = cut.curve.bounds();
        let (first, last) = (first.max(-5.0), last.min(5.0));
        for index in 0..=16 {
            let point = cut
                .curve
                .eval_point(first + (last - first) * f64::from(index) / 16.0);
            let off = swept
                .project_point(point, tol)
                .distance
                .max(plane.project_point(point, tol).distance);
            assert!(
                off < 1e-7,
                "{name}: cut point {point:?} is {off} off the surfaces"
            );
        }
    }
    assert_eq!(
        plane.intersect(swept, tol),
        result,
        "{name}: reversed order"
    );
    let (wrapped, wrapped_plane) = (swept.clone().into(), Surface::from(*plane));
    assert_eq!(
        wrapped.intersect(&wrapped_plane, tol),
        result,
        "{name}: through the enum"
    );
}

fn is_spline(profile: &Curve) -> bool {
    matches!(profile, Curve::BSpline(_) | Curve::Bezier(_))
}

#[test]
fn extrusion_plane_cuts_match_golden_fixture() {
    for case in cut_fixture().extrusions {
        let profile = case.profile.decode();
        let spline = is_spline(&profile);
        let surface = ExtrusionSurface::new(profile, vector3(case.direction)).unwrap();
        assert_cut(&case.name, &surface, &plane(&case.plane), spline, &case.cut);
    }
}

#[test]
fn revolution_plane_cuts_match_golden_fixture() {
    for case in cut_fixture().revolutions {
        let axis = Axis::new(point3(case.axis.origin), vector3(case.axis.direction)).unwrap();
        let profile = case.profile.decode();
        let spline = is_spline(&profile);
        let surface = RevolutionSurface::new(profile, axis).unwrap();
        assert_cut(&case.name, &surface, &plane(&case.plane), spline, &case.cut);
    }
}

#[test]
fn extrusions_match_golden_fixture() {
    let (extrusions, _) = fixture();
    assert!(!extrusions.is_empty());
    for case in &extrusions {
        let Surface::Extrusion(surface) = case.surface.decode() else {
            panic!("{}: not an extrusion", case.name)
        };
        assert_surface_samples(&surface, &case.samples);
        assert_surface_projections(&case.name, &surface, &case.projections);
        // Through the enum as well: the dispatch reaches the same code.
        let through_enum = Surface::from(surface.clone());
        assert_surface_samples(&through_enum, &case.samples);
        assert_surface_projections(&case.name, &through_enum, &case.projections);

        let patch = surface
            .to_bspline_patch(
                profile_window(surface.profile(), (-1000.0, 1000.0)),
                (-20.0, 20.0),
            )
            .unwrap_or_else(|e| panic!("{}: {e}", case.name));
        assert_patch(&case.name, &surface, &patch, &case.projections);
    }
}

#[test]
fn revolutions_match_golden_fixture() {
    let (_, revolutions) = fixture();
    assert!(!revolutions.is_empty());
    for case in &revolutions {
        let Surface::Revolution(surface) = case.surface.decode() else {
            panic!("{}: not a revolution", case.name)
        };
        assert_surface_samples(&surface, &case.samples);
        assert_surface_projections(&case.name, &surface, &case.projections);
        let through_enum = Surface::from(surface.clone());
        assert_surface_samples(&through_enum, &case.samples);
        assert_surface_projections(&case.name, &through_enum, &case.projections);

        let patch = surface
            .to_bspline_patch((0.0, TAU), profile_window(surface.profile(), (-10.0, 10.0)))
            .unwrap_or_else(|e| panic!("{}: {e}", case.name));
        assert_patch(&case.name, &surface, &patch, &case.projections);
    }
}

#[test]
fn circle_swept_cases_match_their_elementary_surfaces() {
    let (extrusions, revolutions) = fixture();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();

    let name = "circle along its axis (a cylinder)";
    let case = extrusions.iter().find(|c| c.name == name).unwrap();
    assert_surface_samples(&cylinder, &case.samples);
    assert_surface_projections(name, &cylinder, &case.projections);

    let name = "line parallel to the axis (a cylinder)";
    let case = revolutions.iter().find(|c| c.name == name).unwrap();
    assert_surface_samples(&cylinder, &case.samples);
    assert_surface_projections(name, &cylinder, &case.projections);

    let name = "circle meridian (a torus)";
    let case = revolutions.iter().find(|c| c.name == name).unwrap();
    assert_surface_samples(&torus, &case.samples);
    assert_surface_projections(name, &torus, &case.projections);
}
