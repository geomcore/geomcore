//! Golden replay for `ExtrusionSurface` and `RevolutionSurface`
//! (`fixtures/swept_surfaces.json`), to 1e-7: evaluation, derivatives and
//! point projections against the reference kernel. Each case is also
//! checked against its exact B-spline patch, and the circle-swept cases
//! against the `Cylinder` and `Torus` they sweep. The plane cuts
//! (`fixtures/swept_plane_cuts.json`) are replayed through the router.

mod common;

use std::f64::consts::TAU;

use common::{
    FrameJson, assert_close, assert_point3, assert_vector3, frame3, load, point3, vector3,
};
use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve, BSplineSurface, BezierCurve, Circle, Curve, Cylinder, Ellipse,
    ExtrusionSurface, Hyperbola, IntersectionMethod, Line, Parabola, Plane, Point,
    RevolutionSurface, Surface, SurfaceSurfaceIntersection, Tolerance, Torus, Vector,
};
use serde::Deserialize;

#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
enum ProfileJson {
    Line {
        origin: [f64; 3],
        direction: [f64; 3],
    },
    Circle {
        frame: FrameJson,
        radius: f64,
    },
    Ellipse {
        frame: FrameJson,
        major_radius: f64,
        minor_radius: f64,
    },
    Parabola {
        frame: FrameJson,
        focal: f64,
    },
    Hyperbola {
        frame: FrameJson,
        major_radius: f64,
        minor_radius: f64,
    },
    Bspline {
        degree: usize,
        poles: Vec<[f64; 3]>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
    },
    Bezier {
        degree: usize,
        poles: Vec<[f64; 3]>,
        weights: Option<Vec<f64>>,
    },
}

#[derive(Deserialize)]
struct Sample {
    u: f64,
    v: f64,
    point: [f64; 3],
    d1u: [f64; 3],
    d1v: [f64; 3],
    d2u: [f64; 3],
    d2v: [f64; 3],
    d2uv: [f64; 3],
}

#[derive(Deserialize)]
struct Projection {
    point: [f64; 3],
    /// `None` for a query on the axis of a revolution, whose closest
    /// points fill a circle.
    closest: Option<[f64; 3]>,
    distance: f64,
}

#[derive(Deserialize)]
struct ExtrusionCase {
    name: String,
    profile: ProfileJson,
    direction: [f64; 3],
    samples: Vec<Sample>,
    projections: Vec<Projection>,
}

#[derive(Deserialize)]
struct AxisJson {
    origin: [f64; 3],
    direction: [f64; 3],
}

#[derive(Deserialize)]
struct RevolutionCase {
    name: String,
    profile: ProfileJson,
    axis: AxisJson,
    samples: Vec<Sample>,
    projections: Vec<Projection>,
}

#[derive(Deserialize)]
struct Fixture {
    extrusions: Vec<ExtrusionCase>,
    revolutions: Vec<RevolutionCase>,
}

#[derive(Deserialize)]
struct PlaneJson {
    origin: [f64; 3],
    normal: [f64; 3],
}

#[derive(Deserialize)]
struct ExtrusionCut {
    name: String,
    profile: ProfileJson,
    direction: [f64; 3],
    plane: PlaneJson,
    /// Points on the cut; they need not cover it.
    cut: Vec<[f64; 3]>,
}

#[derive(Deserialize)]
struct RevolutionCut {
    name: String,
    profile: ProfileJson,
    axis: AxisJson,
    plane: PlaneJson,
    cut: Vec<[f64; 3]>,
}

#[derive(Deserialize)]
struct CutFixture {
    extrusions: Vec<ExtrusionCut>,
    revolutions: Vec<RevolutionCut>,
}

fn fixture() -> Fixture {
    serde_json::from_value(load("swept_surfaces.json")).unwrap()
}

fn points(poles: &[[f64; 3]]) -> Vec<Point> {
    poles.iter().map(|&p| point3(p)).collect()
}

fn profile(json: &ProfileJson) -> Curve {
    match json {
        ProfileJson::Line { origin, direction } => Line::new(point3(*origin), vector3(*direction))
            .unwrap()
            .into(),
        ProfileJson::Circle { frame, radius } => {
            Circle::from_frame(frame3(frame), *radius).unwrap().into()
        }
        ProfileJson::Ellipse {
            frame,
            major_radius,
            minor_radius,
        } => Ellipse::from_frame(frame3(frame), *major_radius, *minor_radius)
            .unwrap()
            .into(),
        ProfileJson::Parabola { frame, focal } => {
            Parabola::from_frame(frame3(frame), *focal).unwrap().into()
        }
        ProfileJson::Hyperbola {
            frame,
            major_radius,
            minor_radius,
        } => Hyperbola::from_frame(frame3(frame), *major_radius, *minor_radius)
            .unwrap()
            .into(),
        ProfileJson::Bspline {
            degree,
            poles,
            weights,
            knots,
            multiplicities,
        } => {
            let (poles, knots, multiplicities) =
                (points(poles), knots.clone(), multiplicities.clone());
            match weights {
                Some(weights) => BSplineCurve::new_rational(
                    *degree,
                    poles,
                    weights.clone(),
                    knots,
                    multiplicities,
                    false,
                ),
                None => BSplineCurve::new(*degree, poles, knots, multiplicities, false),
            }
            .unwrap()
            .into()
        }
        ProfileJson::Bezier {
            degree,
            poles,
            weights,
        } => match weights {
            Some(weights) => BezierCurve::new_rational(*degree, points(poles), weights.clone()),
            None => BezierCurve::new(*degree, points(poles)),
        }
        .unwrap()
        .into(),
    }
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

#[track_caller]
fn assert_samples(surface: &impl ParametricSurface, samples: &[Sample]) {
    for s in samples {
        assert_point3(surface.eval_point(s.u, s.v), s.point);
        assert_vector3(surface.eval_derivative(s.u, s.v, 1, 0), s.d1u);
        assert_vector3(surface.eval_derivative(s.u, s.v, 0, 1), s.d1v);
        assert_vector3(surface.eval_derivative(s.u, s.v, 2, 0), s.d2u);
        assert_vector3(surface.eval_derivative(s.u, s.v, 0, 2), s.d2v);
        assert_vector3(surface.eval_derivative(s.u, s.v, 1, 1), s.d2uv);
    }
}

#[track_caller]
fn assert_projections(name: &str, surface: &impl ParametricSurface, projections: &[Projection]) {
    for p in projections {
        let found = surface.project_point(point3(p.point), Tolerance::DEFAULT);
        assert!(
            (found.distance - p.distance).abs() <= 1e-7 * p.distance.max(1.0),
            "{name}: projection of {:?} at distance {} != {}",
            p.point,
            found.distance,
            p.distance
        );
        assert_point3(
            surface.eval_point(found.u, found.v),
            [found.point.x, found.point.y, found.point.z],
        );
        if let Some(closest) = p.closest {
            assert_point3(found.point, closest);
        }
    }
}

/// The patch lies on the surface, and projects each fixture query at the
/// fixture's distance.
#[track_caller]
fn assert_patch(
    name: &str,
    surface: &impl ParametricSurface,
    patch: &BSplineSurface,
    projections: &[Projection],
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

fn is_spline(profile: &ProfileJson) -> bool {
    matches!(
        profile,
        ProfileJson::Bspline { .. } | ProfileJson::Bezier { .. }
    )
}

#[test]
fn extrusion_plane_cuts_match_golden_fixture() {
    for case in cut_fixture().extrusions {
        let surface =
            ExtrusionSurface::new(profile(&case.profile), vector3(case.direction)).unwrap();
        let spline = is_spline(&case.profile);
        assert_cut(&case.name, &surface, &plane(&case.plane), spline, &case.cut);
    }
}

#[test]
fn revolution_plane_cuts_match_golden_fixture() {
    for case in cut_fixture().revolutions {
        let axis = Axis::new(point3(case.axis.origin), vector3(case.axis.direction)).unwrap();
        let surface = RevolutionSurface::new(profile(&case.profile), axis).unwrap();
        let spline = is_spline(&case.profile);
        assert_cut(&case.name, &surface, &plane(&case.plane), spline, &case.cut);
    }
}

#[test]
fn extrusions_match_golden_fixture() {
    for case in fixture().extrusions {
        let surface = ExtrusionSurface::new(profile(&case.profile), vector3(case.direction))
            .unwrap_or_else(|e| panic!("{}: {e}", case.name));
        assert_samples(&surface, &case.samples);
        assert_projections(&case.name, &surface, &case.projections);
        // Through the enum as well: the dispatch reaches the same code.
        let through_enum = Surface::from(surface.clone());
        assert_samples(&through_enum, &case.samples);
        assert_projections(&case.name, &through_enum, &case.projections);

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
    for case in fixture().revolutions {
        let axis = Axis::new(point3(case.axis.origin), vector3(case.axis.direction)).unwrap();
        let surface = RevolutionSurface::new(profile(&case.profile), axis)
            .unwrap_or_else(|e| panic!("{}: {e}", case.name));
        assert_samples(&surface, &case.samples);
        assert_projections(&case.name, &surface, &case.projections);
        let through_enum = Surface::from(surface.clone());
        assert_samples(&through_enum, &case.samples);
        assert_projections(&case.name, &through_enum, &case.projections);

        let patch = surface
            .to_bspline_patch((0.0, TAU), profile_window(surface.profile(), (-10.0, 10.0)))
            .unwrap_or_else(|e| panic!("{}: {e}", case.name));
        assert_patch(&case.name, &surface, &patch, &case.projections);
    }
}

#[test]
fn circle_swept_cases_match_their_elementary_surfaces() {
    let fixture = fixture();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();

    let name = "circle along its axis (a cylinder)";
    let case = fixture.extrusions.iter().find(|c| c.name == name).unwrap();
    assert_samples(&cylinder, &case.samples);
    assert_projections(name, &cylinder, &case.projections);

    let name = "line parallel to the axis (a cylinder)";
    let case = fixture.revolutions.iter().find(|c| c.name == name).unwrap();
    assert_samples(&cylinder, &case.samples);
    assert_projections(name, &cylinder, &case.projections);

    let name = "circle meridian (a torus)";
    let case = fixture.revolutions.iter().find(|c| c.name == name).unwrap();
    assert_samples(&torus, &case.samples);
    assert_projections(name, &torus, &case.projections);
}
