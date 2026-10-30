//! Golden-fixture replay tests: build geomcore values from recorded JSON
//! specs and check parity against recorded expected outputs. Grows with
//! every fixture-replay task.

mod common;

use common::shapes::{
    Curve2DJson, CurveCase, CurveJson, CurveSample, SurfaceCase, SurfaceJson,
    assert_curve_2d_parameters, assert_curve_2d_projections, assert_curve_2d_samples,
    assert_curve_parameters, assert_curve_samples, assert_surface_parameters,
    assert_surface_samples, kind_tagged, tagged_cases,
};
use common::{
    FrameJson, assert_close, assert_point2, assert_point3, assert_vector3, frame3, load, point3,
    vector3,
};
use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve2D, Circle, Cone, Curve, Curve2D, Cylinder, Ellipse, Hyperbola,
    ParametrizeError, Plane, Surface, Tolerance, Transform,
};
use serde::Deserialize;

#[derive(Deserialize)]
struct MirrorAxisSpec {
    origin: [f64; 3],
    dir: [f64; 3],
}

#[derive(Deserialize)]
struct MirrorPlaneSpec {
    frame: FrameJson,
}

#[derive(Deserialize)]
struct TransformCase {
    kind: String,
    spec: serde_json::Value,
    apply_to_points: Vec<[f64; 3]>,
    expected_points: Vec<[f64; 3]>,
    apply_to_vectors: Vec<[f64; 3]>,
    expected_vectors: Vec<[f64; 3]>,
}

#[derive(Deserialize)]
struct TransformsFixture {
    cases: Vec<TransformCase>,
}

fn build_transform(cases: &[TransformCase], index: usize) -> Transform {
    let case = &cases[index];
    match case.kind.as_str() {
        "translation" => {
            #[derive(Deserialize)]
            struct Spec {
                offset: [f64; 3],
            }
            let spec: Spec = serde_json::from_value(case.spec.clone()).unwrap();
            Transform::translation(vector3(spec.offset))
        }
        "rotation" => {
            #[derive(Deserialize)]
            struct Spec {
                axis_origin: [f64; 3],
                axis_dir: [f64; 3],
                angle: f64,
            }
            let spec: Spec = serde_json::from_value(case.spec.clone()).unwrap();
            let axis = Axis::new(point3(spec.axis_origin), vector3(spec.axis_dir)).unwrap();
            Transform::rotation(axis, spec.angle)
        }
        "scale" => {
            #[derive(Deserialize)]
            struct Spec {
                center: [f64; 3],
                factor: f64,
            }
            let spec: Spec = serde_json::from_value(case.spec.clone()).unwrap();
            Transform::scaling(point3(spec.center), spec.factor)
        }
        "mirror_point" => {
            #[derive(Deserialize)]
            struct Spec {
                point: [f64; 3],
            }
            let spec: Spec = serde_json::from_value(case.spec.clone()).unwrap();
            Transform::mirror_point(point3(spec.point))
        }
        "mirror_axis" => {
            let spec: MirrorAxisSpec = serde_json::from_value(case.spec.clone()).unwrap();
            let axis = Axis::new(point3(spec.origin), vector3(spec.dir)).unwrap();
            Transform::mirror_axis(axis)
        }
        "mirror_plane" => {
            let spec: MirrorPlaneSpec = serde_json::from_value(case.spec.clone()).unwrap();
            Transform::mirror_plane(frame3(&spec.frame))
        }
        "compose" => {
            #[derive(Deserialize)]
            struct Spec {
                compose: [usize; 2],
            }
            let spec: Spec = serde_json::from_value(case.spec.clone()).unwrap();
            let first = build_transform(cases, spec.compose[0]);
            let second = build_transform(cases, spec.compose[1]);
            first.then(second)
        }
        other => panic!("unknown transform kind: {other}"),
    }
}

#[test]
fn transforms_match_golden_fixture() {
    let fixture: TransformsFixture = serde_json::from_value(load("transforms.json")).unwrap();

    for (i, case) in fixture.cases.iter().enumerate() {
        let transform = build_transform(&fixture.cases, i);

        assert_eq!(
            case.apply_to_points.len(),
            case.expected_points.len(),
            "case {i} ({}): points/expected_points length mismatch",
            case.kind
        );
        for (p, expected) in case.apply_to_points.iter().zip(&case.expected_points) {
            let actual = transform.apply_point(point3(*p));
            assert_point3(actual, *expected);
        }

        assert_eq!(
            case.apply_to_vectors.len(),
            case.expected_vectors.len(),
            "case {i} ({}): vectors/expected_vectors length mismatch",
            case.kind
        );
        for (v, expected) in case.apply_to_vectors.iter().zip(&case.expected_vectors) {
            let actual = transform.apply_vector(vector3(*v));
            assert_vector3(actual, *expected);
        }
    }
}

#[derive(Deserialize)]
struct ThreePointCase {
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    frame: FrameJson,
    radius: f64,
}

#[derive(Deserialize)]
struct CurvesAnalyticFixture {
    three_point_cases: Vec<ThreePointCase>,
}

/// Replays the `kind` curves listed under `key` in
/// `curves_analytic.json`: their samples, third derivatives included, and
/// inversions.
fn replay_curves(key: &str, kind: &str) {
    let cases: Vec<CurveCase<CurveJson, [f64; 3]>> =
        tagged_cases(&load("curves_analytic.json"), key, kind);
    assert!(!cases.is_empty(), "no {key}");
    for case in &cases {
        assert!(!case.samples.is_empty() && !case.parameter_of.is_empty());
        assert!(case.samples.iter().all(|sample| sample.d3.is_some()));
        let curve = case.curve.decode();
        assert_curve_samples(&curve, &case.samples);
        assert_curve_parameters(&curve, &case.parameter_of);
    }
}

/// Replays the `kind` curves listed under `key` in
/// `curves_analytic_2d.json`: their samples, third derivatives included,
/// inversions and projections, closest parameters and points included.
/// Returns the curves.
fn replay_curves_2d(key: &str, kind: &str) -> Vec<Curve2D> {
    let cases: Vec<CurveCase<Curve2DJson, [f64; 2]>> =
        tagged_cases(&load("curves_analytic_2d.json"), key, kind);
    assert!(!cases.is_empty(), "no {key}");
    cases
        .iter()
        .map(|case| {
            assert!(!case.samples.is_empty() && !case.projections.is_empty());
            assert!(case.samples.iter().all(|sample| sample.d3.is_some()));
            assert!(
                case.projections
                    .iter()
                    .all(|projection| { projection.u.is_some() && projection.closest.is_some() })
            );
            let curve = case.curve.decode();
            assert_curve_2d_samples(&curve, &case.samples);
            assert_curve_2d_parameters(&curve, &case.parameter_of);
            assert_curve_2d_projections(&curve, &case.projections);
            curve
        })
        .collect()
}

#[test]
fn lines_match_golden_fixture() {
    replay_curves("lines", "line");
}

#[test]
fn circles_match_golden_fixture() {
    replay_curves("circles", "circle");
}

#[test]
fn circle_three_point_cases_match_golden_fixture() {
    let fixture: CurvesAnalyticFixture =
        serde_json::from_value(load("curves_analytic.json")).unwrap();

    for (i, case) in fixture.three_point_cases.iter().enumerate() {
        let circle = Circle::from_three_points(point3(case.p1), point3(case.p2), point3(case.p3))
            .unwrap_or_else(|e| panic!("case {i}: failed to build Circle: {e}"));

        assert_close(circle.radius(), case.radius);

        let frame = circle.frame();
        assert_point3(frame.origin(), case.frame.origin);
        assert_vector3(frame.x_direction(), case.frame.x_dir);
        assert_vector3(frame.y_direction(), case.frame.y_dir);
        assert_vector3(frame.z_direction(), case.frame.z_dir);
    }
}

#[test]
fn ellipses_match_golden_fixture() {
    replay_curves("ellipses", "ellipse");
}

#[test]
fn ellipses_2d_match_golden_fixture() {
    let ellipses = replay_curves_2d("ellipses", "ellipse");
    let fixture = load("curves_analytic_2d.json");
    let areas = fixture["ellipses"].as_array().unwrap();
    assert_eq!(ellipses.len(), areas.len());
    for (curve, case) in ellipses.iter().zip(areas) {
        let Curve2D::Ellipse(ellipse) = curve else {
            panic!("not an ellipse: {curve:?}");
        };
        assert_close(ellipse.area(), case["area"].as_f64().unwrap());
    }
}

#[test]
fn parabolas_2d_match_golden_fixture() {
    replay_curves_2d("parabolas", "parabola");
}

#[test]
fn hyperbolas_2d_match_golden_fixture() {
    replay_curves_2d("hyperbolas", "hyperbola");
}

#[test]
fn parabolas_match_golden_fixture() {
    replay_curves("parabolas", "parabola");
}

#[test]
fn hyperbolas_match_golden_fixture() {
    replay_curves("hyperbolas", "hyperbola");
}

#[derive(Deserialize)]
struct EllipseCenterTwoPointsCase {
    s1: [f64; 3],
    s2: [f64; 3],
    center: [f64; 3],
    frame: FrameJson,
    major_radius: f64,
    minor_radius: f64,
}

#[derive(Deserialize)]
struct HyperbolaCenterTwoPointsCase {
    s1: [f64; 3],
    s2: [f64; 3],
    center: [f64; 3],
    frame: FrameJson,
    major_radius: f64,
    minor_radius: f64,
}

#[derive(Deserialize)]
struct ConstructionFixture {
    ellipses_center_two_points: Vec<EllipseCenterTwoPointsCase>,
    hyperbolas_center_two_points: Vec<HyperbolaCenterTwoPointsCase>,
}

#[test]
fn ellipses_center_two_points_match_golden_fixture() {
    let fixture: ConstructionFixture = serde_json::from_value(load("construction.json")).unwrap();

    for (i, case) in fixture.ellipses_center_two_points.iter().enumerate() {
        let ellipse =
            Ellipse::from_center_and_points(point3(case.center), point3(case.s1), point3(case.s2))
                .unwrap_or_else(|e| panic!("case {i}: failed to build Ellipse: {e}"));

        assert_close(ellipse.major_radius(), case.major_radius);
        assert_close(ellipse.minor_radius(), case.minor_radius);

        let frame = ellipse.frame();
        assert_point3(frame.origin(), case.frame.origin);
        assert_vector3(frame.x_direction(), case.frame.x_dir);
        assert_vector3(frame.y_direction(), case.frame.y_dir);
        assert_vector3(frame.z_direction(), case.frame.z_dir);
    }
}

#[test]
fn hyperbolas_center_two_points_match_golden_fixture() {
    let fixture: ConstructionFixture = serde_json::from_value(load("construction.json")).unwrap();

    for (i, case) in fixture.hyperbolas_center_two_points.iter().enumerate() {
        let hyperbola = Hyperbola::from_center_and_points(
            point3(case.center),
            point3(case.s1),
            point3(case.s2),
        )
        .unwrap_or_else(|e| panic!("case {i}: failed to build Hyperbola: {e}"));

        assert_close(hyperbola.major_radius(), case.major_radius);
        assert_close(hyperbola.minor_radius(), case.minor_radius);

        let frame = hyperbola.frame();
        assert_point3(frame.origin(), case.frame.origin);
        assert_vector3(frame.x_direction(), case.frame.x_dir);
        assert_vector3(frame.y_direction(), case.frame.y_dir);
        assert_vector3(frame.z_direction(), case.frame.z_dir);
    }
}

fn bspline_cases() -> Vec<CurveCase<CurveJson, [f64; 3]>> {
    tagged_cases(&load("curves_bspline.json"), "cases", "bspline")
}

#[test]
fn bsplines_match_golden_fixture() {
    for case in &bspline_cases() {
        assert!(!case.samples.is_empty(), "{}", case.name);
        assert_curve_samples(&case.curve.decode(), &case.samples);
    }
}

/// The planar shadow of a recorded B-spline: its poles with `z` dropped.
/// Dropping a coordinate commutes with B-spline evaluation (each
/// coordinate is the same weighted blend of its pole coordinates), so the
/// case's samples with `z` dropped are the shadow's goldens.
fn shadow(case: &CurveCase<CurveJson, [f64; 3]>) -> BSplineCurve2D {
    let CurveJson::Bspline {
        degree,
        poles,
        weights,
        knots,
        multiplicities,
        periodic,
    } = &case.curve
    else {
        panic!("{}: not a B-spline", case.name);
    };
    let poles: Vec<geomcore::Point2D> = poles
        .iter()
        .map(|&[x, y, _]| geomcore::Point2D::new(x, y))
        .collect();
    let (knots, multiplicities) = (knots.clone(), multiplicities.clone());
    let result = match weights {
        Some(weights) => BSplineCurve2D::new_rational(
            *degree,
            poles,
            weights.clone(),
            knots,
            multiplicities,
            *periodic,
        ),
        None => BSplineCurve2D::new(*degree, poles, knots, multiplicities, *periodic),
    };
    result.unwrap_or_else(|e| panic!("{}: failed to build BSplineCurve2D: {e}", case.name))
}

#[test]
fn bspline_shadows_match_golden_fixture() {
    for case in &bspline_cases() {
        let curve = shadow(case);
        let samples: Vec<CurveSample<[f64; 2]>> = case
            .samples
            .iter()
            .map(|sample| {
                let drop_z = |[x, y, _]: [f64; 3]| [x, y];
                CurveSample {
                    u: sample.u,
                    point: drop_z(sample.point),
                    d1: drop_z(sample.d1),
                    d2: drop_z(sample.d2),
                    d3: sample.d3.map(drop_z),
                }
            })
            .collect();
        assert_curve_2d_samples(&curve, &samples);
    }
}

// ---- surfaces_analytic.json: the five analytic surface types ----

/// Replays the `kind` surfaces listed under `key` in
/// `surfaces_analytic.json`: their samples and inversions.
fn replay_surfaces(key: &str, kind: &str) {
    let cases: Vec<SurfaceCase> = tagged_cases(&load("surfaces_analytic.json"), key, kind);
    assert!(!cases.is_empty(), "no {key}");
    for case in &cases {
        assert!(!case.samples.is_empty() && !case.parameters_of.is_empty());
        let surface = case.surface.decode();
        assert_surface_samples(&surface, &case.samples);
        assert_surface_parameters(&surface, &case.parameters_of);
    }
}

#[test]
fn planes_match_golden_fixture() {
    replay_surfaces("planes", "plane");
}

#[test]
fn cylinders_match_golden_fixture() {
    replay_surfaces("cylinders", "cylinder");
}

#[test]
fn cones_match_golden_fixture() {
    replay_surfaces("cones", "cone");
}

#[test]
fn spheres_match_golden_fixture() {
    replay_surfaces("spheres", "sphere");
}

#[test]
fn tori_match_golden_fixture() {
    replay_surfaces("tori", "torus");
}

// ---- construction.json: plane/cone/cylinder derived-construction goldens ----

#[derive(Deserialize)]
struct PlaneFromThreePointsCase {
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    frame: FrameJson,
}

#[derive(Deserialize)]
struct PlaneFromCoefficientsCase {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    frame: FrameJson,
}

#[derive(Deserialize)]
struct ConeTwoPointsRadiiCase {
    p1: [f64; 3],
    p2: [f64; 3],
    r1: f64,
    r2: f64,
    frame: FrameJson,
    ref_radius: f64,
    semi_angle: f64,
}

#[derive(Deserialize)]
struct CircleFrameJson {
    frame: FrameJson,
    radius: f64,
}

#[derive(Deserialize)]
struct CylinderFromCircleCase {
    circle: CircleFrameJson,
    frame: FrameJson,
    radius: f64,
}

#[derive(Deserialize)]
struct SurfaceConstructionFixture {
    planes_from_three_points: Vec<PlaneFromThreePointsCase>,
    planes_from_coefficients: Vec<PlaneFromCoefficientsCase>,
    cones_two_points_radii: Vec<ConeTwoPointsRadiiCase>,
    cylinders_from_circle: Vec<CylinderFromCircleCase>,
}

fn assert_frame_matches(frame: geomcore::Frame, expected: &FrameJson) {
    assert_point3(frame.origin(), expected.origin);
    assert_vector3(frame.x_direction(), expected.x_dir);
    assert_vector3(frame.y_direction(), expected.y_dir);
    assert_vector3(frame.z_direction(), expected.z_dir);
}

#[test]
fn planes_from_three_points_match_golden_fixture() {
    let fixture: SurfaceConstructionFixture =
        serde_json::from_value(load("construction.json")).unwrap();

    for (i, case) in fixture.planes_from_three_points.iter().enumerate() {
        let plane = Plane::from_three_points(point3(case.p1), point3(case.p2), point3(case.p3))
            .unwrap_or_else(|e| panic!("case {i}: failed to build Plane: {e}"));
        assert_frame_matches(plane.frame(), &case.frame);
    }
}

#[test]
fn planes_from_coefficients_match_golden_fixture() {
    let fixture: SurfaceConstructionFixture =
        serde_json::from_value(load("construction.json")).unwrap();

    for (i, case) in fixture.planes_from_coefficients.iter().enumerate() {
        let plane = Plane::from_coefficients(case.a, case.b, case.c, case.d)
            .unwrap_or_else(|e| panic!("case {i}: failed to build Plane: {e}"));
        assert_frame_matches(plane.frame(), &case.frame);
    }
}

#[test]
fn cones_from_two_points_and_radii_match_golden_fixture() {
    let fixture: SurfaceConstructionFixture =
        serde_json::from_value(load("construction.json")).unwrap();

    for (i, case) in fixture.cones_two_points_radii.iter().enumerate() {
        let cone =
            Cone::from_two_points_and_radii(point3(case.p1), point3(case.p2), case.r1, case.r2)
                .unwrap_or_else(|e| panic!("case {i}: failed to build Cone: {e}"));
        assert_close(cone.ref_radius(), case.ref_radius);
        assert_close(cone.semi_angle(), case.semi_angle);
        assert_frame_matches(cone.frame(), &case.frame);
    }
}

#[test]
fn cylinders_from_circle_match_golden_fixture() {
    let fixture: SurfaceConstructionFixture =
        serde_json::from_value(load("construction.json")).unwrap();

    for case in &fixture.cylinders_from_circle {
        let circle = Circle::from_frame(frame3(&case.circle.frame), case.circle.radius).unwrap();
        let cylinder = Cylinder::from_circle(&circle);
        assert_close(cylinder.radius(), case.radius);
        assert_frame_matches(cylinder.frame(), &case.frame);
    }
}

// ---- surfaces_bspline.json: tensor-product B-spline surfaces ----

#[test]
fn bspline_surfaces_match_golden_fixture() {
    let cases: Vec<SurfaceCase> = tagged_cases(&load("surfaces_bspline.json"), "cases", "bspline");
    assert!(!cases.is_empty());
    for case in &cases {
        assert!(!case.samples.is_empty(), "{}", case.name);
        assert_surface_samples(&case.surface.decode(), &case.samples);
    }
}

// ---- parametrize.json: analytic curve-on-surface parametrization ----

#[derive(Deserialize)]
struct ParametrizeExpect {
    kind: String,
    origin: Option<[f64; 2]>,
    direction: Option<[f64; 2]>,
    center: Option<[f64; 2]>,
    x_dir: Option<[f64; 2]>,
    y_dir: Option<[f64; 2]>,
    radius: Option<f64>,
    major_radius: Option<f64>,
    minor_radius: Option<f64>,
    focal: Option<f64>,
    degree: Option<usize>,
    periodic: Option<bool>,
    poles: Option<Vec<[f64; 2]>>,
    weights: Option<Vec<f64>>,
    knots: Option<Vec<f64>>,
    mults: Option<Vec<u32>>,
}

#[derive(Deserialize)]
struct ParametrizeNormalized {
    origin: [f64; 2],
    direction: [f64; 2],
}

#[derive(Deserialize)]
struct ParametrizeCase {
    name: String,
    /// A `SurfaceJson` tagged under `kind`.
    surface: serde_json::Value,
    /// A `CurveJson` tagged under `kind`.
    curve: serde_json::Value,
    expect: ParametrizeExpect,
    #[serde(default)]
    in_window: bool,
    #[serde(default)]
    consistency_params: Vec<f64>,
    #[serde(default)]
    consistency_points: Vec<[f64; 3]>,
    normalized: Option<ParametrizeNormalized>,
}

#[derive(Deserialize)]
struct ParametrizeFixture {
    cases: Vec<ParametrizeCase>,
}

fn parametrize_curve(curve: &Curve, surface: &Surface) -> Result<Curve2D, ParametrizeError> {
    let (surface, tol) = (surface.clone(), Tolerance::DEFAULT);
    match curve {
        Curve::Line(c) => c.parametrize_on(surface, tol),
        Curve::Circle(c) => c.parametrize_on(surface, tol),
        Curve::Ellipse(c) => c.parametrize_on(surface, tol),
        Curve::Parabola(c) => c.parametrize_on(surface, tol),
        Curve::Hyperbola(c) => c.parametrize_on(surface, tol),
        Curve::BSpline(c) => c.parametrize_on(surface, tol),
        other => panic!("no parametrization of {other:?}"),
    }
}

/// Asserts a 2D frame against the golden `center`, `x_dir` and `y_dir`.
#[track_caller]
fn assert_expected_frame2(frame: geomcore::Frame2D, expect: &ParametrizeExpect) {
    let center = expect.center.unwrap();
    let x_dir = expect.x_dir.unwrap();
    let y_dir = expect.y_dir.unwrap();
    assert_close(frame.origin().x, center[0]);
    assert_close(frame.origin().y, center[1]);
    assert_close(frame.x_direction().x, x_dir[0]);
    assert_close(frame.x_direction().y, x_dir[1]);
    assert_close(frame.y_direction().x, y_dir[0]);
    assert_close(frame.y_direction().y, y_dir[1]);
}

#[test]
fn parametrize_matches_golden_fixture() {
    let fixture: ParametrizeFixture = serde_json::from_value(load("parametrize.json")).unwrap();

    for case in &fixture.cases {
        let surface = kind_tagged::<SurfaceJson>(&case.surface).decode();
        let curve = kind_tagged::<CurveJson>(&case.curve).decode();
        let result = parametrize_curve(&curve, &surface);

        // not_analytic cases: must report NotAnalytic and stop.
        if case.expect.kind == "not_analytic" {
            assert_eq!(
                result,
                Err(ParametrizeError::NotAnalytic),
                "case {}: expected NotAnalytic",
                case.name
            );
            continue;
        }

        // curve_not_on_surface cases: the curve leaves the surface.
        if case.expect.kind == "curve_not_on_surface" {
            assert_eq!(
                result,
                Err(ParametrizeError::CurveNotOnSurface),
                "case {}: expected CurveNotOnSurface",
                case.name
            );
            continue;
        }

        let pcurve =
            result.unwrap_or_else(|e| panic!("case {}: unexpected error {e:?}", case.name));

        // (a) result kind matches.
        match (case.expect.kind.as_str(), &pcurve) {
            ("line2d", Curve2D::Line(_)) => {}
            ("circle2d", Curve2D::Circle(_)) => {}
            ("ellipse2d", Curve2D::Ellipse(_)) => {}
            ("parabola2d", Curve2D::Parabola(_)) => {}
            ("hyperbola2d", Curve2D::Hyperbola(_)) => {}
            ("bspline2d", Curve2D::BSpline(_)) => {}
            (k, got) => panic!("case {}: expected kind {k}, got {got:?}", case.name),
        }

        // (b) consistency invariant: surface(q(t)) ~= consistency_points[i].
        for (t, expected) in case.consistency_params.iter().zip(&case.consistency_points) {
            let uv = pcurve.eval_point(*t);
            let on_surface = surface.eval_point(uv.x, uv.y);
            assert_point3(on_surface, *expected);
        }

        // (c) q(0) lies inside the canonical parameter window.
        let q0 = pcurve.eval_point(0.0);
        match case.surface["kind"].as_str().unwrap() {
            "cylinder" | "cone" => {
                assert!(
                    (-1e-9..=std::f64::consts::TAU + 1e-9).contains(&q0.x),
                    "case {}: u(0) = {} out of [0, 2pi)",
                    case.name,
                    q0.x
                );
            }
            "sphere" => {
                assert!(
                    (-1e-9..=std::f64::consts::TAU + 1e-9).contains(&q0.x),
                    "case {}: u(0) = {} out of [0, 2pi)",
                    case.name,
                    q0.x
                );
                assert!(
                    (-std::f64::consts::FRAC_PI_2 - 1e-9..=std::f64::consts::FRAC_PI_2 + 1e-9)
                        .contains(&q0.y),
                    "case {}: v(0) = {} out of [-pi/2, pi/2]",
                    case.name,
                    q0.y
                );
            }
            "torus" => {
                assert!(
                    (-1e-9..=std::f64::consts::TAU + 1e-9).contains(&q0.x),
                    "case {}: u(0) = {} out of [0, 2pi)",
                    case.name,
                    q0.x
                );
                assert!(
                    (-1e-9..=std::f64::consts::TAU + 1e-9).contains(&q0.y),
                    "case {}: v(0) = {} out of [0, 2pi)",
                    case.name,
                    q0.y
                );
            }
            _ => {}
        }

        // (d) literal comparison where the raw golden is already canonical.
        // Plane cases have no window, so compare literally always.
        if case.surface["kind"] == "plane" {
            match &pcurve {
                Curve2D::Line(l) => {
                    let o = case.expect.origin.unwrap();
                    let d = case.expect.direction.unwrap();
                    assert_close(l.origin().x, o[0]);
                    assert_close(l.origin().y, o[1]);
                    assert_close(l.direction().x, d[0]);
                    assert_close(l.direction().y, d[1]);
                }
                Curve2D::Circle(c) => {
                    assert_expected_frame2(c.frame(), &case.expect);
                    assert_close(c.radius(), case.expect.radius.unwrap());
                }
                Curve2D::Ellipse(e) => {
                    assert_expected_frame2(e.frame(), &case.expect);
                    assert_close(e.major_radius(), case.expect.major_radius.unwrap());
                    assert_close(e.minor_radius(), case.expect.minor_radius.unwrap());
                }
                Curve2D::Parabola(p) => {
                    assert_expected_frame2(p.frame(), &case.expect);
                    assert_close(p.focal(), case.expect.focal.unwrap());
                }
                Curve2D::Hyperbola(h) => {
                    assert_expected_frame2(h.frame(), &case.expect);
                    assert_close(h.major_radius(), case.expect.major_radius.unwrap());
                    assert_close(h.minor_radius(), case.expect.minor_radius.unwrap());
                }
                Curve2D::BSpline(b) => {
                    let expect = &case.expect;
                    assert_eq!(b.degree(), expect.degree.unwrap(), "case {}", case.name);
                    assert_eq!(b.is_periodic(), expect.periodic.unwrap_or(false));
                    assert_eq!(b.knots(), expect.knots.as_ref().unwrap().as_slice());
                    assert_eq!(
                        b.multiplicities(),
                        expect.mults.as_ref().unwrap().as_slice()
                    );
                    assert_eq!(b.weights(), expect.weights.as_deref());
                    let poles = expect.poles.as_ref().unwrap();
                    assert_eq!(b.poles().len(), poles.len(), "case {}", case.name);
                    for (pole, expected) in b.poles().iter().zip(poles) {
                        assert_point2(*pole, *expected);
                    }
                }
                other => panic!("case {}: unexpected plane image {other:?}", case.name),
            }
        } else if case.expect.kind == "line2d" && case.in_window && case.normalized.is_none() {
            // Periodic-surface line whose raw projection was already in-window:
            // raw == normalized, so compare origin/direction literally.
            if let Curve2D::Line(l) = &pcurve {
                let o = case.expect.origin.unwrap();
                let d = case.expect.direction.unwrap();
                assert_close(l.origin().x, o[0]);
                assert_close(l.origin().y, o[1]);
                assert_close(l.direction().x, d[0]);
                assert_close(l.direction().y, d[1]);
            }
        }

        // (e) sphere meridian cases: compare against the normalized golden.
        if let (Some(norm), Curve2D::Line(l)) = (&case.normalized, &pcurve) {
            assert_close(l.origin().x, norm.origin[0]);
            assert_close(l.origin().y, norm.origin[1]);
            assert_close(l.direction().x, norm.direction[0]);
            assert_close(l.direction().y, norm.direction[1]);
        }
    }
}
