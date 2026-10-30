//! Shapes recorded in JSON fixtures, decoded into a `Curve`, `Curve2D` or
//! `Surface` from a `type` tag, and the checks that replay a fixture's
//! samples, inversions and projections on any shape through the shape
//! traits.

use super::{Frame2DJson, FrameJson, assert_close, assert_point2, assert_point3, assert_vector2};
use super::{assert_vector3, frame2, frame3, point2, point3, vector3};
use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve, BSplineSurface, BezierCurve, Circle, Cone, Curve, Curve2D, Cylinder,
    Ellipse, Ellipse2D, ExtrusionSurface, Hyperbola, Hyperbola2D, Line, Parabola, Parabola2D,
    Plane, Point, RevolutionSurface, Sphere, Surface, Tolerance, Torus,
};
use serde::Deserialize;
use serde::de::DeserializeOwned;

/// The cases listed under `key` in a fixture, each tagged `"type": kind`
/// first, for a fixture that groups its cases by kind instead of tagging
/// them.
pub fn tagged_cases<T: DeserializeOwned>(
    fixture: &serde_json::Value,
    key: &str,
    kind: &str,
) -> Vec<T> {
    let cases = fixture[key]
        .as_array()
        .unwrap_or_else(|| panic!("the fixture has no case list `{key}`"));
    cases
        .iter()
        .map(|case| {
            let mut case = case.clone();
            case["type"] = kind.into();
            serde_json::from_value(case).unwrap_or_else(|e| panic!("a `{key}` case: {e}"))
        })
        .collect()
}

/// A shape recorded with its tag under `kind` rather than `type`.
pub fn kind_tagged<T: DeserializeOwned>(shape: &serde_json::Value) -> T {
    let mut shape = shape.clone();
    shape["type"] = shape["kind"].take();
    serde_json::from_value(shape).unwrap_or_else(|e| panic!("a shape tagged by kind: {e}"))
}

fn points(poles: &[[f64; 3]]) -> Vec<Point> {
    poles.iter().map(|&pole| point3(pole)).collect()
}

/// A curve in space, tagged by its `type`.
#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub enum CurveJson {
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
        #[serde(alias = "mults")]
        multiplicities: Vec<u32>,
        #[serde(default)]
        periodic: bool,
    },
    Bezier {
        degree: usize,
        poles: Vec<[f64; 3]>,
        weights: Option<Vec<f64>>,
    },
}

impl CurveJson {
    /// The recorded curve; panics if it is not a valid one.
    pub fn decode(&self) -> Curve {
        match self {
            CurveJson::Line { origin, direction } => {
                Line::new(point3(*origin), vector3(*direction))
                    .unwrap()
                    .into()
            }
            CurveJson::Circle { frame, radius } => {
                Circle::from_frame(frame3(frame), *radius).unwrap().into()
            }
            CurveJson::Ellipse {
                frame,
                major_radius,
                minor_radius,
            } => Ellipse::from_frame(frame3(frame), *major_radius, *minor_radius)
                .unwrap()
                .into(),
            CurveJson::Parabola { frame, focal } => {
                Parabola::from_frame(frame3(frame), *focal).unwrap().into()
            }
            CurveJson::Hyperbola {
                frame,
                major_radius,
                minor_radius,
            } => Hyperbola::from_frame(frame3(frame), *major_radius, *minor_radius)
                .unwrap()
                .into(),
            CurveJson::Bspline {
                degree,
                poles,
                weights,
                knots,
                multiplicities,
                periodic,
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
                        *periodic,
                    ),
                    None => BSplineCurve::new(*degree, poles, knots, multiplicities, *periodic),
                }
                .unwrap()
                .into()
            }
            CurveJson::Bezier {
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
}

/// A curve in the plane, tagged by its `type`: the kinds the fixtures
/// record.
#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub enum Curve2DJson {
    Ellipse {
        frame: Frame2DJson,
        major_radius: f64,
        minor_radius: f64,
    },
    Parabola {
        frame: Frame2DJson,
        focal: f64,
    },
    Hyperbola {
        frame: Frame2DJson,
        major_radius: f64,
        minor_radius: f64,
    },
}

impl Curve2DJson {
    /// The recorded curve; panics if it is not a valid one.
    pub fn decode(&self) -> Curve2D {
        match self {
            Curve2DJson::Ellipse {
                frame,
                major_radius,
                minor_radius,
            } => Ellipse2D::from_frame(frame2(frame), *major_radius, *minor_radius)
                .unwrap()
                .into(),
            Curve2DJson::Parabola { frame, focal } => Parabola2D::from_frame(frame2(frame), *focal)
                .unwrap()
                .into(),
            Curve2DJson::Hyperbola {
                frame,
                major_radius,
                minor_radius,
            } => Hyperbola2D::from_frame(frame2(frame), *major_radius, *minor_radius)
                .unwrap()
                .into(),
        }
    }
}

/// An axis: an origin and a direction.
#[derive(Deserialize)]
pub struct AxisJson {
    pub origin: [f64; 3],
    pub direction: [f64; 3],
}

/// A surface, tagged by its `type`.
#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub enum SurfaceJson {
    Plane {
        frame: FrameJson,
    },
    Cylinder {
        frame: FrameJson,
        radius: f64,
    },
    Cone {
        frame: FrameJson,
        ref_radius: f64,
        semi_angle: f64,
    },
    Sphere {
        frame: FrameJson,
        radius: f64,
    },
    Torus {
        frame: FrameJson,
        major_radius: f64,
        minor_radius: f64,
    },
    Extrusion {
        profile: CurveJson,
        direction: [f64; 3],
    },
    Revolution {
        profile: CurveJson,
        axis: AxisJson,
    },
    Bspline {
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Vec<[f64; 3]>>,
        weights: Option<Vec<Vec<f64>>>,
        u_knots: Vec<f64>,
        v_knots: Vec<f64>,
        #[serde(alias = "u_mults")]
        u_multiplicities: Vec<u32>,
        #[serde(alias = "v_mults")]
        v_multiplicities: Vec<u32>,
        #[serde(default)]
        u_periodic: bool,
        #[serde(default)]
        v_periodic: bool,
    },
}

impl SurfaceJson {
    /// The recorded surface; panics if it is not a valid one.
    pub fn decode(&self) -> Surface {
        match self {
            SurfaceJson::Plane { frame } => Plane::from_frame(frame3(frame)).into(),
            SurfaceJson::Cylinder { frame, radius } => {
                Cylinder::from_frame(frame3(frame), *radius).unwrap().into()
            }
            SurfaceJson::Cone {
                frame,
                ref_radius,
                semi_angle,
            } => Cone::from_frame(frame3(frame), *semi_angle, *ref_radius)
                .unwrap()
                .into(),
            SurfaceJson::Sphere { frame, radius } => {
                Sphere::from_frame(frame3(frame), *radius).unwrap().into()
            }
            SurfaceJson::Torus {
                frame,
                major_radius,
                minor_radius,
            } => Torus::from_frame(frame3(frame), *major_radius, *minor_radius)
                .unwrap()
                .into(),
            SurfaceJson::Extrusion { profile, direction } => {
                ExtrusionSurface::new(profile.decode(), vector3(*direction))
                    .unwrap()
                    .into()
            }
            SurfaceJson::Revolution { profile, axis } => {
                let axis = Axis::new(point3(axis.origin), vector3(axis.direction)).unwrap();
                RevolutionSurface::new(profile.decode(), axis)
                    .unwrap()
                    .into()
            }
            SurfaceJson::Bspline {
                u_degree,
                v_degree,
                poles,
                weights,
                u_knots,
                v_knots,
                u_multiplicities,
                v_multiplicities,
                u_periodic,
                v_periodic,
            } => {
                let poles: Vec<Vec<Point>> = poles.iter().map(|row| points(row)).collect();
                let (u_knots, v_knots) = (u_knots.clone(), v_knots.clone());
                let (u_multiplicities, v_multiplicities) =
                    (u_multiplicities.clone(), v_multiplicities.clone());
                match weights {
                    Some(weights) => BSplineSurface::new_rational(
                        *u_degree,
                        *v_degree,
                        poles,
                        weights.clone(),
                        u_knots,
                        u_multiplicities,
                        v_knots,
                        v_multiplicities,
                        *u_periodic,
                        *v_periodic,
                    ),
                    None => BSplineSurface::new(
                        *u_degree,
                        *v_degree,
                        poles,
                        u_knots,
                        u_multiplicities,
                        v_knots,
                        v_multiplicities,
                        *u_periodic,
                        *v_periodic,
                    ),
                }
                .unwrap()
                .into()
            }
        }
    }
}

/// A recorded evaluation of a curve: the point and derivatives at `u`,
/// `P` a `[f64; 3]` in space and a `[f64; 2]` in the plane. Not every
/// fixture records the third derivative.
#[derive(Deserialize)]
pub struct CurveSample<P> {
    pub u: f64,
    pub point: P,
    pub d1: P,
    pub d2: P,
    pub d3: Option<P>,
}

/// A recorded inversion: the parameter of a point of the curve.
#[derive(Deserialize)]
pub struct ParameterOf<P> {
    pub point: P,
    pub u: f64,
}

/// A recorded projection: the query point, the parameter and closest
/// point of the curve or surface where the fixture records them, and the
/// distance.
#[derive(Deserialize)]
pub struct Projection<P> {
    pub point: P,
    pub u: Option<f64>,
    /// `None` where the closest points are not unique, as for a query on
    /// the axis of a revolution, whose closest points fill a circle.
    pub closest: Option<P>,
    pub distance: f64,
}

/// A recorded curve with what the fixture records of it.
#[derive(Deserialize)]
pub struct CurveCase<Shape, P> {
    #[serde(default)]
    pub name: String,
    #[serde(flatten)]
    pub curve: Shape,
    #[serde(default)]
    pub samples: Vec<CurveSample<P>>,
    #[serde(default)]
    pub parameter_of: Vec<ParameterOf<P>>,
    #[serde(default)]
    pub projections: Vec<Projection<P>>,
}

/// A recorded evaluation of a surface: the point and derivatives at
/// `(u, v)`.
#[derive(Deserialize)]
pub struct SurfaceSample {
    pub u: f64,
    pub v: f64,
    pub point: [f64; 3],
    pub d1u: [f64; 3],
    pub d1v: [f64; 3],
    pub d2u: [f64; 3],
    pub d2v: [f64; 3],
    pub d2uv: [f64; 3],
}

/// A recorded inversion: the parameters of a point of the surface.
#[derive(Deserialize)]
pub struct ParametersOf {
    pub point: [f64; 3],
    pub u: f64,
    pub v: f64,
}

/// A recorded surface with what the fixture records of it.
#[derive(Deserialize)]
pub struct SurfaceCase {
    #[serde(default)]
    pub name: String,
    #[serde(flatten)]
    pub surface: SurfaceJson,
    #[serde(default)]
    pub samples: Vec<SurfaceSample>,
    #[serde(default)]
    pub parameters_of: Vec<ParametersOf>,
    #[serde(default)]
    pub projections: Vec<Projection<[f64; 3]>>,
}

/// The curve evaluates to every recorded sample.
#[track_caller]
pub fn assert_curve_samples(curve: &impl ParametricCurve, samples: &[CurveSample<[f64; 3]>]) {
    for sample in samples {
        let u = sample.u;
        assert_point3(curve.eval_point(u), sample.point);
        assert_vector3(curve.eval_derivative(u, 1), sample.d1);
        assert_vector3(curve.eval_derivative(u, 2), sample.d2);
        if let Some(d3) = sample.d3 {
            assert_vector3(curve.eval_derivative(u, 3), d3);
        }
    }
}

/// The plane curve evaluates to every recorded sample.
#[track_caller]
pub fn assert_curve_2d_samples(curve: &impl ParametricCurve2D, samples: &[CurveSample<[f64; 2]>]) {
    for sample in samples {
        let u = sample.u;
        assert_point2(curve.eval_point(u), sample.point);
        assert_vector2(curve.eval_derivative(u, 1), sample.d1);
        assert_vector2(curve.eval_derivative(u, 2), sample.d2);
        if let Some(d3) = sample.d3 {
            assert_vector2(curve.eval_derivative(u, 3), d3);
        }
    }
}

/// The curve finds the recorded parameter of every recorded point.
#[track_caller]
pub fn assert_curve_parameters(curve: &impl PointExtrema, inversions: &[ParameterOf<[f64; 3]>]) {
    for inversion in inversions {
        let u = curve.parameter_of(point3(inversion.point), Tolerance::DEFAULT);
        assert_close(u, inversion.u);
    }
}

/// The plane curve finds the recorded parameter of every recorded point.
#[track_caller]
pub fn assert_curve_2d_parameters(
    curve: &impl PointExtrema2D,
    inversions: &[ParameterOf<[f64; 2]>],
) {
    for inversion in inversions {
        let u = curve.parameter_of(point2(inversion.point), Tolerance::DEFAULT);
        assert_close(u, inversion.u);
    }
}

/// The plane curve projects every recorded query onto the recorded
/// parameter, closest point and distance.
#[track_caller]
pub fn assert_curve_2d_projections(
    curve: &impl PointExtrema2D,
    projections: &[Projection<[f64; 2]>],
) {
    for projection in projections {
        let found = curve.project_point(point2(projection.point), Tolerance::DEFAULT);
        if let Some(u) = projection.u {
            assert_close(found.parameter, u);
        }
        if let Some(closest) = projection.closest {
            assert_point2(found.point, closest);
        }
        assert_close(found.distance, projection.distance);
    }
}

/// The surface evaluates to every recorded sample.
#[track_caller]
pub fn assert_surface_samples(surface: &impl ParametricSurface, samples: &[SurfaceSample]) {
    for sample in samples {
        let (u, v) = (sample.u, sample.v);
        assert_point3(surface.eval_point(u, v), sample.point);
        assert_vector3(surface.eval_derivative(u, v, 1, 0), sample.d1u);
        assert_vector3(surface.eval_derivative(u, v, 0, 1), sample.d1v);
        assert_vector3(surface.eval_derivative(u, v, 2, 0), sample.d2u);
        assert_vector3(surface.eval_derivative(u, v, 0, 2), sample.d2v);
        assert_vector3(surface.eval_derivative(u, v, 1, 1), sample.d2uv);
    }
}

/// The surface finds the recorded parameters of every recorded point.
#[track_caller]
pub fn assert_surface_parameters(surface: &impl ParametricSurface, inversions: &[ParametersOf]) {
    for inversion in inversions {
        let (u, v) = surface.parameters_of(point3(inversion.point), Tolerance::DEFAULT);
        assert_close(u, inversion.u);
        assert_close(v, inversion.v);
    }
}

/// The surface projects every recorded query at the recorded distance,
/// onto a point it evaluates to at the found parameters, and onto the
/// recorded closest point where there is one.
#[track_caller]
pub fn assert_surface_projections(
    name: &str,
    surface: &impl ParametricSurface,
    projections: &[Projection<[f64; 3]>],
) {
    for projection in projections {
        let found = surface.project_point(point3(projection.point), Tolerance::DEFAULT);
        assert!(
            (found.distance - projection.distance).abs() <= 1e-7 * projection.distance.max(1.0),
            "{name}: projection of {:?} at distance {} != {}",
            projection.point,
            found.distance,
            projection.distance
        );
        assert_point3(
            surface.eval_point(found.u, found.v),
            [found.point.x, found.point.y, found.point.z],
        );
        if let Some(closest) = projection.closest {
            assert_point3(found.point, closest);
        }
    }
}
