//! Python bindings for the `geomcore` geometric kernel.
//!
//! The extension module mirrors the Rust crate's namespaces: value types live
//! at the package root (`geomcore`), curves under `geomcore.curves`, and
//! surfaces under `geomcore.surfaces`. Fallible operations raise
//! the `GeomcoreError` subclass matching the Rust error type (see `errors`);
//! `GeomcoreError` is a `ValueError`.

mod errors;

use errors::invalid_argument;

use pyo3::prelude::*;
use pyo3::types::PyModule;

use geomcore::{Continuity, PointParametrization, Tolerance};

/// The tolerance given, or the default.
pub(crate) fn tolerance(tolerance: Option<PyTolerance>) -> Tolerance {
    tolerance.map_or(Tolerance::DEFAULT, |given| given.0)
}

/// A [`Continuity`] as Python sees it: the `int` order of `Order(k)`, or
/// `math.inf` for `Smooth`, so `c >= 2` and `min(...)` work directly.
pub(crate) struct ContinuityNumber(pub(crate) Continuity);

impl<'py> IntoPyObject<'py> for ContinuityNumber {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = std::convert::Infallible;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        Ok(match self.0 {
            Continuity::Order(order) => order.into_pyobject(py)?.into_any(),
            Continuity::Smooth => f64::INFINITY.into_pyobject(py)?.into_any(),
        })
    }
}

/// Checks a curve derivative order is at least 1.
pub(crate) fn check_curve_order(order: u32) -> PyResult<()> {
    if order == 0 {
        Err(invalid_argument(
            "order must be >= 1 (use eval_point for order 0)",
        ))
    } else {
        Ok(())
    }
}

/// Checks a surface derivative order satisfies `1 <= du + dv <= 2`.
pub(crate) fn check_surface_order(du: u32, dv: u32) -> PyResult<()> {
    match du.checked_add(dv) {
        Some(1 | 2) => Ok(()),
        Some(0) => Err(invalid_argument(
            "du + dv must be >= 1 (use eval_point for the (0, 0) order)",
        )),
        _ => Err(invalid_argument(format!(
            "order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
        ))),
    }
}

pub(crate) fn parse_parametrization(s: &str) -> PyResult<PointParametrization> {
    match s {
        "chordal" => Ok(PointParametrization::Chordal),
        "centripetal" => Ok(PointParametrization::Centripetal),
        _ => Err(invalid_argument(
            "parametrization must be 'chordal' or 'centripetal'",
        )),
    }
}

mod arguments;
mod bounding_box;
mod curves;
mod discretization;
mod intersect;
mod method_sets;
mod results;
mod surfaces;
mod values;

use arguments::{curve_2d_classes, curve_classes, surface_classes};
use bounding_box::{PyBoundingBox, PyBoundingBox2D};
use curves::py_parametrize_numeric;
use curves::{PyTrimmedCurve, PyTrimmedCurve2D};
use discretization::{PyDiscretization, PyDiscretization2D};
use intersect::py_curve_curve_extrema;
use results::{
    PyCompleteness, PyCompletenessStatus, PyCurveCurveExtrema, PyCurveCurveExtremum,
    PyCurveCurveHit, PyCurveCurveHit2D, PyCurveCurveIntersection, PyCurveCurveIntersection2D,
    PyCurveCurveOverlap, PyCurveCurveRegion, PyCurveProjection, PyCurveProjection2D,
    PyCurveSurfaceHit, PyCurveSurfaceIntersection, PyCurveSurfaceRegion, PyIncompleteReason,
    PyIntersectionCurve, PyIntersectionMethod, PySurfaceCurvatures, PySurfaceProjection,
    PySurfaceSurfaceIntersection, PySurfaceSurfaceRegion,
};
use values::{
    PyAxis, PyAxis2D, PyFrame, PyFrame2D, PyPoint, PyPoint2D, PySide, PyTolerance, PyTransform,
    PyVector, PyVector2D,
};

// ---------------------------------------------------------------------------
// Module wiring
// ---------------------------------------------------------------------------

#[pymodule(name = "geomcore")]
fn geomcore_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();

    errors::register(m)?;
    m.add_class::<PyPoint>()?;
    m.add_class::<PyPoint2D>()?;
    m.add_class::<PyVector>()?;
    m.add_class::<PyVector2D>()?;
    m.add_class::<PyAxis2D>()?;
    m.add_class::<PyAxis>()?;
    m.add_class::<PyFrame2D>()?;
    m.add_class::<PyFrame>()?;
    m.add_class::<PyTolerance>()?;
    m.add_class::<PyTransform>()?;
    m.add_class::<PyBoundingBox>()?;
    m.add_class::<PyBoundingBox2D>()?;
    m.add_class::<PyDiscretization2D>()?;
    m.add_class::<PyDiscretization>()?;
    m.add_class::<PySide>()?;

    m.add_class::<PyIntersectionMethod>()?;
    m.add_class::<PyCompletenessStatus>()?;
    m.add_class::<PyIncompleteReason>()?;
    m.add_class::<PyCompleteness>()?;
    m.add_class::<PyCurveCurveRegion>()?;
    m.add_class::<PyCurveSurfaceRegion>()?;
    m.add_class::<PySurfaceSurfaceRegion>()?;
    m.add_class::<PyCurveCurveHit>()?;
    m.add_class::<PyCurveCurveHit2D>()?;
    m.add_class::<PyCurveSurfaceHit>()?;
    m.add_class::<PyCurveCurveOverlap>()?;
    m.add_class::<PyIntersectionCurve>()?;
    m.add_class::<PyCurveCurveIntersection>()?;
    m.add_class::<PyCurveCurveIntersection2D>()?;
    m.add_class::<PyCurveSurfaceIntersection>()?;
    m.add_class::<PySurfaceSurfaceIntersection>()?;
    m.add_class::<PyCurveProjection>()?;
    m.add_class::<PyCurveProjection2D>()?;
    m.add_class::<PySurfaceProjection>()?;
    m.add_class::<PyCurveCurveExtremum>()?;
    m.add_class::<PyCurveCurveExtrema>()?;
    m.add_class::<PySurfaceCurvatures>()?;

    let curves = PyModule::new(py, "curves")?;
    curve_classes::add_to(&curves)?;
    curve_2d_classes::add_to(&curves)?;
    curves.add_class::<PyTrimmedCurve>()?;
    curves.add_class::<PyTrimmedCurve2D>()?;
    m.add_submodule(&curves)?;

    let surfaces = PyModule::new(py, "surfaces")?;
    surface_classes::add_to(&surfaces)?;
    m.add_submodule(&surfaces)?;

    // Register the submodules in sys.modules so that
    // `from geomcore.curves import Circle` works.
    let sys_modules = py.import("sys")?.getattr("modules")?;
    sys_modules.set_item("geomcore.curves", &curves)?;
    sys_modules.set_item("geomcore.surfaces", &surfaces)?;

    m.add_function(wrap_pyfunction!(py_parametrize_numeric, m)?)?;
    m.add_function(wrap_pyfunction!(py_curve_curve_extrema, m)?)?;

    Ok(())
}
