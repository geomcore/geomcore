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

use geomcore::{Continuity, Tolerance};
use values::ToleranceArgument;

/// The tolerance given, or the default.
pub(crate) fn tolerance(tolerance: Option<ToleranceArgument>) -> Tolerance {
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
            "order must be >= 1 (use point_at for order 0)",
        ))
    } else {
        Ok(())
    }
}

/// The highest total order `du + dv` a core surface evaluates.
pub(crate) trait HighestDerivativeOrder {
    const HIGHEST_DERIVATIVE_ORDER: u32;
}

macro_rules! highest_derivative_order {
    ($order:expr => $($surface:ty),+ $(,)?) => {
        $(impl HighestDerivativeOrder for $surface {
            const HIGHEST_DERIVATIVE_ORDER: u32 = $order;
        })+
    };
}

highest_derivative_order!(u32::MAX =>
    geomcore::Plane, geomcore::Cylinder, geomcore::Cone, geomcore::Sphere, geomcore::Torus,
);
highest_derivative_order!(2 =>
    geomcore::ExtrusionSurface, geomcore::RevolutionSurface, geomcore::BSplineSurface,
    geomcore::BezierSurface,
);

/// Checks a surface derivative order satisfies
/// `1 <= du + dv <= highest`: every order on the plane, cylinder, cone,
/// sphere and torus, 2 on the other surfaces.
pub(crate) fn check_surface_order<S: HighestDerivativeOrder>(
    _surface: &S,
    du: u32,
    dv: u32,
) -> PyResult<()> {
    let highest = S::HIGHEST_DERIVATIVE_ORDER;
    match du.checked_add(dv) {
        Some(0) => Err(invalid_argument(
            "du + dv must be >= 1 (use point_at for the (0, 0) order)",
        )),
        Some(order) if order <= highest => Ok(()),
        _ => Err(invalid_argument(format!(
            "order du={du}, dv={dv} is not supported (du + dv must be <= {highest} on this surface)"
        ))),
    }
}

mod approximation;
mod arguments;
mod arrays;
mod bounding_box;
mod curves;
mod discretization;
mod enums;
mod intersect;
mod method_sets;
mod pickling;
mod results;
mod surfaces;
mod values;

use approximation::{PyCurveApproximation, PyCurveApproximation2D, PySurfaceApproximation};
use arguments::{curve_2d_classes, curve_classes, surface_classes};
use bounding_box::{PyBoundingBox, PyBoundingBox2D};
use curves::{PyTrimmedCurve, PyTrimmedCurve2D};
use discretization::{PyDiscretization, PyDiscretization2D};
use intersect::py_curve_curve_extrema;
use results::{
    PyCompleteness, PyCurveCurveExtrema, PyCurveCurveExtremum, PyCurveCurveHit, PyCurveCurveHit2D,
    PyCurveCurveIntersection, PyCurveCurveIntersection2D, PyCurveCurveOverlap, PyCurveCurveRegion,
    PyCurveProjection, PyCurveProjection2D, PyCurveSurfaceHit, PyCurveSurfaceIntersection,
    PyCurveSurfaceRegion, PyIntersectionCurve, PySurfaceCurvatures, PySurfaceProjection,
    PySurfaceSurfaceIntersection, PySurfaceSurfaceRegion,
};
use surfaces::{py_curve_3d_to_uv_within, py_uv_curve_to_3d_within};
use values::{
    PyAxis, PyAxis2D, PyFrame, PyFrame2D, PyPoint, PyPoint2D, PyTolerance, PyTransform, PyVector,
    PyVector2D,
};

// ---------------------------------------------------------------------------
// Module wiring
// ---------------------------------------------------------------------------

#[pymodule(name = "geomcore")]
fn geomcore_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();

    m.setattr("__version__", env!("CARGO_PKG_VERSION"))?;
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
    enums::register(m)?;
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
    m.add_class::<PyCurveApproximation>()?;
    m.add_class::<PyCurveApproximation2D>()?;
    m.add_class::<PySurfaceApproximation>()?;

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

    m.add_function(wrap_pyfunction!(py_curve_3d_to_uv_within, m)?)?;
    m.add_function(wrap_pyfunction!(py_uv_curve_to_3d_within, m)?)?;
    m.add_function(wrap_pyfunction!(py_curve_curve_extrema, m)?)?;

    Ok(())
}
