//! `geomcore.CurveApproximation`, `geomcore.CurveApproximation2D` and
//! `geomcore.SurfaceApproximation`: what the operations within a deviation
//! return.

use crate::curves::{PyBSplineCurve, PyBSplineCurve2D};
use crate::surfaces::PyBSplineSurface;
use crate::values::{float_repr, py_bool};
use geomcore::{CurveApproximation, CurveApproximation2D, SurfaceApproximation};
use pyo3::prelude::*;

/// A 3D B-spline curve within a deviation of what it approximates, from
/// `uv_curve_to_3d_within` and `BSplineCurve.approximate_within`. The
/// operation says what the curve is and how its deviation is measured.
#[pyclass(name = "CurveApproximation", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveApproximation(pub(crate) CurveApproximation);

#[pymethods]
impl PyCurveApproximation {
    /// The curve.
    #[getter]
    fn curve(&self) -> PyBSplineCurve {
        PyBSplineCurve(self.0.curve().clone())
    }

    /// The largest distance between the curve and what it approximates,
    /// as the operation measures it: a bound when `is_proven()`, the
    /// largest sampled distance otherwise. It is within the asked
    /// deviation, up to rounding.
    #[getter]
    fn deviation(&self) -> f64 {
        self.0.deviation()
    }

    /// Whether `deviation` is a proven bound, up to rounding, rather than
    /// sampled.
    fn is_proven(&self) -> bool {
        self.0.is_proven()
    }

    fn __repr__(&self) -> String {
        let curve = self.0.curve();
        format!(
            "CurveApproximation(degree={}, pole_count={}, deviation={}, proven={})",
            curve.degree(),
            curve.poles().len(),
            float_repr(self.0.deviation()),
            py_bool(self.0.is_proven())
        )
    }
}

/// A 2D B-spline curve within a deviation of what it approximates, from
/// `curve_3d_to_uv_within` and `BSplineCurve2D.approximate_within`: as
/// `CurveApproximation` is for a 3D curve.
#[pyclass(name = "CurveApproximation2D", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveApproximation2D(pub(crate) CurveApproximation2D);

#[pymethods]
impl PyCurveApproximation2D {
    /// The curve.
    #[getter]
    fn curve(&self) -> PyBSplineCurve2D {
        PyBSplineCurve2D(self.0.curve().clone())
    }

    /// The largest distance between the curve and what it approximates,
    /// as the operation measures it: a bound when `is_proven()`, the
    /// largest sampled distance otherwise. It is within the asked
    /// deviation, up to rounding.
    #[getter]
    fn deviation(&self) -> f64 {
        self.0.deviation()
    }

    /// Whether `deviation` is a proven bound, up to rounding, rather than
    /// sampled.
    fn is_proven(&self) -> bool {
        self.0.is_proven()
    }

    fn __repr__(&self) -> String {
        let curve = self.0.curve();
        format!(
            "CurveApproximation2D(degree={}, pole_count={}, deviation={}, proven={})",
            curve.degree(),
            curve.poles().len(),
            float_repr(self.0.deviation()),
            py_bool(self.0.is_proven())
        )
    }
}

/// A B-spline surface within a deviation of what it approximates, from
/// `BSplineSurface.approximate_within`: as `CurveApproximation` is for a
/// curve.
#[pyclass(name = "SurfaceApproximation", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PySurfaceApproximation(pub(crate) SurfaceApproximation);

#[pymethods]
impl PySurfaceApproximation {
    /// The surface.
    #[getter]
    fn surface(&self) -> PyBSplineSurface {
        PyBSplineSurface(self.0.surface().clone())
    }

    /// The largest distance between the surface and what it approximates,
    /// as the operation measures it: a bound when `is_proven()`, the
    /// largest sampled distance otherwise. It is within the asked
    /// deviation, up to rounding.
    #[getter]
    fn deviation(&self) -> f64 {
        self.0.deviation()
    }

    /// Whether `deviation` is a proven bound, up to rounding, rather than
    /// sampled.
    fn is_proven(&self) -> bool {
        self.0.is_proven()
    }

    fn __repr__(&self) -> String {
        let surface = self.0.surface();
        format!(
            "SurfaceApproximation(u_degree={}, v_degree={}, pole_count=({}, {}), deviation={}, proven={})",
            surface.u_degree(),
            surface.v_degree(),
            surface.u_pole_count(),
            surface.v_pole_count(),
            float_repr(self.0.deviation()),
            py_bool(self.0.is_proven())
        )
    }
}
