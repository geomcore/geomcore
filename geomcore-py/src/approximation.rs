//! `geomcore.CurveApproximation`, `geomcore.CurveApproximation2D` and
//! `geomcore.SurfaceApproximation`: what the operations within a deviation
//! return; and `geomcore.CurveOffset2D`, what a 2D curve's `offset`
//! returns, exact or within a deviation.

use crate::arguments::curve_2d_classes;
use crate::curves::{PyBSplineCurve, PyBSplineCurve2D, PyTrimmedCurve2D};
use crate::surfaces::PyBSplineSurface;
use crate::values::{float_repr, py_bool};
use geomcore::{
    BSplineCurve2D, Circle2D, Curve2D, CurveApproximation, CurveApproximation2D, CurveOffset2D,
    Line2D, SurfaceApproximation, TrimmedCurve2D,
};
use pyo3::prelude::*;

/// A 3D B-spline curve within a deviation of what it approximates, from
/// `uv_curve_to_3d_within` and `BSplineCurve.approximate_within`. The
/// operation says what the curve is and how its deviation is measured.
#[pyclass(name = "CurveApproximation", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveApproximation(pub(crate) CurveApproximation);

#[pymethods]
impl PyCurveApproximation {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

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
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

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
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

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

/// The offset curve of a `CurveOffset2D`: a whole curve, or a trimmed
/// piece for the offset of a `TrimmedCurve2D`.
#[derive(Clone, PartialEq)]
pub(crate) enum OffsetCurve2D {
    Whole(Curve2D),
    Trimmed(TrimmedCurve2D),
}

/// The offset of a 2D curve, from `offset`: the offset curve, how far it is
/// from the true offset `c(t) + d n(t)` at the same parameter, whether that
/// distance is proven, and the parameters where the true offset has a cusp.
#[pyclass(name = "CurveOffset2D", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveOffset2D(pub(crate) CurveOffset2D<OffsetCurve2D>);

/// A core offset as the Python result.
pub(crate) trait IntoPyCurveOffset2D {
    fn into_py_offset(self) -> PyCurveOffset2D;
}

/// The result of `offset`, with its curve mapped by `curve`.
fn py_offset<C>(
    offset: CurveOffset2D<C>,
    curve: impl FnOnce(C) -> OffsetCurve2D,
) -> PyCurveOffset2D {
    let (deviation, proven) = (offset.deviation(), offset.is_proven());
    let cusps = offset.cusp_parameters().to_vec();
    PyCurveOffset2D(CurveOffset2D::from_parts(
        curve(offset.into_curve()),
        deviation,
        proven,
        cusps,
    ))
}

macro_rules! whole_offsets {
    ($($curve:ty),+) => {$(
        impl IntoPyCurveOffset2D for CurveOffset2D<$curve> {
            fn into_py_offset(self) -> PyCurveOffset2D {
                py_offset(self, |curve| OffsetCurve2D::Whole(curve.into()))
            }
        }
    )+};
}

whole_offsets!(Line2D, Circle2D, BSplineCurve2D, Curve2D);

impl IntoPyCurveOffset2D for CurveOffset2D<TrimmedCurve2D> {
    fn into_py_offset(self) -> PyCurveOffset2D {
        py_offset(self, OffsetCurve2D::Trimmed)
    }
}

#[pymethods]
impl PyCurveOffset2D {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The offset curve, in the parameter of the source curve: a `Line2D`
    /// for a line, a `Circle2D` for a circle, a `BSplineCurve2D` for the
    /// other curves, and a `TrimmedCurve2D` over the same range for a
    /// trimmed curve.
    #[getter]
    fn curve(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        offset_curve_to_py(py, self.0.curve())
    }

    /// The largest distance between the offset curve's point at `t` and the
    /// true offset's: zero for an exact offset, a bound when `is_proven()`,
    /// the largest sampled distance otherwise. It is within the asked
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

    /// Whether the offset curve is the true offset itself, up to rounding:
    /// a line's or a circle's.
    fn is_exact(&self) -> bool {
        self.0.is_exact()
    }

    /// The parameters on the source curve, in increasing order, where the
    /// true offset has a cusp: where the distance crosses the radius of
    /// curvature on its side, so that the offset stops and turns back.
    /// Between two cusps the offset runs against the curve and forms a
    /// loop, which is kept, not trimmed.
    #[getter]
    fn cusp_parameters(&self) -> Vec<f64> {
        self.0.cusp_parameters().to_vec()
    }

    fn __repr__(&self) -> String {
        format!(
            "CurveOffset2D(exact={}, deviation={}, proven={}, cusp_count={})",
            py_bool(self.0.is_exact()),
            float_repr(self.0.deviation()),
            py_bool(self.0.is_proven()),
            self.0.cusp_parameters().len()
        )
    }
}

/// The Python object of an offset curve.
pub(crate) fn offset_curve_to_py(py: Python<'_>, curve: &OffsetCurve2D) -> PyResult<Py<PyAny>> {
    match curve {
        OffsetCurve2D::Whole(curve) => curve_2d_classes::to_py(py, curve.clone()),
        OffsetCurve2D::Trimmed(curve) => {
            Ok(Py::new(py, PyTrimmedCurve2D(curve.clone()))?.into_any())
        }
    }
}
