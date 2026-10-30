//! The Python classes `geomcore.Discretization` and `Discretization2D`, and
//! the argument handling the curve classes share for their `discretize` and
//! `point_count` methods.

use crate::errors::python_error;
use crate::values::*;
use geomcore::prelude::*;
use geomcore::{Discretization, Discretization2D, DiscretizationCriteria};
use pyo3::prelude::*;

/// A discretized 3D curve: points on the curve in order of parameter,
/// joined by straight segments. `parameters[i]` is the parameter of
/// `points[i]`; segment `i` joins `points[i]` and `points[i + 1]` and
/// carries a proven upper bound on its chord deviation (the largest
/// distance from the curve to the segment) and angle deviation (the
/// largest angle between the curve's tangent and the segment, π where no
/// bound exists).
#[pyclass(name = "Discretization", module = "geomcore", frozen, eq)]
#[derive(PartialEq)]
pub(crate) struct PyDiscretization(pub(crate) Discretization);

#[pymethods]
impl PyDiscretization {
    /// The parameter of each point, increasing.
    #[getter]
    fn parameters(&self) -> Vec<f64> {
        self.0.parameters.clone()
    }

    /// The points on the curve.
    #[getter]
    fn points(&self) -> Vec<PyPoint> {
        self.0.points.iter().copied().map(PyPoint).collect()
    }

    /// The chord deviation of each segment.
    #[getter]
    fn chord_deviations(&self) -> Vec<f64> {
        self.0.chord_deviations.clone()
    }

    /// The angle deviation of each segment, in radians.
    #[getter]
    fn angle_deviations(&self) -> Vec<f64> {
        self.0.angle_deviations.clone()
    }

    fn __repr__(&self) -> String {
        format!("Discretization({} points)", self.0.points.len())
    }
}

/// A discretized 2D curve, such as a curve in a surface's parameter
/// plane: the 2D counterpart of `Discretization`, its chord deviations
/// being distances in the plane.
#[pyclass(name = "Discretization2D", module = "geomcore", frozen, eq)]
#[derive(PartialEq)]
pub(crate) struct PyDiscretization2D(pub(crate) Discretization2D);

#[pymethods]
impl PyDiscretization2D {
    /// The parameter of each point, increasing.
    #[getter]
    fn parameters(&self) -> Vec<f64> {
        self.0.parameters.clone()
    }

    /// The points on the curve.
    #[getter]
    fn points(&self) -> Vec<PyPoint2D> {
        self.0.points.iter().copied().map(PyPoint2D).collect()
    }

    /// The chord deviation of each segment.
    #[getter]
    fn chord_deviations(&self) -> Vec<f64> {
        self.0.chord_deviations.clone()
    }

    /// The angle deviation of each segment, in radians.
    #[getter]
    fn angle_deviations(&self) -> Vec<f64> {
        self.0.angle_deviations.clone()
    }

    fn __repr__(&self) -> String {
        format!("Discretization2D({} points)", self.0.points.len())
    }
}

/// The criteria from the keyword arguments of `discretize`.
pub(crate) fn discretization_criteria(
    max_chord_deviation: f64,
    max_angle_deviation: Option<f64>,
) -> PyResult<DiscretizationCriteria> {
    let criteria = DiscretizationCriteria::new(max_chord_deviation).map_err(python_error)?;
    match max_angle_deviation {
        Some(angle) => criteria
            .with_max_angle_deviation(angle)
            .map_err(python_error),
        None => Ok(criteria),
    }
}

pub(crate) fn discretize_3d(
    curve: &impl Discretize,
    max_chord_deviation: f64,
    max_angle_deviation: Option<f64>,
) -> PyResult<PyDiscretization> {
    let criteria = discretization_criteria(max_chord_deviation, max_angle_deviation)?;
    Ok(PyDiscretization(curve.discretize(criteria)))
}

pub(crate) fn discretize_2d(
    curve: &impl Discretize2D,
    max_chord_deviation: f64,
    max_angle_deviation: Option<f64>,
) -> PyResult<PyDiscretization2D> {
    let criteria = discretization_criteria(max_chord_deviation, max_angle_deviation)?;
    Ok(PyDiscretization2D(curve.discretize(criteria)))
}

pub(crate) fn point_count_3d(
    curve: &impl Discretize,
    max_chord_deviation: f64,
    max_angle_deviation: Option<f64>,
) -> PyResult<usize> {
    let criteria = discretization_criteria(max_chord_deviation, max_angle_deviation)?;
    Ok(curve.point_count(criteria))
}

pub(crate) fn point_count_2d(
    curve: &impl Discretize2D,
    max_chord_deviation: f64,
    max_angle_deviation: Option<f64>,
) -> PyResult<usize> {
    let criteria = discretization_criteria(max_chord_deviation, max_angle_deviation)?;
    Ok(curve.point_count(criteria))
}
