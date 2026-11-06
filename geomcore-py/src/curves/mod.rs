//! The curve classes, one file per class as in the core crate's
//! `curves/`.

mod bezier;
mod bezier_2d;
mod bspline;
mod bspline_2d;
mod circle;
mod circle_2d;
mod ellipse;
mod ellipse_2d;
mod hyperbola;
mod hyperbola_2d;
mod line;
mod line_2d;
mod parabola;
mod parabola_2d;
mod trimmed;

pub(crate) use bezier::PyBezierCurve;
pub(crate) use bezier_2d::PyBezierCurve2D;
pub(crate) use bspline::PyBSplineCurve;
pub(crate) use bspline_2d::PyBSplineCurve2D;
pub(crate) use circle::PyCircle;
pub(crate) use circle_2d::PyCircle2D;
pub(crate) use ellipse::PyEllipse;
pub(crate) use ellipse_2d::PyEllipse2D;
pub(crate) use hyperbola::PyHyperbola;
pub(crate) use hyperbola_2d::PyHyperbola2D;
pub(crate) use line::PyLine;
pub(crate) use line_2d::PyLine2D;
pub(crate) use parabola::PyParabola;
pub(crate) use parabola_2d::PyParabola2D;
pub(crate) use trimmed::{PyTrimmedCurve, PyTrimmedCurve2D};

/// The repr of a B-spline curve: its degree, pole count and flags, and
/// its knots and multiplicities while there are at most
/// `REPR_KNOT_LIMIT` knots.
fn bspline_repr(
    class: &str,
    (degree, pole_count): (usize, usize),
    (rational, periodic): (bool, bool),
    knots: &[f64],
    multiplicities: &[u32],
) -> String {
    const REPR_KNOT_LIMIT: usize = 8;
    let knots = if knots.len() <= REPR_KNOT_LIMIT {
        let values: Vec<String> = knots
            .iter()
            .map(|&k| crate::values::float_repr(k))
            .collect();
        format!(
            ", knots=[{}], multiplicities={multiplicities:?}",
            values.join(", ")
        )
    } else {
        String::new()
    };
    format!(
        "{class}(degree={degree}, pole_count={pole_count}, rational={}, periodic={}{knots})",
        crate::values::py_bool(rational),
        crate::values::py_bool(periodic),
    )
}

/// A Python pole index, negative counting from the end, as an index into
/// `count` poles, or an `IndexError`.
fn pole_index(index: isize, count: usize) -> pyo3::PyResult<usize> {
    let resolved = if index < 0 {
        index + count as isize
    } else {
        index
    };
    if (0..count as isize).contains(&resolved) {
        Ok(resolved as usize)
    } else {
        Err(pyo3::exceptions::PyIndexError::new_err(format!(
            "pole index {index} is out of range for {count} poles"
        )))
    }
}
