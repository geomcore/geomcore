//! Fitting B-spline surfaces to grids of points: least-squares
//! approximation, with a given number of poles or within a deviation.
//!
//! On a full grid the tensor-product least-squares fit separates: each
//! grid row is fitted along `v`, then each column of those fits' poles
//! along `u`, every row (column) on the same parameters and knots, so one
//! factored system per direction serves them all (see
//! [`crate::BSplineCurve::approximate`] for the curve fit).

use crate::curves::BSplineConstructionError;
use crate::curves::fitting::{
    BasisTable, approximation_knots, averaged_knots, check_degree, fit_within_fewest,
    interpolation, least_squares, shared_parameters,
};
use crate::{BSplineSurface, Point, PointParametrization};

/// A grid of points `points[i][j]`, `i` along `u` and `j` along `v`,
/// checked rectangular, with its parameters in both directions.
struct Grid<'a> {
    points: &'a [Vec<Point>],
    rows: usize,
    cols: usize,
    u_params: Vec<f64>,
    v_params: Vec<f64>,
}

impl<'a> Grid<'a> {
    fn new(
        points: &'a [Vec<Point>],
        parametrization: PointParametrization,
    ) -> Result<Self, BSplineConstructionError> {
        let rows = points.len();
        let cols = points.first().map_or(0, Vec::len);
        if rows < 2 || cols < 2 {
            return Err(BSplineConstructionError::TooFewPoints);
        }
        if points.iter().any(|row| row.len() != cols) {
            return Err(BSplineConstructionError::PoleCountMismatch);
        }
        // Along `u`, the columns are the lines of points; along `v`, the rows.
        let u_params = shared_parameters(cols, rows, |j, i| points[i][j], parametrization)?;
        let v_params = shared_parameters(rows, cols, |i, j| points[i][j], parametrization)?;
        Ok(Grid {
            points,
            rows,
            cols,
            u_params,
            v_params,
        })
    }

    /// The least-squares poles along `v`, `count` per grid row, on `basis`
    /// (`rows x count`, row-major).
    fn fit_rows(&self, basis: &BasisTable, count: usize) -> Option<Vec<Point>> {
        least_squares(basis, count, self.rows, |i, j| self.points[i][j])
    }

    /// The least-squares poles along `u`, `count` per column of
    /// `row_poles` (`rows x columns`, row-major), on `basis`, as the
    /// surface's grid (`count x columns`, row-major).
    fn fit_columns(
        &self,
        basis: &BasisTable,
        count: usize,
        row_poles: &[Point],
        columns: usize,
    ) -> Option<Vec<Point>> {
        let by_column = least_squares(basis, count, columns, |l, i| row_poles[i * columns + l])?;
        Some(transposed(&by_column, columns, count))
    }

    /// The deviation of each row of the grid along `v`: of grid point
    /// `(i, j)` from the point at `v` parameter `j` of the curve with
    /// `basis` and the poles of row `i` in `row_poles` (`count` per row),
    /// the most over the rows, into `deviations[j]`.
    fn row_deviations(
        &self,
        basis: &BasisTable,
        row_poles: &[Point],
        count: usize,
        deviations: &mut [f64],
    ) {
        deviations.fill(0.0);
        for (row, line) in row_poles.chunks_exact(count).zip(self.points) {
            for (j, deviation) in deviations.iter_mut().enumerate() {
                *deviation = deviation.max(basis.point(j, row).distance(line[j]));
            }
        }
    }

    /// The deviation of each grid point `(i, j)` from the point at its
    /// parameters of the surface with the `u` and `v` bases
    /// `(u_basis, v_basis)` and the poles `poles` (row-major, `columns`
    /// per row), the most over `j`, into `deviations[i]`.
    fn deviations(
        &self,
        (u_basis, v_basis): (&BasisTable, &BasisTable),
        poles: &[Point],
        columns: usize,
        deviations: &mut [f64],
    ) {
        // The poles of the `v` curve at each `u` parameter, then its points.
        let mut curve = vec![Point::ORIGIN; columns];
        for (i, deviation) in deviations.iter_mut().enumerate() {
            let (first, values) = u_basis.at(i);
            for (l, pole) in curve.iter_mut().enumerate() {
                let mut x = [0.0; 3];
                for (a, &value) in values.iter().enumerate() {
                    let q = poles[(first + a) * columns + l];
                    x[0] += value * q.x;
                    x[1] += value * q.y;
                    x[2] += value * q.z;
                }
                *pole = Point::new(x[0], x[1], x[2]);
            }
            *deviation = (0..self.cols)
                .map(|j| v_basis.point(j, &curve).distance(self.points[i][j]))
                .fold(0.0, f64::max);
        }
    }
}

/// `values` (`lines` lines of `count`) as `count` lines of `lines`.
fn transposed(values: &[Point], lines: usize, count: usize) -> Vec<Point> {
    let mut out = Vec::with_capacity(values.len());
    out.extend((0..count).flat_map(|k| (0..lines).map(move |l| values[l * count + k])));
    out
}

impl BSplineSurface {
    /// Least-squares B-spline surface fit through the grid of points
    /// `points[i][j]` (`i` along `u`, `j` along `v`), with
    /// `u_pole_count x v_pole_count` poles.
    ///
    /// Each direction's parameters average those of its lines of points
    /// (the columns for `u`, the rows for `v`), as
    /// [`crate::BSplineCurve::approximate`] gives a curve's, over the
    /// lines that have a length, and its knots follow those parameters as
    /// the curve fit's do. The fit
    /// approximates rather than hits the points.
    ///
    /// # Errors
    ///
    /// [`BSplineConstructionError::TooFewPoints`] for fewer than two rows
    /// or columns, [`BSplineConstructionError::PoleCountMismatch`] for
    /// rows of differing lengths or when the points leave a pole
    /// undetermined, [`BSplineConstructionError::InvalidDegree`] outside
    /// `1..=25`, [`BSplineConstructionError::DegreeTooHigh`] when a
    /// degree exceeds its pole count minus one,
    /// [`BSplineConstructionError::TooManyPoles`] for more poles than
    /// points in a direction, or
    /// [`BSplineConstructionError::ConfusedPoints`] when a direction's
    /// parameters coincide.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point, PointParametrization};
    /// let points: Vec<Vec<Point>> = (0..10)
    ///     .map(|i| {
    ///         (0..10)
    ///             .map(|j| {
    ///                 let (x, y) = (i as f64 / 9.0, j as f64 / 9.0);
    ///                 Point::new(x, y, x * x - y * y)
    ///             })
    ///             .collect()
    ///     })
    ///     .collect();
    /// let surface =
    ///     BSplineSurface::approximate(&points, 2, 2, 3, 3, PointParametrization::Chordal)
    ///         .unwrap();
    /// // A quadratic Bézier patch holds `x² - y²` exactly.
    /// let p = surface.eval_point(0.5, 0.5);
    /// assert!((p.z - (p.x * p.x - p.y * p.y)).abs() < 1e-12);
    /// ```
    pub fn approximate(
        points: &[Vec<Point>],
        u_degree: usize,
        v_degree: usize,
        u_pole_count: usize,
        v_pole_count: usize,
        parametrization: PointParametrization,
    ) -> Result<BSplineSurface, BSplineConstructionError> {
        check_degree(u_degree)?;
        check_degree(v_degree)?;
        let grid = Grid::new(points, parametrization)?;
        for (degree, count, available) in [
            (u_degree, u_pole_count, grid.rows),
            (v_degree, v_pole_count, grid.cols),
        ] {
            if degree > count.saturating_sub(1) {
                return Err(BSplineConstructionError::DegreeTooHigh);
            }
            if count > available {
                return Err(BSplineConstructionError::TooManyPoles);
            }
        }
        let v_flat = approximation_knots(&grid.v_params, v_pole_count, v_degree);
        let v_basis = BasisTable::new(&v_flat, v_degree, &grid.v_params);
        let u_flat = approximation_knots(&grid.u_params, u_pole_count, u_degree);
        let u_basis = BasisTable::new(&u_flat, u_degree, &grid.u_params);
        let poles = grid
            .fit_rows(&v_basis, v_pole_count)
            .and_then(|rows| grid.fit_columns(&u_basis, u_pole_count, &rows, v_pole_count))
            .ok_or(BSplineConstructionError::PoleCountMismatch)?;
        Ok(BSplineSurface::from_flat_knots(
            (u_degree, u_flat),
            (v_degree, v_flat),
            poles,
        ))
    }

    /// Least-squares B-spline surface fit through the grid of points
    /// `points[i][j]` (as [`BSplineSurface::approximate`]), with as few
    /// poles as keep every point within `max_deviation` of the surface.
    ///
    /// The deviation is measured from each point to the surface at the
    /// point's parameters (which bounds the distance to the surface).
    /// The `v` knots come first, until every row fits within 0.8 of
    /// `max_deviation`; then the `u` knots, until the whole grid fits
    /// within `max_deviation`. Each direction keeps the fewer poles of
    /// knots refined as [`crate::BSplineCurve::approximate_within`]
    /// refines a curve's and of the fewest poles that fit on knots spread
    /// along the parameters, as [`BSplineSurface::approximate`] places
    /// them. A direction no least-squares fit brings close enough
    /// interpolates its lines of points instead, one pole per point, as
    /// with a `max_deviation` of zero.
    ///
    /// # Errors
    ///
    /// As [`BSplineSurface::approximate`], a degree exceeding its
    /// direction's points minus one being
    /// [`BSplineConstructionError::DegreeTooHigh`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineSurface, Point, PointParametrization};
    /// let points: Vec<Vec<Point>> = (0..20)
    ///     .map(|i| {
    ///         (0..20)
    ///             .map(|j| {
    ///                 let (x, y) = (i as f64 / 19.0, j as f64 / 19.0);
    ///                 Point::new(x, y, (3.0 * x).sin() * (2.0 * y).cos())
    ///             })
    ///             .collect()
    ///     })
    ///     .collect();
    /// let surface = BSplineSurface::approximate_within(
    ///     &points, 3, 3, 1e-3, PointParametrization::Chordal,
    /// )
    /// .unwrap();
    /// assert!(surface.u_pole_count() < 20 && surface.v_pole_count() < 20);
    /// ```
    pub fn approximate_within(
        points: &[Vec<Point>],
        u_degree: usize,
        v_degree: usize,
        max_deviation: f64,
        parametrization: PointParametrization,
    ) -> Result<BSplineSurface, BSplineConstructionError> {
        check_degree(u_degree)?;
        check_degree(v_degree)?;
        let grid = Grid::new(points, parametrization)?;
        if u_degree > grid.rows - 1 || v_degree > grid.cols - 1 {
            return Err(BSplineConstructionError::DegreeTooHigh);
        }
        // Along `v`: every row within most of the deviation. Fewer `v`
        // poles cost few along `u`, whose fit is checked on the whole grid.
        let (q, v_params) = (v_degree, &grid.v_params);
        let fitted_rows = fit_within_fewest(
            v_params,
            q,
            0.8 * max_deviation,
            |flat, basis, deviations| {
                let count = flat.len() - q - 1;
                let poles = grid.fit_rows(basis, count)?;
                grid.row_deviations(basis, &poles, count, deviations);
                Some(poles)
            },
        );
        let (v_flat, row_poles) = fitted_rows.unwrap_or_else(|| {
            let flat = averaged_knots(v_params, q);
            let basis = BasisTable::new(&flat, q, v_params);
            (flat, interpolation(&basis, grid.rows, |i, j| points[i][j]))
        });
        let v_basis = BasisTable::new(&v_flat, q, v_params);
        let columns = row_poles.len() / grid.rows;
        // Along `u`: the whole grid within the deviation. Interpolating
        // the rows' poles leaves only the rows' deviation, within it.
        let (p, u_params) = (u_degree, &grid.u_params);
        let fitted = fit_within_fewest(u_params, p, max_deviation, |flat, basis, deviations| {
            let poles = grid.fit_columns(basis, flat.len() - p - 1, &row_poles, columns)?;
            grid.deviations((basis, &v_basis), &poles, columns, deviations);
            Some(poles)
        });
        let (u_flat, poles) = fitted.unwrap_or_else(|| {
            let flat = averaged_knots(u_params, p);
            let basis = BasisTable::new(&flat, p, u_params);
            let by_column = interpolation(&basis, columns, |l, i| row_poles[i * columns + l]);
            (flat, transposed(&by_column, columns, grid.rows))
        });
        Ok(BSplineSurface::from_flat_knots(
            (u_degree, u_flat),
            (v_degree, v_flat),
            poles,
        ))
    }
}

#[cfg(test)]
mod tests;
