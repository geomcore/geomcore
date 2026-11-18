//! The grid every tessellation is: rows of parameters in `u` and in `v`,
//! a vertex at each pair, two triangles per cell. Counting and filling it
//! is the same for every surface; only its planning differs.

use super::{TessellationBuffers, TessellationDeviation, TessellationError, TessellationSize};
use crate::geometry::surfaces::ParametricSurface;
use crate::geometry::surfaces::parametric_surface::{direct_normal, reference_widths};
use crate::{Direction, Point, Point2D, Tolerance};

/// Fails unless a grid of `columns` by `rows` vertices fits `u32`
/// indices: checked before the grid's parameters are allocated.
pub(super) fn check_count(columns: usize, rows: usize) -> Result<(), TessellationError> {
    columns
        .checked_mul(rows)
        .filter(|&count| u32::try_from(count).is_ok())
        .map(|_| ())
        .ok_or(TessellationError::TooManyVertices)
}

/// A planned tessellation: the parameters of its vertices, which of their
/// rows collapse to a point, and the deviation the plan proves.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Grid {
    /// The `u` of each column of vertices, nondecreasing: two equal
    /// neighbours split the grid there, with no cell between them.
    pub(super) us: Vec<f64>,
    /// The `v` of each row of vertices, as `us`.
    pub(super) vs: Vec<f64>,
    /// Whether each column lies within the confusion tolerance of a point.
    pub(super) collapsed_columns: Vec<bool>,
    /// Whether each row lies within the confusion tolerance of a point.
    pub(super) collapsed_rows: Vec<bool>,
    /// Whether the `u` range spans a whole period: the last column then
    /// takes the first one's points.
    pub(super) u_closed: bool,
    /// As `u_closed`, for the rows.
    pub(super) v_closed: bool,
    /// The deviation the plan proves, with every collapsed row's vertices
    /// moved onto its first one. A seam's last row takes the first row's
    /// points, which differ from its own by rounding only.
    pub(super) deviation: f64,
    /// Whether `deviation` is proven.
    pub(super) proven: bool,
}

impl Grid {
    /// A grid with no collapsed row or column and no seam.
    pub(super) fn open(us: Vec<f64>, vs: Vec<f64>, deviation: f64, proven: bool) -> Grid {
        Grid {
            collapsed_columns: vec![false; us.len()],
            collapsed_rows: vec![false; vs.len()],
            us,
            vs,
            u_closed: false,
            v_closed: false,
            deviation,
            proven,
        }
    }

    /// The vertex index of column `i`, row `j`.
    fn index(&self, i: usize, j: usize) -> usize {
        j * self.us.len() + i
    }

    /// The triangles of cell `(i, j)`, between columns `i` and `i + 1` and
    /// rows `j` and `j + 1`: none across a split, and without the one a
    /// collapsed row or column leaves without area. Its corners are
    /// `a = (i, j)`, `b = (i + 1, j)`, `c = (i + 1, j + 1)`,
    /// `d = (i, j + 1)`, so both triangles turn as `u` then `v` do.
    fn cell_triangles(&self, i: usize, j: usize) -> [Option<[usize; 3]>; 2] {
        if self.us[i] >= self.us[i + 1] || self.vs[j] >= self.vs[j + 1] {
            return [None, None];
        }
        let (a, b) = (self.index(i, j), self.index(i + 1, j));
        let (c, d) = (self.index(i + 1, j + 1), self.index(i, j + 1));
        // `abc` loses its area when `a` meets `b` or `b` meets `c`, `acd`
        // when `c` meets `d` or `d` meets `a`.
        let first = !(self.collapsed_rows[j] || self.collapsed_columns[i + 1]);
        let second = !(self.collapsed_rows[j + 1] || self.collapsed_columns[i]);
        [first.then_some([a, b, c]), second.then_some([a, c, d])]
    }

    /// The cells, column-major within each row of cells.
    fn cells(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        let columns = self.us.len() - 1;
        (0..self.vs.len() - 1).flat_map(move |j| (0..columns).map(move |i| (i, j)))
    }

    /// The number of vertices and triangles.
    pub(super) fn size(&self) -> Result<TessellationSize, TessellationError> {
        check_count(self.us.len(), self.vs.len())?;
        let vertices = self.us.len() * self.vs.len();
        let triangles = self
            .cells()
            .map(|(i, j)| self.cell_triangles(i, j).iter().flatten().count())
            .sum();
        Ok(TessellationSize {
            vertices,
            triangles,
        })
    }

    /// Writes the tessellation of `surface` into `buffers` and returns the
    /// deviation reached.
    pub(super) fn fill<S: ParametricSurface + ?Sized>(
        &self,
        surface: &S,
        tol: Tolerance,
        buffers: TessellationBuffers<'_>,
    ) -> Result<TessellationDeviation, TessellationError> {
        let size = self.size()?;
        let TessellationBuffers {
            vertices,
            normals,
            uvs,
            triangles,
        } = buffers;
        assert_eq!(vertices.len(), size.vertices, "vertices: wrong length");
        assert_eq!(normals.len(), size.vertices, "normals: wrong length");
        assert_eq!(uvs.len(), size.vertices, "uvs: wrong length");
        assert_eq!(triangles.len(), size.triangles, "triangles: wrong length");
        let [u_width, v_width] = reference_widths(surface, tol);
        for (j, &v) in self.vs.iter().enumerate() {
            for (i, &u) in self.us.iter().enumerate() {
                let k = self.index(i, j);
                let (point, [su, sv]) = surface.first_derivatives(u, v);
                vertices[k] = point;
                uvs[k] = Point2D::new(u, v);
                normals[k] = match direct_normal(su, sv, u_width, v_width, tol) {
                    Some(normal) => normal,
                    None => self.recovered_normal(surface, i, j, tol)?,
                };
            }
        }
        let columns = self.us.len();
        let place = |vertices: &mut [Point], k: usize, point: Point| vertices[k] = point;
        for j in (0..self.vs.len()).filter(|&j| self.collapsed_rows[j]) {
            let point = vertices[self.index(0, j)];
            for i in 1..columns {
                place(vertices, self.index(i, j), point);
            }
        }
        for i in (0..columns).filter(|&i| self.collapsed_columns[i]) {
            let point = vertices[self.index(i, 0)];
            for j in 1..self.vs.len() {
                place(vertices, self.index(i, j), point);
            }
        }
        if self.u_closed {
            for j in 0..self.vs.len() {
                let (first, last) = (self.index(0, j), self.index(columns - 1, j));
                place(vertices, last, vertices[first]);
                normals[last] = normals[first];
            }
        }
        if self.v_closed {
            let last_row = self.vs.len() - 1;
            for i in 0..columns {
                let (first, last) = (self.index(i, 0), self.index(i, last_row));
                place(vertices, last, vertices[first]);
                normals[last] = normals[first];
            }
        }
        let mut written = 0;
        for (i, j) in self.cells() {
            for triangle in self.cell_triangles(i, j).into_iter().flatten() {
                // The size check bounds every index by `u32::MAX`.
                triangles[written] = triangle.map(|index| index as u32);
                written += 1;
            }
        }
        Ok(TessellationDeviation {
            deviation: self.deviation,
            proven: self.proven,
        })
    }

    /// The normal at vertex `(i, j)` where its partials give none: the
    /// surface's own recovery, else the normal beside the vertex, halfway
    /// into a cell it belongs to along `v`, then along `u`, then at the
    /// cell's middle. A cone's apex has no normal of its own, and takes
    /// its meridian's on the side of its cells.
    fn recovered_normal<S: ParametricSurface + ?Sized>(
        &self,
        surface: &S,
        i: usize,
        j: usize,
        tol: Tolerance,
    ) -> Result<Direction, TessellationError> {
        let (u, v) = (self.us[i], self.vs[j]);
        let (beside_u, beside_v) = (beside(&self.us, i), beside(&self.vs, j));
        [(u, v), (u, beside_v), (beside_u, v), (beside_u, beside_v)]
            .into_iter()
            .find_map(|(u, v)| surface.normal(u, v, tol))
            .ok_or(TessellationError::UndefinedNormal { u, v })
    }
}

/// The parameter halfway from `values[k]` into a cell it bounds: toward
/// the next value when that one is larger, else toward the previous one;
/// `values[k]` itself when neither differs.
fn beside(values: &[f64], k: usize) -> f64 {
    let value = values[k];
    match (values.get(k + 1), k.checked_sub(1).map(|p| values[p])) {
        (Some(&next), _) if next > value => 0.5 * (value + next),
        (_, Some(previous)) if previous < value => 0.5 * (value + previous),
        _ => value,
    }
}
