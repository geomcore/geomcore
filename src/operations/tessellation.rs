//! Tessellation: triangles on a surface over a rectangle of its
//! parameters, every point of the surface there within a deviation of
//! the triangles and every point of the triangles within it of the
//! surface (see `CONTEXT.md`, ADR 0056). The trait, its criteria, buffers
//! and results; the grid each family plans is in the files below.

mod bspline;
mod grid;
mod revolution;
mod swept;

use grid::Grid;

use crate::geometry::surfaces::{ParametricSurface, check_range};
use crate::{
    BSplineSurface, BezierSurface, Cone, Cylinder, Direction, ElementarySurface, ExtrusionSurface,
    Plane, Point, Point2D, RevolutionSurface, Sphere, Surface, Tolerance, Torus, TrimmedSurface,
    TrimmedSurfaceConstructionError,
};
use std::fmt;

/// Error returned when a surface cannot be tessellated as asked.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum TessellationError {
    /// The maximum deviation is not positive and finite.
    InvalidDeviation,
    /// The surface is unbounded in `u` (a plane, an extrusion of a line)
    /// and no `u` range was given.
    UnboundedU,
    /// The surface is unbounded in `v` (a plane, a cylinder or a cone
    /// along its axis, an extrusion) and no `v` range was given.
    UnboundedV,
    /// The `u` range is not a range of the surface, for the reason given:
    /// as [`TrimmedSurface::new`] checks it.
    InvalidURange(TrimmedSurfaceConstructionError),
    /// The `v` range is not a range of the surface, for the reason given.
    InvalidVRange(TrimmedSurfaceConstructionError),
    /// The tessellation has more vertices than a `u32` triangle index
    /// reaches: the deviation is far below the surface's size.
    TooManyVertices,
    /// No normal is defined at the vertex at these parameters nor beside
    /// it: the surface there is smaller than the confusion tolerance.
    UndefinedNormal {
        /// The vertex's `u`.
        u: f64,
        /// The vertex's `v`.
        v: f64,
    },
}

impl fmt::Display for TessellationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TessellationError::InvalidDeviation => {
                f.write_str("the maximum deviation must be positive and finite")
            }
            TessellationError::UnboundedU => {
                f.write_str("the surface is unbounded in u: give a u range")
            }
            TessellationError::UnboundedV => {
                f.write_str("the surface is unbounded in v: give a v range")
            }
            TessellationError::InvalidURange(reason) => write!(f, "invalid u range: {reason}"),
            TessellationError::InvalidVRange(reason) => write!(f, "invalid v range: {reason}"),
            TessellationError::TooManyVertices => {
                f.write_str("the tessellation has more vertices than a u32 index reaches")
            }
            TessellationError::UndefinedNormal { u, v } => {
                write!(f, "no normal is defined at or beside ({u}, {v})")
            }
        }
    }
}

impl std::error::Error for TessellationError {}

/// What a tessellation must meet, and where: a maximum deviation, and the
/// parameter ranges to cover. A direction without a range covers the
/// surface's bounds, which must then be finite.
///
/// # Examples
///
/// ```
/// use geomcore::TessellationCriteria;
///
/// let criteria = TessellationCriteria::new(1e-3)?.with_v_range((0.0, 2.0));
/// assert_eq!(criteria.max_deviation(), 1e-3);
/// assert_eq!(criteria.u_range(), None);
/// assert_eq!(criteria.v_range(), Some((0.0, 2.0)));
/// # Ok::<(), geomcore::TessellationError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TessellationCriteria {
    max_deviation: f64,
    u_range: Option<(f64, f64)>,
    v_range: Option<(f64, f64)>,
}

impl TessellationCriteria {
    /// Criteria asking for a deviation of at most `max_deviation` over the
    /// surface's bounds.
    ///
    /// Any positive deviation is accepted, but the vertices cannot lie
    /// closer to the surface than their coordinates' rounding, about 1e-15
    /// times the surface's distance from the origin. A finer request only
    /// adds triangles, about ten times as many for each tenfold on a
    /// curved surface;
    /// [`tessellation_size`](Tessellate::tessellation_size) gives the
    /// counts before anything is allocated.
    ///
    /// # Errors
    ///
    /// [`TessellationError::InvalidDeviation`] unless `max_deviation` is
    /// positive and finite.
    pub fn new(max_deviation: f64) -> Result<Self, TessellationError> {
        if !(max_deviation.is_finite() && max_deviation > 0.0) {
            return Err(TessellationError::InvalidDeviation);
        }
        Ok(TessellationCriteria {
            max_deviation,
            u_range: None,
            v_range: None,
        })
    }

    /// The same criteria over the `u` range `(first, last)`, checked when
    /// a surface is tessellated: a periodic direction takes any start, so a
    /// range may cross its seam, but no more than one period.
    pub fn with_u_range(self, u_range: (f64, f64)) -> Self {
        TessellationCriteria {
            u_range: Some(u_range),
            ..self
        }
    }

    /// The same criteria over the `v` range `(first, last)`, as
    /// [`with_u_range`](Self::with_u_range).
    pub fn with_v_range(self, v_range: (f64, f64)) -> Self {
        TessellationCriteria {
            v_range: Some(v_range),
            ..self
        }
    }

    /// The largest deviation the tessellation may have.
    pub fn max_deviation(&self) -> f64 {
        self.max_deviation
    }

    /// The `u` range asked for, if any.
    pub fn u_range(&self) -> Option<(f64, f64)> {
        self.u_range
    }

    /// The `v` range asked for, if any.
    pub fn v_range(&self) -> Option<(f64, f64)> {
        self.v_range
    }
}

/// How many vertices and triangles a tessellation has: the lengths of the
/// buffers [`tessellate_into`](Tessellate::tessellate_into) fills.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TessellationSize {
    /// The number of vertices, normals and UV pairs.
    pub vertices: usize,
    /// The number of triangles.
    pub triangles: usize,
}

/// The caller's slices a tessellation is written into, each as long as
/// [`tessellation_size`](Tessellate::tessellation_size) says.
#[derive(Debug)]
pub struct TessellationBuffers<'a> {
    /// The vertices: points on the surface.
    pub vertices: &'a mut [Point],
    /// The surface's unit normal at each vertex.
    pub normals: &'a mut [Direction],
    /// The surface parameters `(u, v)` of each vertex, as `x` and `y`.
    pub uvs: &'a mut [Point2D],
    /// The triangles, three vertex indices each, counterclockwise seen
    /// from the side the normals point to.
    pub triangles: &'a mut [[u32; 3]],
}

/// The deviation a tessellation reached: every point of the surface over
/// the ranges lies within `deviation` of the triangles, and every point of
/// the triangles within it of the surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TessellationDeviation {
    /// The deviation, at most the one asked for.
    pub deviation: f64,
    /// Whether `deviation` is a proven bound. It is, but on a swept
    /// surface whose profile has no proven discretization (a quartic
    /// branch), where it is measured on samples.
    pub proven: bool,
}

/// A tessellated surface: vertices on it, joined by triangles. See
/// [`Tessellate`].
#[derive(Debug, Clone, PartialEq, Default)]
#[non_exhaustive]
pub struct Tessellation {
    /// The vertices: points on the surface.
    pub vertices: Vec<Point>,
    /// The surface's unit normal at each vertex.
    pub normals: Vec<Direction>,
    /// The surface parameters `(u, v)` of each vertex, as `x` and `y`.
    pub uvs: Vec<Point2D>,
    /// The triangles, three vertex indices each, counterclockwise seen
    /// from the side the normals point to.
    pub triangles: Vec<[u32; 3]>,
    /// The deviation reached (see [`TessellationDeviation`]).
    pub deviation: f64,
    /// Whether [`deviation`](Self::deviation) is a proven bound.
    pub proven: bool,
}

impl Tessellation {
    /// The tessellation made of its fields: for rebuilding a stored one,
    /// as Python's `pickle` does, since the struct is `non_exhaustive`.
    /// The fields are kept as they are, not checked.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction, Point, Point2D, Tessellation};
    /// let rebuilt = Tessellation::from_parts(
    ///     vec![Point::ORIGIN; 3],
    ///     vec![Direction::Z; 3],
    ///     vec![Point2D::ORIGIN; 3],
    ///     vec![[0, 1, 2]],
    ///     0.0,
    ///     true,
    /// );
    /// assert_eq!(rebuilt.triangles, [[0, 1, 2]]);
    /// ```
    pub fn from_parts(
        vertices: Vec<Point>,
        normals: Vec<Direction>,
        uvs: Vec<Point2D>,
        triangles: Vec<[u32; 3]>,
        deviation: f64,
        proven: bool,
    ) -> Self {
        Tessellation {
            vertices,
            normals,
            uvs,
            triangles,
            deviation,
            proven,
        }
    }
}

/// Tessellation of a surface over a rectangle of its parameters: a grid
/// of vertices on the surface, each cell of the grid two triangles, every
/// point of the surface within the deviation of the triangles and every
/// point of the triangles within it of the surface.
///
/// The grid's spacing comes from a bound on the deviation, worked out
/// before any vertex is computed (ADR 0056):
///
/// - a plane is one cell;
/// - a cylinder, cone, sphere and torus take equal steps around the axis
///   and along the meridian, their count from the chordal deviation of a
///   circle (its sagitta, `r (1 − cos(θ/2))` for an arc of angle `θ`);
/// - an extrusion takes the points of its profile's discretization across
///   and one step along the sweep, a surface of revolution the points of
///   its profile's discretization along it and equal steps around the
///   axis;
/// - a B-spline or Bézier surface takes equal steps in each knot span,
///   their count from the second differences of the span's Bézier control
///   points, as a B-spline curve's discretization does (ADR 0003).
///
/// The deviation is proven, but on a swept surface whose profile is a
/// quartic branch, whose discretization is measured on samples. Where a
/// whole row of the grid collapses to a point (a sphere's pole, a cone's
/// apex), its vertices are made one point and the triangles that would
/// have no area are left out; a cone's apex inside the `v` range splits
/// the grid there into two rows of vertices, one for each nappe's
/// normals. A range spanning a whole period gives a last row or column of
/// vertices at the same points as the first, with their own parameters.
///
/// The work is done in two steps, so a caller controls memory and threads
/// (ADR 0004): [`tessellation_size`](Self::tessellation_size) counts the
/// vertices and triangles, and [`tessellate_into`](Self::tessellate_into)
/// writes exactly that many into slices the caller owns.
/// [`tessellate`](Self::tessellate) wraps both and allocates the result.
/// The same surface, criteria and tolerance always give the same
/// tessellation. The tolerance decides which rows collapse (within
/// `tol.confusion` of a point) and the normals, as
/// [`normal`](ParametricSurface::normal) gives them.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{
///     Cylinder, Direction, Point, Point2D, TessellationBuffers, TessellationCriteria, Tolerance,
///     Vector,
/// };
///
/// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
/// // The axis is unbounded: the v range is needed.
/// let criteria = TessellationCriteria::new(1e-3).unwrap().with_v_range((0.0, 2.0));
/// let tol = Tolerance::DEFAULT;
///
/// // Count, then fill buffers the caller owns...
/// let size = cylinder.tessellation_size(criteria, tol).unwrap();
/// let mut vertices = vec![Point::ORIGIN; size.vertices];
/// let mut normals = vec![Direction::Z; size.vertices];
/// let mut uvs = vec![Point2D::ORIGIN; size.vertices];
/// let mut triangles = vec![[0; 3]; size.triangles];
/// let reached = cylinder
///     .tessellate_into(
///         criteria,
///         tol,
///         TessellationBuffers {
///             vertices: &mut vertices,
///             normals: &mut normals,
///             uvs: &mut uvs,
///             triangles: &mut triangles,
///         },
///     )
///     .unwrap();
/// assert!(reached.deviation <= 1e-3 && reached.proven);
///
/// // ...or let geomcore allocate the result.
/// let tessellation = cylinder.tessellate(criteria, tol).unwrap();
/// assert_eq!(tessellation.vertices, vertices);
/// ```
pub trait Tessellate {
    /// The number of vertices and triangles of the tessellation.
    ///
    /// # Errors
    ///
    /// [`TessellationError::UnboundedU`] or
    /// [`UnboundedV`](TessellationError::UnboundedV) for an unbounded
    /// direction without a range,
    /// [`InvalidURange`](TessellationError::InvalidURange) or
    /// [`InvalidVRange`](TessellationError::InvalidVRange) for a range
    /// that is not one of the surface, and
    /// [`TooManyVertices`](TessellationError::TooManyVertices).
    fn tessellation_size(
        &self,
        criteria: TessellationCriteria,
        tol: Tolerance,
    ) -> Result<TessellationSize, TessellationError>;

    /// Writes the tessellation into caller-owned slices and returns the
    /// deviation reached. It allocates only the grid's parameters, a few
    /// numbers per row and column, never per vertex.
    ///
    /// # Errors
    ///
    /// As [`tessellation_size`](Self::tessellation_size), and
    /// [`TessellationError::UndefinedNormal`] where no normal is defined
    /// at a vertex nor beside it, the buffers then partly written.
    ///
    /// # Panics
    ///
    /// Panics if a slice has the wrong length.
    fn tessellate_into(
        &self,
        criteria: TessellationCriteria,
        tol: Tolerance,
        buffers: TessellationBuffers<'_>,
    ) -> Result<TessellationDeviation, TessellationError>;

    /// Counts, allocates and fills the tessellation.
    ///
    /// # Errors
    ///
    /// As [`tessellate_into`](Self::tessellate_into).
    fn tessellate(
        &self,
        criteria: TessellationCriteria,
        tol: Tolerance,
    ) -> Result<Tessellation, TessellationError>;
}

/// The range of one direction: the one asked for, checked against the
/// direction's bounds and period, or the bounds when finite. A periodic
/// range a rounding longer than the period is cut to one period.
fn resolve_range(
    requested: Option<(f64, f64)>,
    bounds: (f64, f64),
    period: Option<f64>,
    unbounded: TessellationError,
    invalid: fn(TrimmedSurfaceConstructionError) -> TessellationError,
) -> Result<(f64, f64), TessellationError> {
    let (first, last) = match requested {
        None if bounds.0.is_finite() && bounds.1.is_finite() => return Ok(bounds),
        None => return Err(unbounded),
        Some(range) => range,
    };
    check_range((first, last), bounds, period).map_err(invalid)?;
    Ok(match period {
        Some(period) if last - first > period => (first, first + period),
        _ => (first, last),
    })
}

/// The ranges `criteria` asks for on a surface with these bounds and
/// periods.
fn resolve_ranges<S: ParametricSurface + ?Sized>(
    surface: &S,
    criteria: &TessellationCriteria,
) -> Result<[(f64, f64); 2], TessellationError> {
    Ok([
        resolve_range(
            criteria.u_range,
            surface.u_bounds(),
            surface.u_period(),
            TessellationError::UnboundedU,
            TessellationError::InvalidURange,
        )?,
        resolve_range(
            criteria.v_range,
            surface.v_bounds(),
            surface.v_period(),
            TessellationError::UnboundedV,
            TessellationError::InvalidVRange,
        )?,
    ])
}

/// The grid of `surface` over `ranges`, from the routing hooks: every
/// surface answers one of them.
fn plan<S: ParametricSurface + ?Sized>(
    surface: &S,
    ranges: [(f64, f64); 2],
    max_deviation: f64,
    tol: Tolerance,
) -> Result<Grid, TessellationError> {
    if let Some(bspline) = surface.as_bspline() {
        bspline::grid(bspline, ranges, max_deviation, tol)
    } else if let Some(elementary) = surface.as_elementary() {
        revolution::elementary_grid(&elementary, ranges, max_deviation, tol)
    } else if let Some(extrusion) = surface.as_extrusion() {
        swept::extrusion_grid(extrusion, ranges, max_deviation)
    } else if let Some(revolution) = surface.as_revolution() {
        swept::revolution_grid(revolution, ranges, max_deviation, tol)
    } else {
        unreachable!("every tessellated surface type answers one routing hook")
    }
}

/// Counts, allocates and fills the tessellation of `surface` over `grid`.
fn allocated<S: ParametricSurface + ?Sized>(
    surface: &S,
    grid: &Grid,
    tol: Tolerance,
) -> Result<Tessellation, TessellationError> {
    let size = grid.size()?;
    let mut tessellation = Tessellation {
        vertices: vec![Point::ORIGIN; size.vertices],
        normals: vec![Direction::Z; size.vertices],
        uvs: vec![Point2D::ORIGIN; size.vertices],
        triangles: vec![[0; 3]; size.triangles],
        deviation: 0.0,
        proven: true,
    };
    let reached = grid.fill(
        surface,
        tol,
        TessellationBuffers {
            vertices: &mut tessellation.vertices,
            normals: &mut tessellation.normals,
            uvs: &mut tessellation.uvs,
            triangles: &mut tessellation.triangles,
        },
    )?;
    tessellation.deviation = reached.deviation;
    tessellation.proven = reached.proven;
    Ok(tessellation)
}

/// Implements [`Tessellate`] for whole surfaces, over their own bounds.
macro_rules! tessellate_surfaces {
    ($($surface:ty),*) => {$(
        impl Tessellate for $surface {
            fn tessellation_size(
                &self,
                criteria: TessellationCriteria,
                tol: Tolerance,
            ) -> Result<TessellationSize, TessellationError> {
                let ranges = resolve_ranges(self, &criteria)?;
                plan(self, ranges, criteria.max_deviation, tol)?.size()
            }

            fn tessellate_into(
                &self,
                criteria: TessellationCriteria,
                tol: Tolerance,
                buffers: TessellationBuffers<'_>,
            ) -> Result<TessellationDeviation, TessellationError> {
                let ranges = resolve_ranges(self, &criteria)?;
                plan(self, ranges, criteria.max_deviation, tol)?.fill(self, tol, buffers)
            }

            fn tessellate(
                &self,
                criteria: TessellationCriteria,
                tol: Tolerance,
            ) -> Result<Tessellation, TessellationError> {
                let ranges = resolve_ranges(self, &criteria)?;
                allocated(self, &plan(self, ranges, criteria.max_deviation, tol)?, tol)
            }
        }
    )*};
}

tessellate_surfaces!(
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    ExtrusionSurface,
    RevolutionSurface,
    BSplineSurface,
    BezierSurface,
    Surface,
    ElementarySurface
);

/// Over its own ranges by default; a range asked for must lie inside
/// them. The basis plans the grid.
impl<S: ParametricSurface> Tessellate for TrimmedSurface<S> {
    fn tessellation_size(
        &self,
        criteria: TessellationCriteria,
        tol: Tolerance,
    ) -> Result<TessellationSize, TessellationError> {
        let ranges = resolve_ranges(self, &criteria)?;
        plan(self.basis(), ranges, criteria.max_deviation, tol)?.size()
    }

    fn tessellate_into(
        &self,
        criteria: TessellationCriteria,
        tol: Tolerance,
        buffers: TessellationBuffers<'_>,
    ) -> Result<TessellationDeviation, TessellationError> {
        let ranges = resolve_ranges(self, &criteria)?;
        plan(self.basis(), ranges, criteria.max_deviation, tol)?.fill(self.basis(), tol, buffers)
    }

    fn tessellate(
        &self,
        criteria: TessellationCriteria,
        tol: Tolerance,
    ) -> Result<Tessellation, TessellationError> {
        let ranges = resolve_ranges(self, &criteria)?;
        let grid = plan(self.basis(), ranges, criteria.max_deviation, tol)?;
        allocated(self.basis(), &grid, tol)
    }
}
