//! B-spline (NURBS) surfaces in 3D: tensor-product de Boor evaluation with
//! first derivatives, rational (weighted) poles, and periodic surfaces
//! in either direction, a thin wrapper over [`crate::bspline::surface`].
//!
//! Construction stores the pole grid row-major, with a separate weight per
//! pole if rational, and each direction's knot vector with its flattened
//! knot sequence, so evaluation never re-derives them.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::bspline::curve;
use crate::bspline::editing::{self, Edit, Elevation, Insertion, Lanes};
use crate::bspline::surface as math;
use crate::curves::{Continuity, Side};
use crate::math::boxtree::{BoxTree, Lazy, box_distance};
use crate::math::newton_min_2d;
use crate::projection::SurfaceProjection;
use crate::surfaces::iso::{SurfaceIsoError, check_iso_param};
use crate::surfaces::{BezierSurface, ParametricSurface};
use crate::{BSplineCurve, Point, Tolerance, Transform, Vector};

pub use crate::curves::BSplineConstructionError;
use crate::curves::{
    BSplineSplitError, DegreeElevationError, KnotInsertionError, KnotRemovalError,
};

/// A parametric direction of a surface.
#[derive(Clone, Copy)]
enum Direction {
    U,
    V,
}

/// A tensor-product B-spline (optionally rational, optionally periodic in
/// either direction) surface in 3D, evaluated by de Boor's algorithm applied
/// once per parametric direction.
///
/// The pole grid `poles[i_u][i_v]` defines the control net; `i_u` indexes the
/// `u` direction and `i_v` the `v` direction. Each direction carries its own
/// degree, `knots`/`multiplicities`, and periodic flag, validated exactly as a
/// [`crate::BSplineCurve`] direction is. A rational surface additionally
/// carries a `weight` per pole and is evaluated in homogeneous coordinates,
/// then projected back to Euclidean space via the quotient rule.
///
/// Derivatives up to second order are supported; see
/// [`BSplineSurface::eval_derivative`].
///
/// # Examples
///
/// A bilinear (degree-1 × degree-1) patch over the unit square:
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::surfaces::BSplineSurface;
/// use geomcore::Point;
///
/// let poles = vec![
///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
/// ];
/// let surface = BSplineSurface::new(
///     1,
///     1,
///     poles,
///     vec![0.0, 1.0],
///     vec![2, 2],
///     vec![0.0, 1.0],
///     vec![2, 2],
///     false,
///     false,
/// )
/// .unwrap();
///
/// assert_eq!(surface.eval_point(0.5, 0.5), Point::new(0.5, 0.5, 0.25));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BSplineSurface {
    u_knots: curve::KnotVector,
    v_knots: curve::KnotVector,
    /// The pole grid, row-major `[i_u][i_v]`: a row per `u` pole.
    poles: Vec<Point>,
    /// One weight per pole, laid out as `poles`, if rational.
    weights: Option<Vec<f64>>,
    /// Bounding-box hierarchy of the knot-span patches, built on the first
    /// projection.
    patches: Lazy<PatchTree>,
}

/// The knot-span boundaries of each direction, whether the surface is a
/// crease or ends there ([`curve::KnotVector::crease_breaks`]), each
/// span's Bezier matrix ([`curve::KnotVector::span_bezier_matrices`]), and
/// the box hierarchy over the patches (`i_u * (v_breaks.len() - 1) + i_v`).
#[derive(Clone)]
struct PatchTree {
    u_breaks: Vec<f64>,
    v_breaks: Vec<f64>,
    u_creases: Vec<bool>,
    v_creases: Vec<bool>,
    u_spans: (Vec<usize>, Vec<f64>),
    v_spans: (Vec<usize>, Vec<f64>),
    tree: BoxTree,
}

impl PatchTree {
    /// The parameter rectangle `[u0, u1, v0, v1]` of patch `k` and whether
    /// each of its edges is a crease or a boundary.
    fn patch(&self, k: usize) -> ([f64; 4], [bool; 4]) {
        let n_vspans = self.v_breaks.len() - 1;
        let (a, b) = (k / n_vspans, k % n_vspans);
        (
            [
                self.u_breaks[a],
                self.u_breaks[a + 1],
                self.v_breaks[b],
                self.v_breaks[b + 1],
            ],
            [
                self.u_creases[a],
                self.u_creases[a + 1],
                self.v_creases[b],
                self.v_creases[b + 1],
            ],
        )
    }

    /// The number of patches.
    fn len(&self) -> usize {
        (self.u_breaks.len() - 1) * (self.v_breaks.len() - 1)
    }

    /// Distances below this much of the best one found are not worth
    /// telling apart: a relative `1e-12` of the surface's extent and of
    /// its distance to `point`.
    fn resolution(&self, point: Point) -> f64 {
        let b = self.tree.bounds();
        let diagonal =
            ((b[3] - b[0]).powi(2) + (b[4] - b[1]).powi(2) + (b[5] - b[2]).powi(2)).sqrt();
        1e-12 * (diagonal + box_distance(point, b))
    }
}

impl BSplineSurface {
    /// Creates a non-rational (polynomial) tensor-product B-spline surface.
    ///
    /// `poles[i_u][i_v]` is the control net: `i_u` indexes the `u` direction
    /// (`poles.len()` rows) and `i_v` the `v` direction (row length).
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConstructionError::PoleCountMismatch`] if the pole grid
    /// is empty or ragged (rows of differing length), since a well-formed
    /// tensor-product net must be rectangular. Otherwise returns a
    /// [`BSplineConstructionError`] if either direction's `degree`, `knots`,
    /// `multiplicities`, and pole count (grid rows for `u`, columns for `v`)
    /// are not mutually consistent (see the error variants for the exact
    /// conditions checked, applied per direction).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::surfaces::BSplineSurface;
    /// use geomcore::Point;
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// let surface = BSplineSurface::new(
    ///     1, 1, poles,
    ///     vec![0.0, 1.0], vec![2, 2],
    ///     vec![0.0, 1.0], vec![2, 2],
    ///     false, false,
    /// )
    /// .unwrap();
    /// assert_eq!(surface.eval_point(0.0, 0.0), Point::new(0.0, 0.0, 0.0));
    /// ```
    #[allow(clippy::too_many_arguments)] // mirrors the two-direction knot/degree data model
    pub fn new(
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Vec<Point>>,
        u_knots: Vec<f64>,
        u_multiplicities: Vec<u32>,
        v_knots: Vec<f64>,
        v_multiplicities: Vec<u32>,
        u_periodic: bool,
        v_periodic: bool,
    ) -> Result<BSplineSurface, BSplineConstructionError> {
        Self::build(
            u_degree,
            v_degree,
            poles,
            None,
            u_knots,
            u_multiplicities,
            v_knots,
            v_multiplicities,
            u_periodic,
            v_periodic,
        )
    }

    /// Creates a rational (NURBS) tensor-product B-spline surface with one
    /// weight per pole, laid out `weights[i_u][i_v]` matching `poles`.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConstructionError::WeightCountMismatch`] if the weight
    /// grid does not match the pole grid shape, or
    /// [`BSplineConstructionError::NonPositiveWeight`] if any weight is `<= 0`.
    /// Otherwise validates the pole grid and both directions exactly as
    /// [`BSplineSurface::new`] does.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::surfaces::BSplineSurface;
    /// use geomcore::Point;
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// let weights = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
    /// let surface = BSplineSurface::new_rational(
    ///     1, 1, poles, weights,
    ///     vec![0.0, 1.0], vec![2, 2],
    ///     vec![0.0, 1.0], vec![2, 2],
    ///     false, false,
    /// )
    /// .unwrap();
    /// // All-unit weights reproduce the polynomial patch.
    /// assert_eq!(surface.eval_point(0.5, 0.5), Point::new(0.5, 0.5, 0.25));
    /// ```
    #[allow(clippy::too_many_arguments)] // mirrors the two-direction knot/degree data model
    pub fn new_rational(
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Vec<Point>>,
        weights: Vec<Vec<f64>>,
        u_knots: Vec<f64>,
        u_multiplicities: Vec<u32>,
        v_knots: Vec<f64>,
        v_multiplicities: Vec<u32>,
        u_periodic: bool,
        v_periodic: bool,
    ) -> Result<BSplineSurface, BSplineConstructionError> {
        Self::build(
            u_degree,
            v_degree,
            poles,
            Some(weights),
            u_knots,
            u_multiplicities,
            v_knots,
            v_multiplicities,
            u_periodic,
            v_periodic,
        )
    }

    /// A polynomial surface from flat knots consistent with its pole grid
    /// (row-major, a row per `u` pole), clamped in both directions:
    /// fitting results skip revalidation.
    pub(crate) fn from_flat_knots(
        (u_degree, u_flat): (usize, Vec<f64>),
        (v_degree, v_flat): (usize, Vec<f64>),
        poles: Vec<Point>,
    ) -> BSplineSurface {
        let u_knots = curve::KnotVector::from_flat(u_degree, u_flat, false);
        let v_knots = curve::KnotVector::from_flat(v_degree, v_flat, false);
        debug_assert_eq!(u_knots.pole_count() * v_knots.pole_count(), poles.len());
        BSplineSurface {
            u_knots,
            v_knots,
            poles,
            weights: None,
            patches: Lazy::default(),
        }
    }

    #[allow(clippy::too_many_arguments)] // internal shared builder
    fn build(
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Vec<Point>>,
        weights: Option<Vec<Vec<f64>>>,
        u_knots: Vec<f64>,
        u_multiplicities: Vec<u32>,
        v_knots: Vec<f64>,
        v_multiplicities: Vec<u32>,
        u_periodic: bool,
        v_periodic: bool,
    ) -> Result<BSplineSurface, BSplineConstructionError> {
        // Rectangular, non-empty grid. A ragged or empty net is a shape error,
        // reported as PoleCountMismatch.
        let n_u = poles.len();
        if n_u == 0 || poles[0].is_empty() {
            return Err(BSplineConstructionError::PoleCountMismatch);
        }
        let n_v = poles[0].len();
        if poles.iter().any(|row| row.len() != n_v) {
            return Err(BSplineConstructionError::PoleCountMismatch);
        }

        // Weight grid, if present, must match the pole grid shape, and be
        // strictly positive.
        if let Some(ws) = &weights {
            if ws.len() != n_u || ws.iter().any(|row| row.len() != n_v) {
                return Err(BSplineConstructionError::WeightCountMismatch);
            }
            if ws.iter().flatten().any(|&w| w <= 0.0) {
                return Err(BSplineConstructionError::NonPositiveWeight);
            }
        }

        // Validate each direction: u over grid rows, v over grid columns.
        let u_knots = curve::KnotVector::new(u_degree, u_knots, u_multiplicities, u_periodic, n_u)?;
        let v_knots = curve::KnotVector::new(v_degree, v_knots, v_multiplicities, v_periodic, n_v)?;

        Ok(BSplineSurface {
            u_knots,
            v_knots,
            poles: poles.concat(),
            weights: weights.map(|w| w.concat()),
            patches: Lazy::default(),
        })
    }

    /// Returns the surface's degree in the `u` direction.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert_eq!(surface.u_degree(), 1);
    /// ```
    pub fn u_degree(&self) -> usize {
        self.u_knots.degree()
    }

    /// Returns the surface's degree in the `v` direction.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert_eq!(surface.v_degree(), 1);
    /// ```
    pub fn v_degree(&self) -> usize {
        self.v_knots.degree()
    }

    /// Returns whether the surface is periodic in the `u` direction.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert!(!surface.is_u_periodic());
    /// ```
    pub fn is_u_periodic(&self) -> bool {
        self.u_knots.is_periodic()
    }

    /// Returns whether the surface is periodic in the `v` direction.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert!(!surface.is_v_periodic());
    /// ```
    pub fn is_v_periodic(&self) -> bool {
        self.v_knots.is_periodic()
    }

    /// Whether the isoparametric curves at the two bounds of `direction`
    /// coincide. Both share the knot vector across `direction`, so they
    /// are the same curve exactly when their poles match and their weights
    /// are proportional; that holds for clamped and unclamped knots alike.
    fn boundaries_meet(&self, direction: Direction, tol: Tolerance) -> bool {
        let (first, last) = match direction {
            Direction::U => self.u_knots.bounds(),
            Direction::V => self.v_knots.bounds(),
        };
        let start = self.iso(direction, first);
        let end = self.iso(direction, last);
        let poles_meet = start
            .poles()
            .iter()
            .zip(end.poles())
            .all(|(a, b)| a.distance(*b) <= tol.confusion);
        poles_meet
            && match (start.weights(), end.weights()) {
                (None, None) => true,
                (Some(a), Some(b)) => {
                    let scale = a[0] / b[0];
                    a.iter()
                        .zip(b)
                        .all(|(&wa, &wb)| (wa - wb * scale).abs() <= tol.confusion * wa)
                }
                // One boundary has uniform weights, the other does not.
                _ => false,
            }
    }

    /// Returns whether the surface is rational (has per-pole weights).
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert!(!surface.is_rational());
    /// ```
    pub fn is_rational(&self) -> bool {
        self.weights.is_some()
    }

    /// Returns the control net row by row: pole `(i_u, i_v)` at
    /// `i_u * v_pole_count() + i_v`, a row per `u` pole.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert_eq!(surface.poles()[surface.v_pole_count() + 1], Point::new(1.0, 1.0, 1.0));
    /// ```
    pub fn poles(&self) -> &[Point] {
        &self.poles
    }

    /// Returns the per-pole weights, laid out as [`BSplineSurface::poles`],
    /// or `None` if the surface is not rational.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert_eq!(surface.weights(), None);
    /// ```
    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    /// Returns the number of poles in the `u` direction: the rows of the net.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert_eq!(surface.u_pole_count(), 2);
    /// ```
    pub fn u_pole_count(&self) -> usize {
        self.u_knots.pole_count()
    }

    /// Returns the number of poles in the `v` direction: the poles of a row.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert_eq!(surface.v_pole_count(), 2);
    /// ```
    pub fn v_pole_count(&self) -> usize {
        self.v_knots.pole_count()
    }

    /// Returns the (non-flattened) knot values of the `u` direction.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert_eq!(surface.u_knots(), &[0.0, 1.0]);
    /// ```
    pub fn u_knots(&self) -> &[f64] {
        self.u_knots.knots()
    }

    /// Returns the knot multiplicities of the `u` direction.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert_eq!(surface.u_multiplicities(), &[2, 2]);
    /// ```
    pub fn u_multiplicities(&self) -> &[u32] {
        self.u_knots.multiplicities()
    }

    /// Returns the (non-flattened) knot values of the `v` direction.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert_eq!(surface.v_knots(), &[0.0, 1.0]);
    /// ```
    pub fn v_knots(&self) -> &[f64] {
        self.v_knots.knots()
    }

    /// Returns the knot multiplicities of the `v` direction.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::Point;
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// assert_eq!(surface.v_multiplicities(), &[2, 2]);
    /// ```
    pub fn v_multiplicities(&self) -> &[u32] {
        self.v_knots.multiplicities()
    }

    /// Returns the same surface with the knot `u` inserted `times` times in
    /// the `u` direction: the net gains `times` rows of poles, every point
    /// and derivative staying the same. Every column of the net is refined
    /// in the same pass.
    ///
    /// Rational and periodic directions behave as for
    /// [`crate::BSplineCurve::insert_knot`].
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::insert_knot`], against the `u` bounds,
    /// degree and multiplicities.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 1, poles, vec![0.0, 1.0], vec![3, 3], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let refined = surface.insert_u_knot(0.5, 2).unwrap();
    /// assert_eq!(refined.u_multiplicities(), &[3, 2, 3]);
    /// assert_eq!(refined.u_pole_count(), 5);
    /// let (p, q) = (surface.eval_point(0.3, 0.6), refined.eval_point(0.3, 0.6));
    /// assert!(p.distance(q) < 1e-12);
    /// ```
    pub fn insert_u_knot(&self, u: f64, times: u32) -> Result<BSplineSurface, KnotInsertionError> {
        let copies = [u; curve::MAX_DEGREE + 1];
        self.refine_u_knots(&copies[..(times as usize).min(copies.len())])
    }

    /// Returns the same surface with the knot `v` inserted `times` times in
    /// the `v` direction: every row of the net gains `times` poles, every
    /// point and derivative staying the same. Every row is refined in the
    /// same pass.
    ///
    /// Rational and periodic directions behave as for
    /// [`crate::BSplineCurve::insert_knot`].
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::insert_knot`], against the `v` bounds,
    /// degree and multiplicities.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..2)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     1, 2, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![3, 3], false, false,
    /// )
    /// .unwrap();
    /// let refined = surface.insert_v_knot(0.25, 1).unwrap();
    /// assert_eq!(refined.v_knots(), &[0.0, 0.25, 1.0]);
    /// assert_eq!(refined.v_pole_count(), 4);
    /// let (p, q) = (surface.eval_point(0.3, 0.6), refined.eval_point(0.3, 0.6));
    /// assert!(p.distance(q) < 1e-12);
    /// ```
    pub fn insert_v_knot(&self, v: f64, times: u32) -> Result<BSplineSurface, KnotInsertionError> {
        let copies = [v; curve::MAX_DEGREE + 1];
        self.refine_v_knots(&copies[..(times as usize).min(copies.len())])
    }

    /// Returns the same surface with all of `knots` inserted in the `u`
    /// direction in one pass, each value once, a repeated value as many
    /// times as it appears: the net gains one row of poles per value,
    /// every point and derivative staying the same.
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::refine_knots`], against the `u` bounds,
    /// degree and multiplicities.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineSurface, KnotInsertionError, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 1, poles, vec![0.0, 1.0], vec![3, 3], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let refined = surface.refine_u_knots(&[0.25, 0.5, 0.5]).unwrap();
    /// assert_eq!(refined.u_knots(), &[0.0, 0.25, 0.5, 1.0]);
    /// assert_eq!(refined.u_multiplicities(), &[3, 1, 2, 3]);
    /// assert_eq!(
    ///     surface.refine_u_knots(&[0.5, 1.5]),
    ///     Err(KnotInsertionError::KnotOutOfRange)
    /// );
    /// ```
    pub fn refine_u_knots(&self, knots: &[f64]) -> Result<BSplineSurface, KnotInsertionError> {
        if knots.is_empty() {
            return Ok(self.clone());
        }
        self.edited(Direction::U, &Insertion(knots))
    }

    /// Returns the same surface with all of `knots` inserted in the `v`
    /// direction in one pass, each value once, a repeated value as many
    /// times as it appears: every row of the net gains one pole per value,
    /// every point and derivative staying the same.
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::refine_knots`], against the `v` bounds,
    /// degree and multiplicities.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..2)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     1, 2, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![3, 3], false, false,
    /// )
    /// .unwrap();
    /// let refined = surface.refine_v_knots(&[0.2, 0.4, 0.6, 0.8]).unwrap();
    /// assert_eq!(refined.v_pole_count(), 7);
    /// let (p, q) = (surface.eval_point(0.5, 0.45), refined.eval_point(0.5, 0.45));
    /// assert!(p.distance(q) < 1e-12);
    /// ```
    pub fn refine_v_knots(&self, knots: &[f64]) -> Result<BSplineSurface, KnotInsertionError> {
        if knots.is_empty() {
            return Ok(self.clone());
        }
        self.edited(Direction::V, &Insertion(knots))
    }

    /// Returns the same surface with its `u` degree raised by `amount`: it
    /// evaluates to the same points and derivatives everywhere, up to
    /// rounding.
    ///
    /// Raises every column of the net as
    /// [`crate::BSplineCurve::elevate_degree`] raises a curve, all
    /// columns in one sweep over the `u` knot spans: every `u`
    /// multiplicity rises by `amount` (an unclamped non-periodic direction
    /// comes back clamped at its bounds), and the `v` direction is
    /// unchanged. An `amount` of 0 returns an identical copy.
    ///
    /// # Errors
    ///
    /// [`DegreeElevationError::DegreeTooHigh`] if the raised degree would
    /// exceed 25.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// let surface = BSplineSurface::new(
    ///     1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let raised = surface.elevate_u_degree(2).unwrap();
    /// assert_eq!((raised.u_degree(), raised.v_degree()), (3, 1));
    /// let (p, q) = (raised.eval_point(0.3, 0.7), surface.eval_point(0.3, 0.7));
    /// assert!(p.distance(q) < 1e-12);
    /// ```
    pub fn elevate_u_degree(&self, amount: usize) -> Result<BSplineSurface, DegreeElevationError> {
        self.edited(Direction::U, &Elevation(amount))
    }

    /// Returns the same surface with its `v` degree raised by `amount`: the
    /// counterpart of [`BSplineSurface::elevate_u_degree`] along `v`, every
    /// row of the net raised in one sweep over the `v` knot spans.
    ///
    /// # Errors
    ///
    /// [`DegreeElevationError::DegreeTooHigh`] if the raised degree would
    /// exceed 25.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// let surface = BSplineSurface::new(
    ///     1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let raised = surface.elevate_v_degree(1).unwrap();
    /// assert_eq!((raised.u_degree(), raised.v_degree()), (1, 2));
    /// let (p, q) = (raised.eval_point(0.3, 0.7), surface.eval_point(0.3, 0.7));
    /// assert!(p.distance(q) < 1e-12);
    /// ```
    pub fn elevate_v_degree(&self, amount: usize) -> Result<BSplineSurface, DegreeElevationError> {
        self.edited(Direction::V, &Elevation(amount))
    }

    /// Returns the surface with the knot `u` removed up to `times` times in
    /// the `u` direction, and how many times it was: each removal takes a
    /// row of poles away, as long as every point of the surface stays
    /// within `max_deviation` of the point of this one at the same
    /// parameters. A removal must keep every column of the net within the
    /// deviation allowed, so all columns lose the knot the same number of
    /// times.
    ///
    /// The deviation bound, rational and periodic directions behave as for
    /// [`crate::BSplineCurve::remove_knot`].
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::remove_knot`], against the `u` bounds
    /// and knots.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 1, poles, vec![0.0, 1.0], vec![3, 3], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let refined = surface.insert_u_knot(0.5, 2).unwrap();
    /// let (simpler, removed) = refined.remove_u_knot(0.5, 2, 1e-9).unwrap();
    /// assert_eq!(removed, 2);
    /// assert_eq!(simpler.u_knots(), surface.u_knots());
    /// assert_eq!(simpler.u_pole_count(), 3);
    /// ```
    pub fn remove_u_knot(
        &self,
        u: f64,
        times: u32,
        max_deviation: f64,
    ) -> Result<(BSplineSurface, u32), KnotRemovalError> {
        let copies = [u; curve::MAX_DEGREE + 1];
        let copies = &copies[..(times as usize).min(copies.len())];
        if copies.is_empty() {
            editing::check_removal(&self.u_knots, &[u], max_deviation)?;
        }
        let (surface, removed) = self.remove_u_knots(copies, max_deviation)?;
        Ok((surface, removed as u32))
    }

    /// Returns the surface with the knot `v` removed up to `times` times in
    /// the `v` direction, and how many times it was: each removal takes a
    /// pole away from every row, as long as every point of the surface
    /// stays within `max_deviation` of the point of this one at the same
    /// parameters. A removal must keep every row of the net within the
    /// deviation allowed, so all rows lose the knot the same number of
    /// times.
    ///
    /// The deviation bound, rational and periodic directions behave as for
    /// [`crate::BSplineCurve::remove_knot`].
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::remove_knot`], against the `v` bounds
    /// and knots.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..2)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     1, 2, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![3, 3], false, false,
    /// )
    /// .unwrap();
    /// let refined = surface.insert_v_knot(0.25, 1).unwrap();
    /// let (simpler, removed) = refined.remove_v_knot(0.25, 1, 1e-9).unwrap();
    /// assert_eq!(removed, 1);
    /// assert_eq!(simpler.v_knots(), surface.v_knots());
    ///
    /// // A ridge along `v = 1` stays.
    /// let poles: Vec<Vec<Point>> = (0..2)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (j % 2) as f64)).collect())
    ///     .collect();
    /// let ridge = BSplineSurface::new(
    ///     1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0, 2.0], vec![2, 1, 2], false, false,
    /// )
    /// .unwrap();
    /// let (same, removed) = ridge.remove_v_knot(1.0, 1, 1e-3).unwrap();
    /// assert_eq!(removed, 0);
    /// assert_eq!(same, ridge);
    /// ```
    pub fn remove_v_knot(
        &self,
        v: f64,
        times: u32,
        max_deviation: f64,
    ) -> Result<(BSplineSurface, u32), KnotRemovalError> {
        let copies = [v; curve::MAX_DEGREE + 1];
        let copies = &copies[..(times as usize).min(copies.len())];
        if copies.is_empty() {
            editing::check_removal(&self.v_knots, &[v], max_deviation)?;
        }
        let (surface, removed) = self.remove_v_knots(copies, max_deviation)?;
        Ok((surface, removed as u32))
    }

    /// Returns the surface with each of `knots` removed once in the `u`
    /// direction, in order, a value repeated as many times as it appears,
    /// and how many copies were removed in all; every removal keeps the
    /// surface within `max_deviation` of this one, however many are made.
    ///
    /// See [`crate::BSplineCurve::remove_knots`] for the order and
    /// [`BSplineSurface::remove_u_knot`] for the columns of the net.
    ///
    /// # Errors
    ///
    /// As [`BSplineSurface::remove_u_knot`], for any value of `knots`; no
    /// knot is removed then.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 1, poles, vec![0.0, 1.0], vec![3, 3], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let refined = surface.refine_u_knots(&[0.25, 0.5, 0.5]).unwrap();
    /// let (simpler, removed) = refined.remove_u_knots(&[0.5, 0.25, 0.5], 1e-9).unwrap();
    /// assert_eq!(removed, 3);
    /// assert_eq!(simpler.u_knots(), surface.u_knots());
    /// ```
    pub fn remove_u_knots(
        &self,
        knots: &[f64],
        max_deviation: f64,
    ) -> Result<(BSplineSurface, usize), KnotRemovalError> {
        self.removed(Direction::U, knots, max_deviation)
    }

    /// Returns the surface with each of `knots` removed once in the `v`
    /// direction, in order, a value repeated as many times as it appears,
    /// and how many copies were removed in all; every removal keeps the
    /// surface within `max_deviation` of this one, however many are made.
    ///
    /// See [`crate::BSplineCurve::remove_knots`] for the order and
    /// [`BSplineSurface::remove_v_knot`] for the rows of the net.
    ///
    /// # Errors
    ///
    /// As [`BSplineSurface::remove_v_knot`], for any value of `knots`; no
    /// knot is removed then.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineSurface, KnotRemovalError, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..2)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     1, 2, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![3, 3], false, false,
    /// )
    /// .unwrap();
    /// let refined = surface.refine_v_knots(&[0.2, 0.4, 0.6, 0.8]).unwrap();
    /// let (simpler, removed) = refined.remove_v_knots(&[0.2, 0.4, 0.6, 0.8], 1e-9).unwrap();
    /// assert_eq!(removed, 4);
    /// assert_eq!(simpler.v_pole_count(), 3);
    /// assert_eq!(
    ///     refined.remove_v_knots(&[0.2, 1.0], 1e-9),
    ///     Err(KnotRemovalError::KnotOutOfRange)
    /// );
    /// ```
    pub fn remove_v_knots(
        &self,
        knots: &[f64],
        max_deviation: f64,
    ) -> Result<(BSplineSurface, usize), KnotRemovalError> {
        self.removed(Direction::V, knots, max_deviation)
    }

    /// Returns the polynomial (or rational) patch of each pair of knot
    /// spans within the bounds as a [`BezierSurface`]: `patches[i][j]`
    /// covers the `i`-th `u` span and the `j`-th `v` span, and keeps them
    /// as its bounds.
    ///
    /// The net is cut into strips of `u` spans, then each strip into
    /// patches of `v` spans, as [`crate::BSplineCurve::to_bezier_pieces`]
    /// cuts a curve. A periodic direction is decomposed from its lower
    /// bound. [`BezierSurface::to_bspline`] gives a patch as a B-spline
    /// surface.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..4)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 2, poles, vec![0.0, 1.0, 2.0], vec![3, 1, 3], vec![0.0, 1.0], vec![3, 3], false, false,
    /// )
    /// .unwrap();
    /// let patches = surface.to_bezier_patches();
    /// assert_eq!((patches.len(), patches[0].len()), (2, 1));
    /// assert_eq!(patches[1][0].u_bounds(), (1.0, 2.0));
    /// let (p, q) = (surface.eval_point(1.5, 0.5), patches[1][0].eval_point(1.5, 0.5));
    /// assert!(p.distance(q) < 1e-12);
    /// ```
    pub fn to_bezier_patches(&self) -> Vec<Vec<BezierSurface>> {
        let (u, v) = (&self.u_knots, &self.v_knots);
        let mut rows: Vec<Vec<BezierSurface>> = Vec::with_capacity(u.knots().len() - 1);
        editing::surface_patches(
            u,
            v,
            &self.poles,
            self.weights.as_deref(),
            |i, (u0, u1), (v0, v1), poles, weights| {
                if i == rows.len() {
                    rows.push(Vec::with_capacity(v.knots().len() - 1));
                }
                let patch = BSplineSurface {
                    u_knots: curve::KnotVector::bezier(u.degree(), u0, u1),
                    v_knots: curve::KnotVector::bezier(v.degree(), v0, v1),
                    poles,
                    weights,
                    patches: Lazy::default(),
                };
                rows[i].push(BezierSurface::from_surface(patch));
            },
        );
        rows
    }

    /// The single Bézier patch of the given degrees over `u_range` x
    /// `v_range`, its `(u_degree + 1) (v_degree + 1)` poles row-major and
    /// their weights, if any: the caller has checked the counts.
    pub(crate) fn bezier_patch(
        (u_degree, u_range): (usize, (f64, f64)),
        (v_degree, v_range): (usize, (f64, f64)),
        poles: Vec<Point>,
        weights: Option<Vec<f64>>,
    ) -> BSplineSurface {
        debug_assert_eq!(poles.len(), (u_degree + 1) * (v_degree + 1));
        BSplineSurface {
            u_knots: curve::KnotVector::bezier(u_degree, u_range.0, u_range.1),
            v_knots: curve::KnotVector::bezier(v_degree, v_range.0, v_range.1),
            poles,
            weights,
            patches: Lazy::default(),
        }
    }

    /// Whether the surface is one clamped, non-periodic span in both
    /// directions: its poles are the Bézier control points of the whole
    /// surface.
    pub(crate) fn is_bezier(&self) -> bool {
        [&self.u_knots, &self.v_knots].iter().all(|knots| {
            let end = knots.degree() as u32 + 1;
            !knots.is_periodic() && knots.multiplicities() == [end, end]
        })
    }

    /// Returns the sub-patch of the surface over `u_range` and `v_range`,
    /// clamped and non-periodic in both directions, whose bounds are
    /// exactly the two ranges and which evaluates as this surface there.
    ///
    /// Each direction is cut as [`crate::BSplineCurve::extract`] cuts a
    /// curve, every row or column of the net in the same pass, the `v`
    /// direction only over the rows the `u` cut kept; a range already
    /// equal to the bounds of a clamped direction leaves it alone. A
    /// periodic direction accepts a range across its seam.
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::extract`], for either range against
    /// its direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, BSplineSplitError, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 2, poles, vec![0.0, 1.0], vec![3, 3], vec![0.0, 1.0], vec![3, 3], false, false,
    /// )
    /// .unwrap();
    /// let patch = surface.extract((0.25, 0.5), (0.0, 0.75)).unwrap();
    /// assert_eq!(patch.u_knots(), &[0.25, 0.5]);
    /// assert_eq!(patch.v_knots(), &[0.0, 0.75]);
    /// let (p, q) = (surface.eval_point(0.3, 0.6), patch.eval_point(0.3, 0.6));
    /// assert!(p.distance(q) < 1e-12);
    /// assert_eq!(
    ///     surface.extract((0.5, 0.25), (0.0, 1.0)),
    ///     Err(BSplineSplitError::EmptyRange)
    /// );
    /// ```
    pub fn extract(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<BSplineSurface, BSplineSplitError> {
        self.extract_ranges(Some(u_range), Some(v_range))
    }

    /// Returns the piece of the surface over `[first, last]` in `u`,
    /// clamped and non-periodic in `u`, the `v` direction unchanged; see
    /// [`BSplineSurface::extract`].
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::extract`], against the `u` bounds.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 1, poles, vec![0.0, 1.0], vec![3, 3], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let piece = surface.extract_u(0.5, 1.0).unwrap();
    /// assert_eq!(piece.u_knots(), &[0.5, 1.0]);
    /// assert_eq!(piece.v_knots(), surface.v_knots());
    /// ```
    pub fn extract_u(&self, first: f64, last: f64) -> Result<BSplineSurface, BSplineSplitError> {
        self.extract_ranges(Some((first, last)), None)
    }

    /// Returns the piece of the surface over `[first, last]` in `v`,
    /// clamped and non-periodic in `v`, the `u` direction unchanged; see
    /// [`BSplineSurface::extract`].
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::extract`], against the `v` bounds.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..2)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     1, 2, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![3, 3], false, false,
    /// )
    /// .unwrap();
    /// let piece = surface.extract_v(0.0, 0.4).unwrap();
    /// assert_eq!(piece.v_knots(), &[0.0, 0.4]);
    /// let (p, q) = (surface.eval_point(0.5, 0.2), piece.eval_point(0.5, 0.2));
    /// assert!(p.distance(q) < 1e-12);
    /// ```
    pub fn extract_v(&self, first: f64, last: f64) -> Result<BSplineSurface, BSplineSplitError> {
        self.extract_ranges(None, Some((first, last)))
    }

    /// Splits the surface at `u` into the pieces over `[first, u]` and
    /// `[u, last]` of its `u` bounds, clamped and non-periodic in `u`, the
    /// `v` direction unchanged; they evaluate as this surface over their
    /// ranges and meet along the `u` isoparametric curve at `u`. See
    /// [`crate::BSplineCurve::split_at`].
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::split_at`], against the `u` bounds.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 1, poles, vec![0.0, 1.0], vec![3, 3], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let (low, high) = surface.split_u_at(0.5).unwrap();
    /// // The last row of `low` is the first of `high`.
    /// assert_eq!(&low.poles()[low.poles().len() - 2..], &high.poles()[..2]);
    /// let (p, q) = (surface.eval_point(0.7, 0.3), high.eval_point(0.7, 0.3));
    /// assert!(p.distance(q) < 1e-12);
    /// ```
    pub fn split_u_at(
        &self,
        u: f64,
    ) -> Result<(BSplineSurface, BSplineSurface), BSplineSplitError> {
        editing::check_split(&self.u_knots, u)?;
        let (first, last) = self.u_knots.bounds();
        let [low, high] = self.cut(Direction::U, &[first, u, last]);
        Ok((low, high))
    }

    /// Splits the surface at `v` into the pieces over `[first, v]` and
    /// `[v, last]` of its `v` bounds, clamped and non-periodic in `v`, the
    /// `u` direction unchanged; see [`BSplineSurface::split_u_at`].
    ///
    /// # Errors
    ///
    /// As [`crate::BSplineCurve::split_at`], against the `v` bounds.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..2)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     1, 2, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![3, 3], false, false,
    /// )
    /// .unwrap();
    /// let (low, high) = surface.split_v_at(0.25).unwrap();
    /// assert_eq!(low.v_knots(), &[0.0, 0.25]);
    /// assert_eq!(high.v_knots(), &[0.25, 1.0]);
    /// ```
    pub fn split_v_at(
        &self,
        v: f64,
    ) -> Result<(BSplineSurface, BSplineSurface), BSplineSplitError> {
        editing::check_split(&self.v_knots, v)?;
        let (first, last) = self.v_knots.bounds();
        let [low, high] = self.cut(Direction::V, &[first, v, last]);
        Ok((low, high))
    }

    /// Returns the surface with its `u` direction traversed the other way:
    /// its point at `(u, v)` is this surface's point at
    /// `(first + last - u, v)` over the same bounds, so `Su` changes sign.
    /// The rows of the net are listed in reverse order and the `u` knots
    /// mirrored, as for [`crate::BSplineCurve::reversed`]; a periodic
    /// direction stays periodic.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 1, poles, vec![0.0, 1.0], vec![3, 3], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let reversed = surface.reversed_u();
    /// assert_eq!(reversed.eval_point(0.0, 0.5), surface.eval_point(1.0, 0.5));
    /// assert_eq!(reversed.reversed_u(), surface);
    /// ```
    pub fn reversed_u(&self) -> BSplineSurface {
        let row = self.v_knots.pole_count();
        let order: Vec<usize> = (0..self.u_knots.pole_count())
            .flat_map(|i| {
                let old = editing::reversed_slot(&self.u_knots, i);
                old * row..(old + 1) * row
            })
            .collect();
        BSplineSurface {
            u_knots: editing::reversed_knots(&self.u_knots),
            v_knots: self.v_knots.clone(),
            poles: gathered(&self.poles, &order),
            weights: self.weights.as_ref().map(|w| gathered(w, &order)),
            patches: Lazy::default(),
        }
    }

    /// Returns the surface with its `v` direction traversed the other way:
    /// its point at `(u, v)` is this surface's point at
    /// `(u, first + last - v)`, so `Sv` changes sign; see
    /// [`BSplineSurface::reversed_u`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..2)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     1, 2, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![3, 3], false, false,
    /// )
    /// .unwrap();
    /// let reversed = surface.reversed_v();
    /// assert_eq!(reversed.eval_point(0.5, 0.25), surface.eval_point(0.5, 0.75));
    /// assert_eq!(reversed.reversed_v(), surface);
    /// ```
    pub fn reversed_v(&self) -> BSplineSurface {
        let row = self.v_knots.pole_count();
        let columns: Vec<usize> = (0..row)
            .map(|j| editing::reversed_slot(&self.v_knots, j))
            .collect();
        let order: Vec<usize> = (0..self.u_knots.pole_count())
            .flat_map(|i| columns.iter().map(move |&j| i * row + j))
            .collect();
        BSplineSurface {
            u_knots: self.u_knots.clone(),
            v_knots: editing::reversed_knots(&self.v_knots),
            poles: gathered(&self.poles, &order),
            weights: self.weights.as_ref().map(|w| gathered(w, &order)),
            patches: Lazy::default(),
        }
    }

    /// Maps `u` onto the [`reversed_u`](Self::reversed_u) surface: the
    /// parameter of the same point there, `first + last - u`. The map is
    /// its own inverse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 1, poles, vec![1.0, 4.0], vec![3, 3], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let u = surface.u_reversed_parameter(1.5);
    /// assert_eq!(u, 3.5);
    /// assert_eq!(surface.reversed_u().eval_point(u, 0.5), surface.eval_point(1.5, 0.5));
    /// ```
    pub fn u_reversed_parameter(&self, u: f64) -> f64 {
        let (first, last) = self.u_param_range();
        first + last - u
    }

    /// Maps `v` onto the [`reversed_v`](Self::reversed_v) surface:
    /// `first + last - v`, as [`u_reversed_parameter`](Self::u_reversed_parameter)
    /// for `u`.
    pub fn v_reversed_parameter(&self, v: f64) -> f64 {
        let (first, last) = self.v_param_range();
        first + last - v
    }

    /// The `u` isoparametric curve at `u`: `u` fixed, running along `v`
    /// with the `v` degree, knots and periodicity, so
    /// `u_iso(u).eval_point(t)` is `eval_point(u, t)` for every `t`.
    ///
    /// Each pole is a column of the net blended by the `u` basis functions
    /// at `u`, in homogeneous coordinates if rational. The curve is
    /// rational unless its weights all come out equal, as for a surface
    /// rational in `u` alone. On a collapsed edge (a pole of a surface of
    /// revolution) every pole coincides and the curve is that point.
    ///
    /// # Errors
    ///
    /// Returns [`SurfaceIsoError::NonFinite`] for a NaN or infinite `u`
    /// and [`SurfaceIsoError::OutOfRange`] for a `u` outside the bounds of
    /// a non-periodic `u` direction; a periodic one wraps it into the
    /// period.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 1, poles, vec![0.0, 1.0], vec![3, 3], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let iso = surface.u_iso(0.25).unwrap();
    /// assert_eq!(iso.degree(), 1);
    /// assert!(iso.eval_point(0.5).distance(surface.eval_point(0.25, 0.5)) < 1e-12);
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<BSplineCurve, SurfaceIsoError> {
        let u = check_iso_param(u, Some(&self.u_knots))?;
        Ok(self.iso(Direction::U, u))
    }

    /// The `v` isoparametric curve at `v`: `v` fixed, running along `u`
    /// with the `u` degree, knots and periodicity, so
    /// `v_iso(v).eval_point(t)` is `eval_point(t, v)` for every `t`. See
    /// [`BSplineSurface::u_iso`].
    ///
    /// # Errors
    ///
    /// Returns [`SurfaceIsoError::NonFinite`] for a NaN or infinite `v`
    /// and [`SurfaceIsoError::OutOfRange`] for a `v` outside the bounds of
    /// a non-periodic `v` direction; a periodic one wraps it into the
    /// period.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point};
    ///
    /// let poles: Vec<Vec<Point>> = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BSplineSurface::new(
    ///     2, 1, poles, vec![0.0, 1.0], vec![3, 3], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let iso = surface.v_iso(0.5).unwrap();
    /// assert_eq!(iso.degree(), 2);
    /// assert!(iso.eval_point(0.25).distance(surface.eval_point(0.25, 0.5)) < 1e-12);
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<BSplineCurve, SurfaceIsoError> {
        let v = check_iso_param(v, Some(&self.v_knots))?;
        Ok(self.iso(Direction::V, v))
    }

    /// The isoparametric curve with `direction` fixed at `t`, a checked
    /// parameter: every lane of the net across `direction` blended by its
    /// basis functions at `t`.
    fn iso(&self, direction: Direction, t: f64) -> BSplineCurve {
        let (fixed, along) = match direction {
            Direction::U => (&self.u_knots, &self.v_knots),
            Direction::V => (&self.v_knots, &self.u_knots),
        };
        let n_v = self.v_knots.pole_count();
        let index = |i: usize, j: usize| match direction {
            Direction::U => i * n_v + j,
            Direction::V => j * n_v + i,
        };
        let mut basis = [0.0f64; curve::LARGE];
        let first = fixed.basis(t, &mut basis);
        let basis = &basis[..=fixed.degree()];
        let net = self.pole_view();

        let mut poles = Vec::with_capacity(along.pole_count());
        let mut weights = Vec::with_capacity(along.pole_count());
        for j in 0..along.pole_count() {
            let mut sum = [0.0f64; 4];
            for (k, &b) in basis.iter().enumerate() {
                let h = net.homogeneous(index(fixed.pole_index(first, k), j));
                for c in 0..4 {
                    sum[c] += b * h[c];
                }
            }
            if self.weights.is_some() {
                let inverse = 1.0 / sum[3];
                poles.push(Point::new(
                    sum[0] * inverse,
                    sum[1] * inverse,
                    sum[2] * inverse,
                ));
                weights.push(sum[3]);
            } else {
                poles.push(Point::new(sum[0], sum[1], sum[2]));
            }
        }
        let rational = weights.windows(2).any(|w| w[0] != w[1]);
        BSplineCurve::from_validated_parts(
            along.degree(),
            poles,
            rational.then_some(weights),
            along.knots().to_vec(),
            along.multiplicities().to_vec(),
            along.flat().to_vec(),
            along.is_periodic(),
        )
    }

    /// The piece over `u_range` and `v_range`, either left whole if
    /// `None`: see [`BSplineSurface::extract`].
    fn extract_ranges(
        &self,
        u_range: Option<(f64, f64)>,
        v_range: Option<(f64, f64)>,
    ) -> Result<BSplineSurface, BSplineSplitError> {
        if let Some((first, last)) = u_range {
            editing::check_piece(&self.u_knots, first, last)?;
        }
        if let Some((first, last)) = v_range {
            editing::check_piece(&self.v_knots, first, last)?;
        }
        // A range that is the bounds of a clamped direction leaves it whole.
        let cut = |knots, range: Option<(f64, f64)>| {
            range.filter(|&(first, last)| !editing::is_clamped_to(knots, first, last))
        };
        let u_piece = cut(&self.u_knots, u_range).map(|(first, last)| {
            let [piece] = self.cut(Direction::U, &[first, last]);
            piece
        });
        let surface = u_piece.as_ref().unwrap_or(self);
        match cut(&self.v_knots, v_range) {
            Some((first, last)) => {
                let [piece] = surface.cut(Direction::V, &[first, last]);
                Ok(piece)
            }
            None => Ok(u_piece.unwrap_or_else(|| self.clone())),
        }
    }

    /// The surface's pieces between the increasing `cuts` along
    /// `direction`, checked against its bounds.
    fn cut<const N: usize>(&self, direction: Direction, cuts: &[f64]) -> [BSplineSurface; N] {
        let (knots, lanes) = match direction {
            Direction::U => (&self.u_knots, Lanes::columns(self.v_knots.pole_count())),
            Direction::V => (&self.v_knots, Lanes::rows(self.u_knots.pole_count())),
        };
        let pieces = editing::cut(knots, &self.poles, self.weights.as_deref(), lanes, cuts);
        let pieces: Vec<BSplineSurface> = pieces
            .into_iter()
            .map(|piece| {
                let (u_knots, v_knots) = match direction {
                    Direction::U => (piece.knots, self.v_knots.clone()),
                    Direction::V => (self.u_knots.clone(), piece.knots),
                };
                BSplineSurface {
                    u_knots,
                    v_knots,
                    poles: piece.points,
                    weights: piece.weights,
                    patches: Lazy::default(),
                }
            })
            .collect();
        pieces
            .try_into()
            .unwrap_or_else(|_| unreachable!("a piece between each two cuts"))
    }

    /// The surface with each of `knots` removed once in `direction`, and
    /// how many copies were: see [`BSplineSurface::remove_u_knots`].
    fn removed(
        &self,
        direction: Direction,
        knots: &[f64],
        max_deviation: f64,
    ) -> Result<(BSplineSurface, usize), KnotRemovalError> {
        let (direction_knots, lanes) = match direction {
            Direction::U => (&self.u_knots, Lanes::columns(self.v_knots.pole_count())),
            Direction::V => (&self.v_knots, Lanes::rows(self.u_knots.pole_count())),
        };
        editing::check_removal(direction_knots, knots, max_deviation)?;
        let removal = editing::remove_knots(
            direction_knots,
            &self.poles,
            self.weights.as_deref(),
            lanes,
            knots,
            max_deviation,
        );
        let Some(removed) = removal else {
            return Ok((self.clone(), 0));
        };
        let (u_knots, v_knots) = match direction {
            Direction::U => (removed.knots, self.v_knots.clone()),
            Direction::V => (self.u_knots.clone(), removed.knots),
        };
        let surface = BSplineSurface {
            u_knots,
            v_knots,
            poles: removed.points,
            weights: removed.weights,
            patches: Lazy::default(),
        };
        Ok((surface, removed.count))
    }

    /// The surface with one direction edited by `edit`, every column (in
    /// `u`) or row (in `v`) of the net in the same pass.
    fn edited<E: Edit>(&self, direction: Direction, edit: &E) -> Result<BSplineSurface, E::Error> {
        let (knots, lanes) = match direction {
            Direction::U => (&self.u_knots, Lanes::columns(self.v_knots.pole_count())),
            Direction::V => (&self.v_knots, Lanes::rows(self.u_knots.pole_count())),
        };
        let (knots, poles, weights) =
            editing::edit(edit, knots, &self.poles, self.weights.as_deref(), lanes)?;
        let (u_knots, v_knots) = match direction {
            Direction::U => (knots, self.v_knots.clone()),
            Direction::V => (self.u_knots.clone(), knots),
        };
        Ok(BSplineSurface {
            u_knots,
            v_knots,
            poles,
            weights,
            patches: Lazy::default(),
        })
    }

    /// Returns the image of this surface under `transform`: every control
    /// pole maps point-wise while degrees, knots, multiplicities,
    /// periodicity, and weights are unchanged.
    ///
    /// B-spline evaluation is affinely invariant, so mapping the poles is
    /// exact for every rigid, mirrored, or uniformly scaled transform.
    /// Infallible: pole images always form a valid net over the unchanged
    /// knot vectors, even for a collapsing (zero-scale) transform.
    ///
    /// # Examples
    ///
    /// ```
    /// # use geomcore::prelude::*;
    /// # use geomcore::surfaces::BSplineSurface;
    /// # use geomcore::{Point, Transform, Vector};
    /// # let poles = vec![
    /// #     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    /// #     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// # ];
    /// # let surface = BSplineSurface::new(1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false).unwrap();
    /// let moved = surface.transformed(Transform::translation(Vector::new(0.0, 0.0, 4.0)));
    /// assert_eq!(moved.eval_point(0.5, 0.5), Point::new(0.5, 0.5, 4.25));
    /// ```
    pub fn transformed(&self, transform: Transform) -> BSplineSurface {
        BSplineSurface {
            u_knots: self.u_knots.clone(),
            v_knots: self.v_knots.clone(),
            poles: self
                .poles
                .iter()
                .map(|&p| transform.apply_point(p))
                .collect(),
            weights: self.weights.clone(),
            patches: Lazy::default(),
        }
    }

    /// Evaluates the point and the first derivatives at `(u, v)`, each
    /// direction from its own side: [`Side::Below`] takes the span ending
    /// at the parameter, [`Side::Above`] the span starting there. Both
    /// sides agree strictly inside a span. Plain evaluation takes the
    /// above side, except at the upper bound, where nothing lies above a
    /// non-periodic direction and both sides take the last span. A
    /// periodic direction's seam is a knot like any other: at either
    /// bound, [`Side::Below`] takes the last span of the period and
    /// [`Side::Above`] the first. The two sides of a direction can differ
    /// only in derivatives of higher order than
    /// [`continuity_u`](Self::continuity_u) or
    /// [`continuity_v`](Self::continuity_v) reports there.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineSurface, Point, Side, Vector};
    /// // Two bilinear patches meeting at a crease along u = 1, flat then rising.
    /// let flat = |x: f64| vec![Point::new(x, 0.0, 0.0), Point::new(x, 1.0, 0.0)];
    /// let tilted = |x: f64| vec![Point::new(x, 0.0, 0.0), Point::new(x, 1.0, 1.0)];
    /// let surface = BSplineSurface::new(
    ///     1, 1,
    ///     vec![flat(0.0), flat(1.0), tilted(2.0)],
    ///     vec![0.0, 1.0, 2.0], vec![2, 1, 2], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let below = surface.eval_d1_sided(1.0, 0.5, Side::Below, Side::Above);
    /// let above = surface.eval_d1_sided(1.0, 0.5, Side::Above, Side::Above);
    /// assert_eq!(below.1[0], Vector::new(1.0, 0.0, 0.0));
    /// assert_eq!(above.1[0], Vector::new(1.0, 0.0, 0.5));
    /// ```
    pub fn eval_d1_sided(
        &self,
        u: f64,
        v: f64,
        u_side: Side,
        v_side: Side,
    ) -> (Point, [Vector; 2]) {
        let (s, su, sv) = self.eval_surface_d1_sided(u, v, u_side, v_side);
        (s, [su, sv])
    }

    /// Evaluates the derivative of order `(du, dv)` at `(u, v)`, each
    /// direction from its own side, as
    /// [`eval_d1_sided`](Self::eval_d1_sided).
    ///
    /// Only orders with `du + dv` of 1 or 2 are supported.
    ///
    /// # Panics
    ///
    /// Panics if `du + dv == 0` or `du + dv > 2`, as
    /// [`eval_derivative`](Self::eval_derivative).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineSurface, Point, Side, Vector};
    /// // Two bilinear patches meeting at a crease along u = 1, flat then rising.
    /// let flat = |x: f64| vec![Point::new(x, 0.0, 0.0), Point::new(x, 1.0, 0.0)];
    /// let tilted = |x: f64| vec![Point::new(x, 0.0, 0.0), Point::new(x, 1.0, 1.0)];
    /// let surface = BSplineSurface::new(
    ///     1, 1,
    ///     vec![flat(0.0), flat(1.0), tilted(2.0)],
    ///     vec![0.0, 1.0, 2.0], vec![2, 1, 2], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// // The cross derivative is zero on the flat patch, (0, 0, 1) on the tilted one.
    /// let below = surface.eval_derivative_sided(1.0, 0.5, 1, 1, Side::Below, Side::Above);
    /// let above = surface.eval_derivative_sided(1.0, 0.5, 1, 1, Side::Above, Side::Above);
    /// assert_eq!(below, Vector::new(0.0, 0.0, 0.0));
    /// assert_eq!(above, Vector::new(0.0, 0.0, 1.0));
    /// ```
    pub fn eval_derivative_sided(
        &self,
        u: f64,
        v: f64,
        du: u32,
        dv: u32,
        u_side: Side,
        v_side: Side,
    ) -> Vector {
        match (du, dv) {
            (0, 0) => panic!(
                "eval_derivative: du + dv must be >= 1 (use eval_point for the (0, 0) order)"
            ),
            (1, 0) | (0, 1) => {
                let (_s, su, sv) = self.eval_surface_d1_sided(u, v, u_side, v_side);
                if du == 1 { su } else { sv }
            }
            (2, 0) | (1, 1) | (0, 2) => {
                let [_, _, _, suu, suv, svv] = self.eval_surface_d2_sided(u, v, u_side, v_side);
                match du {
                    2 => suu,
                    1 => suv,
                    _ => svv,
                }
            }
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    /// Evaluates the point, the first and the second derivatives at
    /// `(u, v)`, each direction from its own side, as
    /// [`eval_d1_sided`](Self::eval_d1_sided).
    pub fn eval_d2_sided(
        &self,
        u: f64,
        v: f64,
        u_side: Side,
        v_side: Side,
    ) -> (Point, [Vector; 2], [Vector; 3]) {
        let [s, su, sv, suu, suv, svv] = self.eval_surface_d2_sided(u, v, u_side, v_side);
        (Point::ORIGIN + s, [su, sv], [suu, suv, svv])
    }

    /// All local minima of the distance from `point` to the surface,
    /// ordered by ascending distance, the first being
    /// [`BSplineSurface::project_point`].
    ///
    /// Every knot-span patch is searched as in
    /// [`BSplineSurface::project_point`], but keeping each minimum found
    /// rather than only the closest: a subpatch is dropped only when it
    /// provably holds no local minimum. Each minimum is polished by the
    /// same Newton minimisation within the C1 region around it (C0 crease
    /// knots bound it), so minima on an edge, corner or crease are reached
    /// exactly. Minima within `tol.confusion` of each other merge.
    /// Periodic parameters are reported within their period.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, Point, Tolerance};
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 2.0, 0.0)],
    ///     vec![Point::new(2.0, 0.0, 0.0), Point::new(2.0, 2.0, 0.0)],
    /// ];
    /// let surface = BSplineSurface::new(
    ///     1, 1, poles, vec![0.0, 1.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let extrema = surface.extrema(Point::new(1.0, 1.0, 1.0), Tolerance::DEFAULT);
    /// assert!(!extrema.is_empty());
    /// assert_eq!(extrema[0].distance, 1.0);
    /// ```
    pub fn extrema(&self, point: Point, tol: Tolerance) -> Vec<SurfaceProjection> {
        use crate::projection::snap_distance;
        let global = self.project_point(point, tol);
        let patches = self.patch_tree();
        let mut search = PatchSearch::new(self, point, patches.resolution(point), true);
        for k in 0..patches.len() {
            let (rect, creases) = patches.patch(k);
            if let Some(root) = search.root(patches, k) {
                search.search(root, rect, creases);
            }
        }
        let u_regions = self.u_knots.smooth_regions();
        let v_regions = self.v_knots.smooth_regions();
        let (u_period, v_period) = (self.u_knots.period(), self.v_knots.period());
        // A minimum on a subpatch edge inside a smooth region is where the
        // search stopped, not necessarily a minimum of the surface: polish
        // across subpatches.
        let mut candidates = vec![(global.u, global.v)];
        for (u, v) in search.minima {
            for (ubounds, u) in curve::regions_at(&u_regions, u, u_period) {
                for (vbounds, v) in curve::regions_at(&v_regions, v, v_period) {
                    candidates.push(self.refine_projection(point, u, v, [ubounds, vbounds]));
                }
            }
        }
        let mut out: Vec<SurfaceProjection> = Vec::new();
        let mut kept: Vec<Point> = Vec::new();
        for (u, v) in candidates {
            let p = self.eval_point(u, v);
            if kept.iter().any(|q| q.distance(p) <= tol.confusion) {
                continue;
            }
            kept.push(p);
            out.push(SurfaceProjection {
                u,
                v,
                point: p,
                distance: snap_distance(p.distance(point), tol.confusion),
            });
        }
        // The projection leads, even past a rounding-level tie.
        out[1..].sort_by(|x, y| x.distance.total_cmp(&y.distance));
        out
    }

    /// The shape of a knot-span patch's Bernstein net: `(rows, cols, h)`,
    /// `u_degree + 1` rows of `v_degree + 1` points of `h` values
    /// (homogeneous when `h == 4`).
    pub(crate) fn patch_net_shape(&self) -> (usize, usize, usize) {
        let h = if self.weights.is_some() { 4 } else { 3 };
        (self.u_degree() + 1, self.v_degree() + 1, h)
    }

    /// The bounding-box hierarchy of the knot-span patches (leaf `k` is
    /// patch `k`, see [`Self::patch_net`]), built on first use.
    pub(crate) fn patch_boxes(&self) -> &BoxTree {
        &self.patch_tree().tree
    }

    /// Writes the Bernstein net of knot-span patch `k` to `out` (laid out
    /// as [`math::patch_bezier`]; `scratch` holds as many values) and
    /// returns its parameter rectangle `[u0, u1, v0, v1]`.
    pub(crate) fn patch_net(&self, k: usize, scratch: &mut [f64], out: &mut [f64]) -> [f64; 4] {
        let patches = self.patch_tree();
        self.write_patch_net(patches, k, scratch, out);
        patches.patch(k).0
    }

    fn write_patch_net(&self, patches: &PatchTree, k: usize, scratch: &mut [f64], out: &mut [f64]) {
        let n_vspans = patches.v_breaks.len() - 1;
        let (a, b) = (k / n_vspans, k % n_vspans);
        let (p, q) = (self.u_degree(), self.v_degree());
        let (mu, mv) = ((p + 1) * (p + 1), (q + 1) * (q + 1));
        math::patch_bezier(
            (
                &patches.u_spans.1[a * mu..(a + 1) * mu],
                patches.u_spans.0[a],
                p,
                self.u_knots.pole_count(),
            ),
            (
                &patches.v_spans.1[b * mv..(b + 1) * mv],
                patches.v_spans.0[b],
                q,
                self.v_knots.pole_count(),
            ),
            self.pole_view(),
            scratch,
            out,
        );
    }

    /// The patch bounding-box hierarchy, built on first use.
    fn patch_tree(&self) -> &PatchTree {
        self.patches.get_or_init(|| {
            let u_breaks = self.u_knots.breaks();
            let v_breaks = self.v_knots.breaks();
            let boxes = math::patch_boxes(
                &self.u_knots,
                &self.v_knots,
                self.pole_view(),
                &u_breaks,
                &v_breaks,
            );
            PatchTree {
                tree: BoxTree::new(&boxes, v_breaks.len() - 1),
                u_spans: self.u_knots.span_bezier_matrices(&u_breaks),
                v_spans: self.v_knots.span_bezier_matrices(&v_breaks),
                u_creases: self.u_knots.crease_breaks(&u_breaks),
                v_creases: self.v_knots.crease_breaks(&v_breaks),
                u_breaks,
                v_breaks,
            }
        })
    }

    /// Shared first-derivative evaluation, returning `(S, Su, Sv)`.
    fn eval_surface_d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
        let [s, su, sv, ..] =
            math::surface_derivatives::<1>(&self.u_knots, &self.v_knots, self.pole_view(), u, v);
        (Point::ORIGIN + s, su, sv)
    }

    /// Shared sided first-derivative evaluation: `u_side` selects the `u`
    /// span, `v_side` the `v` span.
    fn eval_surface_d1_sided(
        &self,
        u: f64,
        v: f64,
        u_side: Side,
        v_side: Side,
    ) -> (Point, Vector, Vector) {
        let [s, su, sv, ..] = math::surface_derivatives_sided::<1>(
            &self.u_knots,
            &self.v_knots,
            self.pole_view(),
            u,
            v,
            u_side,
            v_side,
        );
        (Point::ORIGIN + s, su, sv)
    }

    /// Shared second-derivative evaluation, returning
    /// `[S, Su, Sv, Suu, Suv, Svv]` (the point as a vector).
    fn eval_surface_d2(&self, u: f64, v: f64) -> [Vector; 6] {
        math::surface_derivatives::<2>(&self.u_knots, &self.v_knots, self.pole_view(), u, v)
    }

    /// Shared sided second-derivative evaluation, as
    /// [`eval_surface_d1_sided`](Self::eval_surface_d1_sided).
    fn eval_surface_d2_sided(&self, u: f64, v: f64, u_side: Side, v_side: Side) -> [Vector; 6] {
        math::surface_derivatives_sided::<2>(
            &self.u_knots,
            &self.v_knots,
            self.pole_view(),
            u,
            v,
            u_side,
            v_side,
        )
    }

    /// The poles as read by evaluation.
    fn pole_view(&self) -> curve::Poles<'_> {
        curve::Poles {
            points: &self.poles,
            weights: self.weights.as_deref(),
        }
    }

    /// Values per point of a Bernstein net: homogeneous if rational.
    fn net_stride(&self) -> usize {
        if self.is_rational() { 4 } else { 3 }
    }

    /// Box-bounded Newton minimisation of `|S(u,v) - point|^2` from the
    /// seed `(u, v)`, with exact second derivatives, within `bounds` per
    /// direction (`None`: unbounded, for a smooth periodic direction);
    /// periodic parameters are wrapped into their period on return.
    fn refine_projection(
        &self,
        point: Point,
        u: f64,
        v: f64,
        bounds: [Option<(f64, f64)>; 2],
    ) -> (f64, f64) {
        let (u0, u1) = self.u_param_range();
        let (v0, v1) = self.v_param_range();
        let target = Vector::new(point.x, point.y, point.z);
        // Minimise in each subpatch's own `[0, 1]` coordinates. In global
        // parameters, a span 1e-9 wide beside a full-width one makes
        // Newton's Hessian regularisation damp the wide direction, and its
        // step tolerance stop before moving a visible distance in space.
        let span = [
            bounds[0].map(|(lo, hi)| hi - lo),
            bounds[1].map(|(lo, hi)| hi - lo),
        ];
        let local_bounds = [bounds[0].map(|_| (0.0, 1.0)), bounds[1].map(|_| (0.0, 1.0))];
        let to_local = |x: f64, i: usize| match bounds[i] {
            Some((lo, _)) => (x - lo) / span[i].unwrap_or(1.0),
            None => x,
        };
        // `lo + 1.0 * (hi - lo)` may round past `hi`.
        let to_global = |s: f64, i: usize| match bounds[i] {
            Some((lo, hi)) => (lo + s * span[i].unwrap_or(1.0)).min(hi),
            None => s,
        };
        let objective = |s: f64, t: f64| {
            let (du, dv) = (span[0].unwrap_or(1.0), span[1].unwrap_or(1.0));
            let (u, v) = (to_global(s, 0), to_global(t, 1));
            let [p, su, sv, suu, suv, svv] =
                self.eval_surface_d2(curve::inward(u, bounds[0]), curve::inward(v, bounds[1]));
            let (su, sv) = (su * du, sv * dv);
            let (suu, suv, svv) = (suu * du * du, suv * du * dv, svv * dv * dv);
            let r = p - target;
            let huv = su.dot(sv) + suv.dot(r);
            (
                0.5 * r.dot(r),
                [su.dot(r), sv.dot(r)],
                [
                    [su.dot(su) + suu.dot(r), huv],
                    [huv, sv.dot(sv) + svv.dot(r)],
                ],
            )
        };
        let [s, t] = newton_min_2d(
            objective,
            [to_local(u, 0), to_local(v, 1)],
            local_bounds,
            1e-14,
            100,
        );
        let (u, v) = (to_global(s, 0), to_global(t, 1));
        let wrap = |x: f64, periodic: bool, lo: f64, hi: f64| {
            if periodic {
                lo + (x - lo).rem_euclid(hi - lo)
            } else {
                x
            }
        };
        (
            wrap(u, self.is_u_periodic(), u0, u1),
            wrap(v, self.is_v_periodic(), v0, v1),
        )
    }

    /// `u` parameter bounds `(first, last)`: for a clamped direction the active
    /// flat span, for a periodic direction the full period `(knots[0],
    /// knots[last])`.
    fn u_param_range(&self) -> (f64, f64) {
        self.u_knots.bounds()
    }

    /// `v` parameter bounds `(first, last)`, analogous to [`Self::u_param_range`].
    fn v_param_range(&self) -> (f64, f64) {
        self.v_knots.bounds()
    }
}

impl ParametricSurface for BSplineSurface {
    fn eval_point(&self, u: f64, v: f64) -> Point {
        math::surface_point(&self.u_knots, &self.v_knots, self.pole_view(), u, v)
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        use crate::projection::snap_distance;
        let patches = self.patch_tree();
        let mut search = PatchSearch::new(self, point, patches.resolution(point), false);
        // Every patch that may hold the closest point, by their pole boxes
        // then their Bernstein nets (whose corners tighten the bound),
        // searched closest control point first. Only the first keeps its
        // net while waiting, the others' written again when searched.
        let mut roots = Vec::with_capacity(RESERVED_SLOTS);
        let mut first: Option<Subpatch> = None;
        patches.tree.nearest_first(point, |k| {
            if let Some(root) = search.root(patches, k) {
                roots.push((root.key, root.lower, k));
                let freed = if first.is_none_or(|held| root.key < held.key) {
                    first.replace(root)
                } else {
                    Some(root)
                };
                if let Some(freed) = freed {
                    search.free.push(freed.slot);
                }
            }
            search.best.0
        });
        // A stable sort, so that the first held stays first among ties.
        roots.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (index, (key, lower, k)) in roots.into_iter().enumerate() {
            if lower < search.best.0 - search.resolution {
                let slot = match first {
                    Some(held) if index == 0 => held.slot,
                    _ => {
                        let slot = search.slot();
                        search.write_root(patches, k, slot);
                        slot
                    }
                };
                let (rect, creases) = patches.patch(k);
                let root = Subpatch {
                    lower,
                    key,
                    slot,
                    rect,
                    depth: 0,
                    edge: None,
                };
                search.search(root, rect, creases);
            }
        }
        let (distance, u, v) = search.best;
        // The search tracks distances, not points: one evaluation.
        SurfaceProjection {
            u,
            v,
            point: self.eval_point(u, v),
            distance: snap_distance(distance, tol.confusion),
        }
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        match (du, dv) {
            (0, 0) => panic!(
                "eval_derivative: du + dv must be >= 1 (use eval_point for the (0, 0) order)"
            ),
            (1, 0) | (0, 1) => {
                let (_s, su, sv) = self.eval_surface_d1(u, v);
                if du == 1 { su } else { sv }
            }
            (2, 0) | (1, 1) | (0, 2) => {
                let [_, _, _, suu, suv, svv] = self.eval_surface_d2(u, v);
                match du {
                    2 => suu,
                    1 => suv,
                    _ => svv,
                }
            }
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    fn eval_d1(&self, u: f64, v: f64) -> (Point, [Vector; 2]) {
        let (s, su, sv) = self.eval_surface_d1(u, v);
        (s, [su, sv])
    }

    fn eval_d2(&self, u: f64, v: f64) -> (Point, [Vector; 2], [Vector; 3]) {
        let [s, su, sv, suu, suv, svv] = self.eval_surface_d2(u, v);
        (Point::ORIGIN + s, [su, sv], [suu, suv, svv])
    }

    fn u_bounds(&self) -> (f64, f64) {
        self.u_param_range()
    }

    fn v_bounds(&self) -> (f64, f64) {
        self.v_param_range()
    }

    fn u_period(&self) -> Option<f64> {
        self.u_knots.period()
    }

    fn v_period(&self) -> Option<f64> {
        self.v_knots.period()
    }

    fn continuity_u(&self, u: f64, tol: Tolerance) -> Continuity {
        curve::continuity_at(
            self.u_degree(),
            self.u_knots(),
            self.u_multiplicities(),
            self.u_param_range(),
            self.is_u_periodic(),
            u,
            tol,
        )
    }

    fn continuity_v(&self, v: f64, tol: Tolerance) -> Continuity {
        curve::continuity_at(
            self.v_degree(),
            self.v_knots(),
            self.v_multiplicities(),
            self.v_param_range(),
            self.is_v_periodic(),
            v,
            tol,
        )
    }

    fn min_continuity_u(&self) -> Continuity {
        curve::min_continuity(
            self.u_degree(),
            self.u_knots(),
            self.u_multiplicities(),
            self.u_param_range(),
            self.is_u_periodic(),
        )
    }

    fn min_continuity_v(&self) -> Continuity {
        curve::min_continuity(
            self.v_degree(),
            self.v_knots(),
            self.v_multiplicities(),
            self.v_param_range(),
            self.is_v_periodic(),
        )
    }

    fn is_u_closed(&self, tol: Tolerance) -> bool {
        self.is_u_periodic() || self.boundaries_meet(Direction::U, tol)
    }

    fn is_v_closed(&self, tol: Tolerance) -> bool {
        self.is_v_periodic() || self.boundaries_meet(Direction::V, tol)
    }

    fn as_bspline(&self) -> Option<&BSplineSurface> {
        Some(self)
    }
}

/// Subdivisions per knot-span patch before the search settles for Newton
/// from the closest remaining subpatches: reached only where the distance
/// is flat to the search resolution. Once one patch of a closest-point
/// search reaches it, the later patches get an eighth of it.
const MAX_SPLITS: usize = 256;

/// Net slots reserved when a search starts: a closest point search holds
/// a few patches and their pending subpatches, about this many.
const RESERVED_SLOTS: usize = 8;

/// Newton refinements once [`MAX_SPLITS`] is reached.
const MAX_FALLBACKS: usize = 8;
/// Subdivisions from a patch down to one subpatch: halving below this,
/// parameter steps reach rounding.
const MAX_DEPTH: u32 = 80;

/// Branch and bound over the Bernstein subpatches of knot-span patches for
/// the minima of the distance to a point (see
/// [`BSplineSurface::project_point`]).
struct PatchSearch<'a> {
    surface: &'a BSplineSurface,
    point: Point,
    /// Net shape: `rows x cols` points of `h` values, `size` in all.
    rows: usize,
    cols: usize,
    h: usize,
    size: usize,
    /// Storage of the pending subpatches' nets, `size` values per slot.
    nets: Vec<f64>,
    /// Offsets of unused slots in `nets`.
    free: Vec<usize>,
    pending: BinaryHeap<Subpatch>,
    /// Euclidean control points of the subpatch at hand.
    points: Vec<[f64; 3]>,
    /// A copy of the net being split.
    scratch: Vec<f64>,
    /// `(distance, u, v)` of the closest point so far.
    best: (f64, f64, f64),
    /// Subpatches closer than this to `best` need no search.
    resolution: f64,
    /// Whether every local minimum is wanted, not only the closest.
    all: bool,
    /// The minima found, when `all`.
    minima: Vec<(f64, f64)>,
    /// Subdivisions allowed per patch.
    budget: usize,
}

/// A subpatch awaiting its visit: its net at offset `slot`, its parameter
/// rectangle `[u0, u1, v0, v1]`, a lower bound `lower` on its distance,
/// and the visiting order `key`, its closest control point's distance. Once
/// its minima are known to lie on its edge `side` of `direction`, `edge`
/// holds `(direction, side)` and the net is that edge curve's.
#[derive(Clone, Copy)]
struct Subpatch {
    lower: f64,
    key: f64,
    slot: usize,
    rect: [f64; 4],
    depth: u32,
    edge: Option<(usize, usize)>,
}

impl PartialEq for Subpatch {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Subpatch {}

impl PartialOrd for Subpatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Subpatch {
    /// Reversed, so a [`BinaryHeap`] pops the lowest key first.
    fn cmp(&self, other: &Self) -> Ordering {
        other.key.total_cmp(&self.key)
    }
}

/// What a subpatch's derivative bounds say about its minima.
enum Verdict {
    /// No local minimum of the surface: the gradient keeps a sign towards
    /// a smooth edge, or the distance curves down everywhere.
    None,
    /// At the corner `[s, t]` (each `0` or `1`), the gradient pinning it
    /// to two crease or boundary edges.
    Corner([usize; 2]),
    /// On edge `side` (`0` or `1`) of direction `direction` (`0` for `u`),
    /// pinned there by the gradient, convex along it.
    Edge(usize, usize),
    /// On edge `side` of `direction`, pinned there by the gradient, but
    /// anywhere along it as far as known.
    Pinned(usize, usize),
    /// Anywhere, the distance being convex.
    Convex,
    /// Not known yet: split along `u` if true, else along `v`.
    Unknown(bool),
}

impl<'a> PatchSearch<'a> {
    /// An empty search for the minima of the distance from `point` (every
    /// local one if `all`), subpatches within `resolution` of the best
    /// distance so far being settled.
    fn new(surface: &'a BSplineSurface, point: Point, resolution: f64, all: bool) -> Self {
        let (rows, cols) = (surface.u_degree() + 1, surface.v_degree() + 1);
        let h = surface.net_stride();
        let size = rows * cols * h;
        PatchSearch {
            surface,
            point,
            rows,
            cols,
            h,
            size,
            nets: Vec::with_capacity(RESERVED_SLOTS * size),
            free: Vec::with_capacity(RESERVED_SLOTS),
            pending: BinaryHeap::with_capacity(RESERVED_SLOTS),
            points: Vec::with_capacity(rows * cols),
            scratch: vec![0.0; size],
            best: (f64::INFINITY, 0.0, 0.0),
            resolution,
            all,
            minima: Vec::new(),
            budget: MAX_SPLITS,
        }
    }

    /// The whole patch `k` of `patches` as a subpatch to search, unless it
    /// is farther than the closest point so far.
    fn root(&mut self, patches: &PatchTree, k: usize) -> Option<Subpatch> {
        let slot = self.slot();
        self.write_root(patches, k, slot);
        self.node(slot, patches.patch(k).0, 0, None)
    }

    /// Writes the net of the whole patch `k` of `patches` at `slot`.
    fn write_root(&mut self, patches: &PatchTree, k: usize, slot: usize) {
        self.surface.write_patch_net(
            patches,
            k,
            &mut self.scratch,
            &mut self.nets[slot..slot + self.size],
        );
    }

    /// Searches the patch over `rect`, with the subpatch `root` over all of
    /// it, whose edges are creases or boundaries as flagged in `creases`.
    fn search(&mut self, root: Subpatch, rect: [f64; 4], creases: [bool; 4]) {
        self.pending.push(root);
        let (mut splits, mut fallbacks) = (0, 0);
        while let Some(node) = self.pending.pop() {
            if !self.all && node.lower >= self.best.0 - self.resolution {
                self.free.push(node.slot);
                continue;
            }
            // Only the patch's own edges can be creases.
            let edge = std::array::from_fn(|i| creases[i] && node.rect[i] == rect[i]);
            let verdict = self.verdict(&node, edge);
            match verdict {
                Verdict::None => {}
                Verdict::Corner([s, t]) => {
                    let p = self.points[s * (self.rows - 1) * self.cols + t * (self.cols - 1)];
                    let distance = Point::new(p[0], p[1], p[2]).distance(self.point);
                    self.found(distance, node.rect[s], node.rect[2 + t]);
                }
                Verdict::Edge(direction, side) => {
                    let start = self.closest_control_point(Some((direction, side)));
                    self.minimise(&node, start);
                }
                Verdict::Pinned(direction, side) => {
                    // Search the edge curve instead, keeping the slot.
                    math::collapse_net(
                        &mut self.nets[node.slot..node.slot + self.size],
                        self.rows,
                        self.cols,
                        self.h,
                        direction == 0,
                        side,
                    );
                    self.push(node.slot, node.rect, node.depth, Some((direction, side)));
                    continue;
                }
                Verdict::Convex => {
                    let start = self.closest_control_point(None);
                    self.minimise(&node, start);
                }
                Verdict::Unknown(along_u) => {
                    if node.depth < MAX_DEPTH && splits < self.budget {
                        splits += 1;
                        self.split(&node, along_u);
                        continue;
                    } else if fallbacks < MAX_FALLBACKS {
                        fallbacks += 1;
                        let start = self.closest_control_point(node.edge);
                        self.minimise(&node, start);
                    }
                }
            }
            self.free.push(node.slot);
        }
        // The distance is flat to the resolution here (a point on an axis
        // of revolution): wherever else it is, other patches can only
        // match it, so they get a smaller budget.
        if splits >= self.budget && !self.all {
            self.budget = MAX_SPLITS / 8;
        }
    }

    /// A free net slot.
    fn slot(&mut self) -> usize {
        self.free.pop().unwrap_or_else(|| {
            let slot = self.nets.len();
            self.nets.resize(slot + self.size, 0.0);
            slot
        })
    }

    /// Queues the subpatch with the net at `slot` (see [`Self::node`]).
    fn push(&mut self, slot: usize, rect: [f64; 4], depth: u32, edge: Option<(usize, usize)>) {
        if let Some(node) = self.node(slot, rect, depth, edge) {
            self.pending.push(node);
        }
    }

    /// The subpatch with the net at `slot`, unless only the closest point
    /// is wanted and the subpatch is farther than the closest so far; its
    /// corners, on the surface, may become that point.
    fn node(
        &mut self,
        slot: usize,
        rect: [f64; 4],
        depth: u32,
        edge: Option<(usize, usize)>,
    ) -> Option<Subpatch> {
        let (rows, cols) = (self.rows, self.cols);
        if self.all {
            // Every subpatch counts, in any order.
            return Some(Subpatch {
                lower: 0.0,
                key: 0.0,
                slot,
                rect,
                depth,
                edge,
            });
        }
        math::net_points(&self.nets[slot..slot + self.size], self.h, &mut self.points);
        for (s, t) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
            let p = self.points[s * (rows - 1) * cols + t * (cols - 1)];
            let distance = Point::new(p[0], p[1], p[2]).distance(self.point);
            if distance < self.best.0 {
                // A collapsed net's corners sit on its edge.
                let [s, t] = match edge {
                    Some((0, side)) => [side, t],
                    Some((_, side)) => [s, side],
                    None => [s, t],
                };
                self.best = (distance, rect[s], rect[2 + t]);
            }
        }
        let q = [self.point.x, self.point.y, self.point.z];
        let lower = math::net_distance_bound(&self.points, rows, cols, q);
        if lower >= self.best.0 - self.resolution {
            self.free.push(slot);
            return None;
        }
        let key = self
            .points
            .iter()
            .map(|x| (0..3).map(|c| (x[c] - q[c]).powi(2)).sum::<f64>())
            .fold(f64::INFINITY, f64::min)
            .sqrt();
        Some(Subpatch {
            lower,
            key,
            slot,
            rect,
            depth,
            edge,
        })
    }

    /// Classifies `node` given its crease or boundary edges `creases`
    /// (`[u0, u1, v0, v1]`), leaving its Euclidean control points in
    /// `points`.
    fn verdict(&mut self, node: &Subpatch, creases: [bool; 4]) -> Verdict {
        let (rows, cols, h) = (self.rows, self.cols, self.h);
        let net = &self.nets[node.slot..node.slot + self.size];
        math::net_points(net, h, &mut self.points);
        let q = [self.point.x, self.point.y, self.point.z];
        // The side each direction's gradient pins the minimum to; a
        // collapsed net stays on its edge.
        let pinned = |gradient: [math::Interval; 2]| -> [Option<usize>; 2] {
            std::array::from_fn(|direction| match node.edge {
                Some((d, side)) if d == direction => Some(side),
                _ => {
                    let g = gradient[direction];
                    if g.lo > 0.0 {
                        Some(0)
                    } else if g.hi < 0.0 {
                        Some(1)
                    } else {
                        None
                    }
                }
            })
        };
        // Visiting every patch, most of which only the gradient rules out:
        // try it alone first.
        if self.all {
            let gradient =
                math::distance_bounds(net, &self.points, rows, cols, h, q, false).gradient;
            let smooth = |direction: usize, side: usize| !creases[2 * direction + side];
            if pinned(gradient)
                .iter()
                .enumerate()
                .any(|(direction, side)| side.is_some_and(|side| smooth(direction, side)))
            {
                return Verdict::None;
            }
        }
        let bounds = math::distance_bounds(net, &self.points, rows, cols, h, q, true);
        let [fss, fst, ftt] = bounds.hessian;
        let pinned = pinned(bounds.gradient);
        for direction in 0..2 {
            let crease = |side: usize| creases[2 * direction + side];
            // Descending through a smooth edge, or curving down towards
            // two of them: the minima are elsewhere.
            let concave = bounds.hessian[2 * direction].hi < 0.0;
            if pinned[direction].is_some_and(|side| !crease(side))
                || (concave && !crease(0) && !crease(1))
            {
                return Verdict::None;
            }
        }
        if !creases.contains(&true) && nowhere_convex([fss, fst, ftt]) {
            return Verdict::None;
        }
        match pinned {
            [Some(s), Some(t)] => Verdict::Corner([s, t]),
            [Some(s), None] if ftt.lo > 0.0 => Verdict::Edge(0, s),
            [None, Some(t)] if fss.lo > 0.0 => Verdict::Edge(1, t),
            // Along the edge, as a curve.
            [Some(s), None] if node.edge.is_none() => Verdict::Pinned(0, s),
            [None, Some(t)] if node.edge.is_none() => Verdict::Pinned(1, t),
            [Some(_), None] => Verdict::Unknown(false),
            [None, Some(_)] => Verdict::Unknown(true),
            [None, None]
                if fss.lo > 0.0 && ftt.lo > 0.0 && fss.lo * ftt.lo > fst.magnitude().powi(2) =>
            {
                Verdict::Convex
            }
            // Along the direction the patch stretches most.
            _ => Verdict::Unknown(bounds.stretch[0] >= bounds.stretch[1]),
        }
    }

    /// The parameters, as fractions `(s, t)` of the subpatch, of its
    /// control point closest to the point, among those on edge `side` of
    /// `direction` if given.
    fn closest_control_point(&self, edge: Option<(usize, usize)>) -> (f64, f64) {
        let (p, q) = (self.rows - 1, self.cols - 1);
        let target = [self.point.x, self.point.y, self.point.z];
        let distance = |i: usize, j: usize| {
            let x = self.points[i * self.cols + j];
            (0..3).map(|c| (x[c] - target[c]).powi(2)).sum::<f64>()
        };
        let (rows, cols) = match edge {
            Some((0, side)) => (side * p..=side * p, 0..=q),
            Some((_, side)) => (0..=p, side * q..=side * q),
            None => (0..=p, 0..=q),
        };
        let (i, j) = rows
            .flat_map(|i| cols.clone().map(move |j| (i, j)))
            .min_by(|a, b| distance(a.0, a.1).total_cmp(&distance(b.0, b.1)))
            .unwrap_or((0, 0));
        (i as f64 / p as f64, j as f64 / q as f64)
    }

    /// Newton minimisation within the subpatch from the fractions `(s, t)`.
    fn minimise(&mut self, node: &Subpatch, (s, t): (f64, f64)) {
        let r = node.rect;
        let (u, v) = self.surface.refine_projection(
            self.point,
            r[0] + s * (r[1] - r[0]),
            r[2] + t * (r[3] - r[2]),
            [Some((r[0], r[1])), Some((r[2], r[3]))],
        );
        let distance = self.surface.eval_point(u, v).distance(self.point);
        self.found(distance, u, v);
    }

    /// Records a minimum at `(u, v)`, `distance` from the point.
    fn found(&mut self, distance: f64, u: f64, v: f64) {
        if distance < self.best.0 {
            self.best = (distance, u, v);
        }
        if self.all {
            self.minima.push((u, v));
        }
    }

    /// Halves the subpatch along `u` (else `v`) and queues both halves,
    /// the first taking over its slot.
    fn split(&mut self, node: &Subpatch, along_u: bool) {
        self.scratch
            .copy_from_slice(&self.nets[node.slot..node.slot + self.size]);
        let (left, right) = (node.slot, self.slot());
        let (lo, hi) = if left < right {
            let (a, b) = self.nets.split_at_mut(right);
            (&mut a[left..left + self.size], &mut b[..self.size])
        } else {
            let (a, b) = self.nets.split_at_mut(left);
            (&mut b[..self.size], &mut a[right..right + self.size])
        };
        math::halve_net(&self.scratch, self.rows, self.cols, self.h, along_u, lo, hi);
        let r = node.rect;
        let (first, second) = if along_u {
            let mid = 0.5 * (r[0] + r[1]);
            ([r[0], mid, r[2], r[3]], [mid, r[1], r[2], r[3]])
        } else {
            let mid = 0.5 * (r[2] + r[3]);
            ([r[0], r[1], r[2], mid], [r[0], r[1], mid, r[3]])
        };
        self.push(left, first, node.depth + 1, node.edge);
        self.push(right, second, node.depth + 1, node.edge);
    }
}

/// Whether none of the symmetric matrices `[[a, b], [b, c]]` with entries
/// in `[a, b, c]` is positive semidefinite, as the Hessian at a local
/// minimum is: their determinant stays negative, or they all curve down
/// along one direction (the eigenvector of the central matrix's lower
/// eigenvalue).
fn nowhere_convex([a, b, c]: [math::Interval; 3]) -> bool {
    let largest_ac = (a.lo * c.lo)
        .max(a.lo * c.hi)
        .max(a.hi * c.lo)
        .max(a.hi * c.hi);
    let smallest_b = if b.lo > 0.0 {
        b.lo
    } else if b.hi < 0.0 {
        -b.hi
    } else {
        0.0
    };
    if largest_ac < smallest_b * smallest_b {
        return true;
    }
    let (am, bm, cm) = (a.mid(), b.mid(), c.mid());
    let lower = 0.5 * (am + cm) - (0.25 * (am - cm).powi(2) + bm * bm).sqrt();
    let (x, y) = if (am - lower).abs() > (cm - lower).abs() {
        (bm, lower - am)
    } else {
        (lower - cm, bm)
    };
    let along = if x * y >= 0.0 { b.hi } else { b.lo };
    (x != 0.0 || y != 0.0) && x * x * a.hi + 2.0 * x * y * along + y * y * c.hi < 0.0
}

/// The values at `indices`, in order.
fn gathered<T: Copy>(values: &[T], indices: &[usize]) -> Vec<T> {
    indices.iter().map(|&i| values[i]).collect()
}

#[cfg(test)]
mod tests;
