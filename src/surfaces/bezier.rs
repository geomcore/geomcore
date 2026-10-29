//! Bézier surfaces in 3D: the polynomial (or rational) tensor-product
//! patch of one span in each direction, given by a
//! `(u_degree + 1) x (v_degree + 1)` grid of control poles.
//!
//! A [`BezierSurface`] is a clamped, non-periodic
//! [`BSplineSurface`] of a single knot span in `u` and in `v`, held by that
//! type's engine as [`BezierCurve`](crate::BezierCurve) is held by the
//! curve one: evaluation, projection and editing are the B-spline's, so
//! the two cannot disagree. What a Bézier surface does not expose is the
//! knot vectors: it has no interior knots, no periodicity and no knot
//! editing. [`BezierSurface::to_bspline`] gives the B-spline surface when
//! those are wanted.

use crate::bspline::curve::MAX_DEGREE;
use crate::curves::{BSplineSplitError, BezierConstructionError, DegreeElevationError, Side};
use crate::surfaces::{BSplineSurface, ParametricSurface};
use crate::{Continuity, Point, SurfaceProjection, Tolerance, Vector};

/// A Bézier surface (optionally rational) in 3D: a grid of
/// `(u_degree + 1) x (v_degree + 1)` control poles over one parameter
/// rectangle, `[0, 1] x [0, 1]` when built from poles.
///
/// The pole grid `poles[i_u][i_v]` is laid out as that of a
/// [`BSplineSurface`]: `i_u` indexes the `u` direction, a row per `u`
/// pole, and [`BezierSurface::poles`] lists it row-major. The surface
/// passes through the four corner poles and lies in the convex hull of
/// its poles. A rational surface carries a positive weight per pole and is
/// evaluated in homogeneous coordinates. The patches of a B-spline surface
/// ([`BSplineSurface::to_bezier_patches`]) keep the spans they came from
/// as their bounds.
///
/// # Examples
///
/// A biquadratic patch, `(2u, 2v, 4uv)`:
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{BezierSurface, Point};
///
/// let poles = (0..3)
///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
///     .collect();
/// let surface = BezierSurface::new(2, 2, poles).unwrap();
///
/// assert_eq!(surface.eval_point(1.0, 1.0), Point::new(2.0, 2.0, 4.0));
/// let p = surface.eval_point(0.25, 0.5);
/// assert!(p.distance(Point::new(0.5, 1.0, 0.5)) < 1e-14);
/// assert_eq!(surface.u_bounds(), (0.0, 1.0));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BezierSurface(BSplineSurface);

impl BezierSurface {
    /// Creates a non-rational (polynomial) Bézier surface of the given
    /// degrees over `[0, 1] x [0, 1]`.
    ///
    /// `poles[i_u][i_v]` is the control net: `u_degree + 1` rows of
    /// `v_degree + 1` poles, as for [`BSplineSurface::new`].
    ///
    /// # Errors
    ///
    /// [`BezierConstructionError::InvalidDegree`] if either degree is not
    /// in `1..=25`, and [`BezierConstructionError::PoleCountMismatch`] if
    /// the grid is not `u_degree + 1` rows of `v_degree + 1` poles each
    /// (an empty or ragged grid included).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierConstructionError, BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// let surface = BezierSurface::new(1, 1, poles.clone()).unwrap();
    /// assert_eq!(surface.eval_point(0.5, 0.5), Point::new(0.5, 0.5, 0.25));
    /// assert_eq!(
    ///     BezierSurface::new(2, 1, poles),
    ///     Err(BezierConstructionError::PoleCountMismatch)
    /// );
    /// ```
    pub fn new(
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Vec<Point>>,
    ) -> Result<BezierSurface, BezierConstructionError> {
        Self::build(u_degree, v_degree, poles, None)
    }

    /// Creates a rational Bézier surface of the given degrees over
    /// `[0, 1] x [0, 1]` with one weight per pole, laid out
    /// `weights[i_u][i_v]` matching `poles`.
    ///
    /// # Errors
    ///
    /// [`BezierConstructionError::WeightCountMismatch`] if the weight grid
    /// does not have the shape of the pole grid,
    /// [`BezierConstructionError::NonPositiveWeight`] if a weight is
    /// `<= 0`; otherwise as [`BezierSurface::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierSurface, Point};
    ///
    /// // A rational quadratic quarter circle in `u`, a line in `v`.
    /// let w = std::f64::consts::FRAC_1_SQRT_2;
    /// let row = |p: Point| vec![p, Point::new(p.x, p.y, 1.0)];
    /// let poles = vec![
    ///     row(Point::new(1.0, 0.0, 0.0)),
    ///     row(Point::new(1.0, 1.0, 0.0)),
    ///     row(Point::new(0.0, 1.0, 0.0)),
    /// ];
    /// let weights = vec![vec![1.0; 2], vec![w; 2], vec![1.0; 2]];
    /// let surface = BezierSurface::new_rational(2, 1, poles, weights).unwrap();
    /// let p = surface.eval_point(1.0, 0.5);
    /// assert!(p.distance(Point::new(0.0, 1.0, 0.5)) < 1e-12);
    /// ```
    pub fn new_rational(
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Vec<Point>>,
        weights: Vec<Vec<f64>>,
    ) -> Result<BezierSurface, BezierConstructionError> {
        Self::build(u_degree, v_degree, poles, Some(weights))
    }

    /// The weights are checked before the degrees and the pole grid, as
    /// for [`BezierCurve`](crate::BezierCurve).
    fn build(
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Vec<Point>>,
        weights: Option<Vec<Vec<f64>>>,
    ) -> Result<BezierSurface, BezierConstructionError> {
        if let Some(ws) = &weights {
            let same_shape =
                ws.len() == poles.len() && ws.iter().zip(&poles).all(|(w, p)| w.len() == p.len());
            if !same_shape {
                return Err(BezierConstructionError::WeightCountMismatch);
            }
            if ws.iter().flatten().any(|&w| w <= 0.0) {
                return Err(BezierConstructionError::NonPositiveWeight);
            }
        }
        let degrees = 1..=MAX_DEGREE;
        if !degrees.contains(&u_degree) || !degrees.contains(&v_degree) {
            return Err(BezierConstructionError::InvalidDegree);
        }
        let shaped =
            poles.len() == u_degree + 1 && poles.iter().all(|row| row.len() == v_degree + 1);
        if !shaped {
            return Err(BezierConstructionError::PoleCountMismatch);
        }
        let unit = (0.0, 1.0);
        Ok(BezierSurface(BSplineSurface::bezier_patch(
            (u_degree, unit),
            (v_degree, unit),
            poles.concat(),
            weights.map(|w| w.concat()),
        )))
    }

    /// The Bézier surface a B-spline surface is: one clamped,
    /// non-periodic span in each direction, keeping its bounds, poles and
    /// weights.
    ///
    /// # Errors
    ///
    /// [`BezierConstructionError::NotABezier`] if either direction has
    /// interior knots, is periodic, or has end multiplicities other than
    /// `degree + 1`. [`BSplineSurface::insert_u_knot`] and its siblings
    /// make such a surface; [`BSplineSurface::to_bezier_patches`] cuts one
    /// into Bézier surfaces.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineSurface, BezierConstructionError, BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// let spline = BSplineSurface::new(
    ///     1, 1, poles, vec![2.0, 5.0], vec![2, 2], vec![0.0, 1.0], vec![2, 2], false, false,
    /// )
    /// .unwrap();
    /// let surface = BezierSurface::from_bspline(&spline).unwrap();
    /// assert_eq!(surface.u_bounds(), (2.0, 5.0));
    /// assert_eq!(surface.to_bspline(), spline);
    ///
    /// let refined = spline.insert_u_knot(3.0, 1).unwrap();
    /// assert_eq!(
    ///     BezierSurface::from_bspline(&refined),
    ///     Err(BezierConstructionError::NotABezier)
    /// );
    /// ```
    pub fn from_bspline(
        surface: &BSplineSurface,
    ) -> Result<BezierSurface, BezierConstructionError> {
        if surface.is_bezier() {
            Ok(BezierSurface(surface.clone()))
        } else {
            Err(BezierConstructionError::NotABezier)
        }
    }

    /// The surface as a B-spline surface: degrees, poles and weights
    /// unchanged, over the same bounds with end knots of multiplicity
    /// `degree + 1` in both directions.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// let spline = BezierSurface::new(1, 1, poles).unwrap().to_bspline();
    /// assert_eq!(spline.u_knots(), &[0.0, 1.0]);
    /// assert_eq!(spline.v_multiplicities(), &[2, 2]);
    /// ```
    pub fn to_bspline(&self) -> BSplineSurface {
        self.0.clone()
    }

    /// Returns the surface's degree in the `u` direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// assert_eq!(surface.u_degree(), 1);
    /// ```
    pub fn u_degree(&self) -> usize {
        self.0.u_degree()
    }

    /// Returns the surface's degree in the `v` direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// assert_eq!(surface.v_degree(), 1);
    /// ```
    pub fn v_degree(&self) -> usize {
        self.0.v_degree()
    }

    /// Returns the number of pole rows, `u_degree + 1`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// assert_eq!(surface.u_pole_count(), 2);
    /// ```
    pub fn u_pole_count(&self) -> usize {
        self.0.u_pole_count()
    }

    /// Returns the number of poles in a row, `v_degree + 1`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// assert_eq!(surface.v_pole_count(), 2);
    /// ```
    pub fn v_pole_count(&self) -> usize {
        self.0.v_pole_count()
    }

    /// Returns whether the surface is rational (has per-pole weights).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// assert!(!surface.is_rational());
    /// ```
    pub fn is_rational(&self) -> bool {
        self.0.is_rational()
    }

    /// Returns the control poles, row-major: a row of `v_pole_count`
    /// poles per `u` pole, `poles()[i_u * v_pole_count + i_v]`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// let surface = BezierSurface::new(1, 1, poles.clone()).unwrap();
    /// assert_eq!(surface.poles(), poles.concat().as_slice());
    /// ```
    pub fn poles(&self) -> &[Point] {
        self.0.poles()
    }

    /// Returns the per-pole weights, laid out as
    /// [`poles`](Self::poles), or `None` if the surface is not rational.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// assert_eq!(surface.weights(), None);
    /// ```
    pub fn weights(&self) -> Option<&[f64]> {
        self.0.weights()
    }

    /// Returns the piece of the surface either side of `u`, over
    /// `[first, u]` and `[u, last]` of its [`u_bounds`](Self::u_bounds):
    /// Bézier surfaces of the same degrees, evaluating as this surface
    /// over their ranges and meeting along the `u` isoparametric curve at
    /// `u`, by de Casteljau.
    ///
    /// # Errors
    ///
    /// [`BSplineSplitError::NonFinite`] if `u` is NaN or infinite,
    /// [`BSplineSplitError::OutOfRange`] if it lies outside the bounds,
    /// and [`BSplineSplitError::EmptyRange`] if it is one of them.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BezierSurface::new(2, 1, poles).unwrap();
    /// let (low, high) = surface.split_u_at(0.5).unwrap();
    /// assert_eq!((low.u_bounds(), high.u_bounds()), ((0.0, 0.5), (0.5, 1.0)));
    /// assert_eq!(&low.poles()[4..], &high.poles()[..2]);
    /// assert!(high.eval_point(0.7, 0.3).distance(surface.eval_point(0.7, 0.3)) < 1e-12);
    /// ```
    pub fn split_u_at(&self, u: f64) -> Result<(BezierSurface, BezierSurface), BSplineSplitError> {
        self.0
            .split_u_at(u)
            .map(|(low, high)| (Self::from_surface(low), Self::from_surface(high)))
    }

    /// Returns the piece of the surface either side of `v`, as
    /// [`split_u_at`](Self::split_u_at) for `v`.
    ///
    /// # Errors
    ///
    /// As [`split_u_at`](Self::split_u_at), against the `v` bounds.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = (0..2)
    ///     .map(|i| (0..3).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BezierSurface::new(1, 2, poles).unwrap();
    /// let (low, high) = surface.split_v_at(0.25).unwrap();
    /// assert_eq!((low.v_bounds(), high.v_bounds()), ((0.0, 0.25), (0.25, 1.0)));
    /// ```
    pub fn split_v_at(&self, v: f64) -> Result<(BezierSurface, BezierSurface), BSplineSplitError> {
        self.0
            .split_v_at(v)
            .map(|(low, high)| (Self::from_surface(low), Self::from_surface(high)))
    }

    /// Returns the same surface with its `u` degree raised by `amount`: it
    /// evaluates to the same points and derivatives everywhere, up to
    /// rounding, with `amount` more pole rows. A rational surface is
    /// raised in homogeneous coordinates and stays rational. An `amount`
    /// of 0 returns an identical copy.
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
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = (0..2)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let linear = BezierSurface::new(1, 1, poles).unwrap();
    /// let quadratic = linear.elevate_u_degree(1).unwrap();
    /// assert_eq!((quadratic.u_degree(), quadratic.u_pole_count()), (2, 3));
    /// assert!(quadratic.eval_point(0.3, 0.6).distance(linear.eval_point(0.3, 0.6)) < 1e-12);
    /// ```
    pub fn elevate_u_degree(&self, amount: usize) -> Result<BezierSurface, DegreeElevationError> {
        self.0.elevate_u_degree(amount).map(Self::from_surface)
    }

    /// Returns the same surface with its `v` degree raised by `amount`, as
    /// [`elevate_u_degree`](Self::elevate_u_degree) for `v`.
    ///
    /// # Errors
    ///
    /// As [`elevate_u_degree`](Self::elevate_u_degree).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// let quadratic = surface.elevate_v_degree(1).unwrap();
    /// assert_eq!((quadratic.v_degree(), quadratic.v_pole_count()), (2, 3));
    /// assert!(quadratic.eval_point(0.3, 0.6).distance(surface.eval_point(0.3, 0.6)) < 1e-12);
    /// ```
    pub fn elevate_v_degree(&self, amount: usize) -> Result<BezierSurface, DegreeElevationError> {
        self.0.elevate_v_degree(amount).map(Self::from_surface)
    }

    /// Returns the surface with its `u` direction traversed the other way:
    /// its point at `(u, v)` is this surface's point at
    /// `(first + last - u, v)` over the same bounds, so `Su` changes sign.
    /// The rows of the net are listed in reverse order.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = (0..3)
    ///     .map(|i| (0..2).map(|j| Point::new(i as f64, j as f64, (i * j) as f64)).collect())
    ///     .collect();
    /// let surface = BezierSurface::new(2, 1, poles).unwrap();
    /// let reversed = surface.reversed_u();
    /// assert_eq!(reversed.eval_point(0.0, 0.5), surface.eval_point(1.0, 0.5));
    /// assert_eq!(reversed.reversed_u(), surface);
    /// ```
    pub fn reversed_u(&self) -> BezierSurface {
        Self::from_surface(self.0.reversed_u())
    }

    /// Returns the surface with its `v` direction traversed the other way,
    /// as [`reversed_u`](Self::reversed_u) for `v`; `Sv` changes sign.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierSurface, Point};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// let reversed = surface.reversed_v();
    /// assert_eq!(reversed.eval_point(0.5, 0.0), surface.eval_point(0.5, 1.0));
    /// assert_eq!(reversed.reversed_v(), surface);
    /// ```
    pub fn reversed_v(&self) -> BezierSurface {
        Self::from_surface(self.0.reversed_v())
    }

    /// Evaluates the point and the first derivatives at `(u, v)`, each
    /// direction from its own side, as [`BSplineSurface::eval_d1_sided`].
    /// A Bézier surface is one polynomial piece, so both sides agree
    /// everywhere.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierSurface, Point, Side};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// assert_eq!(
    ///     surface.eval_d1_sided(0.5, 0.5, Side::Below, Side::Below),
    ///     surface.eval_d1_sided(0.5, 0.5, Side::Above, Side::Above)
    /// );
    /// ```
    pub fn eval_d1_sided(
        &self,
        u: f64,
        v: f64,
        u_side: Side,
        v_side: Side,
    ) -> (Point, [Vector; 2]) {
        self.0.eval_d1_sided(u, v, u_side, v_side)
    }

    /// Evaluates the derivative of order `(du, dv)` at `(u, v)` from the
    /// given sides, as [`eval_d1_sided`](Self::eval_d1_sided).
    ///
    /// # Panics
    ///
    /// As [`eval_derivative`](Self::eval_derivative).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierSurface, Point, Side, Vector};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// let suv = surface.eval_derivative_sided(0.5, 0.5, 1, 1, Side::Above, Side::Below);
    /// assert_eq!(suv, Vector::new(0.0, 0.0, 1.0));
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
        self.0.eval_derivative_sided(u, v, du, dv, u_side, v_side)
    }

    /// Evaluates the point, the first and the second derivatives at
    /// `(u, v)` from the given sides, as
    /// [`eval_d1_sided`](Self::eval_d1_sided).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierSurface, Point, Side};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// assert_eq!(
    ///     surface.eval_d2_sided(0.5, 0.5, Side::Below, Side::Above),
    ///     surface.eval_d2(0.5, 0.5)
    /// );
    /// ```
    pub fn eval_d2_sided(
        &self,
        u: f64,
        v: f64,
        u_side: Side,
        v_side: Side,
    ) -> (Point, [Vector; 2], [Vector; 3]) {
        self.0.eval_d2_sided(u, v, u_side, v_side)
    }

    /// The local minima of the distance from `point` to the surface,
    /// closest first, as [`BSplineSurface::extrema`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierSurface, Point, Tolerance};
    ///
    /// let poles = vec![
    ///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
    ///     vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    /// ];
    /// // The surface is (u, v, uv).
    /// let surface = BezierSurface::new(1, 1, poles).unwrap();
    /// let extrema = surface.extrema(Point::new(0.5, 0.5, 0.25), Tolerance::DEFAULT);
    /// assert_eq!(extrema[0].distance, 0.0);
    /// ```
    pub fn extrema(&self, point: Point, tol: Tolerance) -> Vec<SurfaceProjection> {
        self.0.extrema(point, tol)
    }

    /// The surface wrapping `patch`, a single clamped span in both
    /// directions (as [`BSplineSurface::to_bezier_patches`] makes).
    pub(crate) fn from_surface(patch: BSplineSurface) -> BezierSurface {
        debug_assert!(patch.is_bezier());
        BezierSurface(patch)
    }
}

impl ParametricSurface for BezierSurface {
    fn eval_point(&self, u: f64, v: f64) -> Point {
        self.0.eval_point(u, v)
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        self.0.project_point(point, tol)
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        self.0.eval_derivative(u, v, du, dv)
    }

    fn eval_d1(&self, u: f64, v: f64) -> (Point, [Vector; 2]) {
        self.0.eval_d1(u, v)
    }

    fn eval_d2(&self, u: f64, v: f64) -> (Point, [Vector; 2], [Vector; 3]) {
        self.0.eval_d2(u, v)
    }

    fn u_bounds(&self) -> (f64, f64) {
        self.0.u_bounds()
    }

    fn v_bounds(&self) -> (f64, f64) {
        self.0.v_bounds()
    }

    fn u_period(&self) -> Option<f64> {
        None
    }

    fn v_period(&self) -> Option<f64> {
        None
    }

    fn continuity_u(&self, _u: f64, _tol: Tolerance) -> Continuity {
        Continuity::Smooth
    }

    fn continuity_v(&self, _v: f64, _tol: Tolerance) -> Continuity {
        Continuity::Smooth
    }

    fn is_u_closed(&self, tol: Tolerance) -> bool {
        self.0.is_u_closed(tol)
    }

    fn is_v_closed(&self, tol: Tolerance) -> bool {
        self.0.is_v_closed(tol)
    }

    fn as_bspline(&self) -> Option<&BSplineSurface> {
        Some(&self.0)
    }
}
