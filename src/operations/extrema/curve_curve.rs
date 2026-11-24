//! The local minima of the distance between two curves
//! ([`CurveCurveExtrema`] in space, [`CurveCurveExtrema2D`] in the plane).

use crate::Direction;
use crate::base::TierPoint;
use crate::geometry::curves::ParametricCurve;
use crate::math::numeric::{newton_min_1d, newton_min_2d};
use crate::math::parameter_domain::ParameterDomain;
use crate::operations::extrema::curve_curve_pieces::curve_curve_seeds;
use crate::operations::extrema::operand::CurveExtremaOperand;
use crate::operations::extrema::{Extremum, parameters_apart};
use crate::operations::intersect::{Completeness, CurveCurveRegion};
use crate::operations::search_domain::curve_domain;
use crate::{Point, Point2D, Tolerance, Vector};
use std::ops::Deref;

/// One local minimum of the distance between two curves of one tier:
/// [`CurveCurveExtremum`] in space, [`CurveCurveExtremum2D`] in the plane
/// (see [`curve_curve_extrema`](crate::curve_curve_extrema)).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenericCurveCurveExtremum<P: TierPoint> {
    /// Parameter on the first curve.
    pub first_parameter: f64,
    /// Parameter on the second curve.
    pub second_parameter: f64,
    /// The point of the first curve at `first_parameter`.
    pub first_point: P,
    /// The point of the second curve at `second_parameter`.
    pub second_point: P,
    /// Distance between the two points.
    pub distance: f64,
}

/// One local minimum of the distance between two 3D curves (see
/// [`curve_curve_extrema`](crate::curve_curve_extrema)).
///
/// # Examples
///
/// ```
/// use geomcore::{Line, Point, Tolerance, Vector, curve_curve_extrema};
/// let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
/// let closest = curve_curve_extrema(&l1, &l2, Tolerance::DEFAULT)[0];
/// assert_eq!(closest.second_point, Point::new(0.0, 0.0, 1.0));
/// ```
pub type CurveCurveExtremum = GenericCurveCurveExtremum<Point>;

/// One local minimum of the distance between two 2D curves (see
/// [`curve_curve_extrema`](crate::curve_curve_extrema)).
///
/// # Examples
///
/// ```
/// use geomcore::{Line2D, Point2D, Tolerance, Vector2D, curve_curve_extrema};
/// let l1 = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// let l2 = Line2D::new(Point2D::new(0.0, 1.0), Vector2D::X).unwrap();
/// let closest = curve_curve_extrema(&l1, &l2, Tolerance::DEFAULT)[0];
/// assert_eq!(closest.second_point, Point2D::new(0.0, 1.0));
/// ```
pub type CurveCurveExtremum2D = GenericCurveCurveExtremum<Point2D>;

impl<P: TierPoint> GenericCurveCurveExtremum<P> {
    /// The extremum between `first_point`, at `first_parameter` on the
    /// first curve, and `second_point`, at `second_parameter` on the
    /// second.
    pub(crate) fn new(
        first_parameter: f64,
        first_point: P,
        second_parameter: f64,
        second_point: P,
    ) -> Self {
        GenericCurveCurveExtremum {
            first_parameter,
            second_parameter,
            first_point,
            second_point,
            distance: first_point.distance(second_point),
        }
    }

    /// The same extremum with the curves' roles exchanged.
    pub(crate) fn swapped(self) -> Self {
        GenericCurveCurveExtremum::new(
            self.second_parameter,
            self.second_point,
            self.first_parameter,
            self.first_point,
        )
    }
}

impl Extremum for CurveCurveExtremum {
    fn distance(&self) -> f64 {
        self.distance
    }

    fn parameters_apart(&self, other: &Self, tol: Tolerance) -> bool {
        parameters_apart(self.first_parameter, other.first_parameter, tol)
            || parameters_apart(self.second_parameter, other.second_parameter, tol)
    }
}

/// The local minima of the distance between two curves of one tier,
/// closest first, and whether they are all of them:
/// [`CurveCurveExtrema`] in space, [`CurveCurveExtrema2D`] in the plane
/// (see [`curve_curve_extrema`](crate::curve_curve_extrema)).
///
/// It reads as the slice of its minima.
#[derive(Debug, Clone, PartialEq)]
pub struct GenericCurveCurveExtrema<P: TierPoint> {
    minima: Vec<GenericCurveCurveExtremum<P>>,
    completeness: Completeness<CurveCurveRegion>,
}

/// The local minima of the distance between two 3D curves, closest
/// first, and whether they are all of them (see
/// [`curve_curve_extrema`](crate::curve_curve_extrema)).
///
/// It reads as the slice of its minima.
///
/// # Examples
///
/// ```
/// use geomcore::{BSplineCurve, Point, Tolerance, curve_curve_extrema};
/// // An arch peaking at (1, 1, 0), under a rail across it at height 1.
/// let arch = vec![
///     Point::new(0.0, 0.0, 0.0),
///     Point::new(1.0, 2.0, 0.0),
///     Point::new(2.0, 0.0, 0.0),
/// ];
/// let arch = BSplineCurve::new(2, arch, vec![0.0, 1.0], vec![3, 3], false).unwrap();
/// let rail = vec![Point::new(1.0, 0.0, 1.0), Point::new(1.0, 2.0, 1.0)];
/// let rail = BSplineCurve::new(1, rail, vec![0.0, 1.0], vec![2, 2], false).unwrap();
/// let minima = curve_curve_extrema(&arch, &rail, Tolerance::DEFAULT);
/// // Both are B-splines: the search proves there is no other minimum.
/// assert!(minima.is_complete());
/// assert_eq!(minima.len(), 1);
/// assert!((minima[0].distance - 1.0).abs() < 1e-12);
/// for minimum in &minima {
///     assert!((minimum.first_parameter - 0.5).abs() < 1e-9);
/// }
/// ```
pub type CurveCurveExtrema = GenericCurveCurveExtrema<Point>;

/// The local minima of the distance between two 2D curves, closest
/// first, and whether they are all of them (see
/// [`curve_curve_extrema`](crate::curve_curve_extrema)).
///
/// It reads as the slice of its minima.
///
/// # Examples
///
/// ```
/// use geomcore::{Circle2D, Point2D, Tolerance, curve_curve_extrema};
/// // Two disjoint circles: the gap between them is the closest minimum.
/// let first = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
/// let second = Circle2D::new(Point2D::new(5.0, 0.0), 2.0).unwrap();
/// let minima = curve_curve_extrema(&first, &second, Tolerance::DEFAULT);
/// let closest = minima.closest().unwrap();
/// assert!((closest.distance - 2.0).abs() < 1e-12);
/// assert!(closest.first_point.distance(Point2D::new(1.0, 0.0)) < 1e-9);
/// assert!(closest.second_point.distance(Point2D::new(3.0, 0.0)) < 1e-9);
/// ```
pub type CurveCurveExtrema2D = GenericCurveCurveExtrema<Point2D>;

impl<P: TierPoint> GenericCurveCurveExtrema<P> {
    /// `minima`, put in order of distance, with their completeness.
    pub(crate) fn new(
        mut minima: Vec<GenericCurveCurveExtremum<P>>,
        completeness: Completeness<CurveCurveRegion>,
    ) -> Self {
        minima.sort_by(|x, y| x.distance.total_cmp(&y.distance));
        GenericCurveCurveExtrema {
            minima,
            completeness,
        }
    }

    /// The extrema made of the parts their getters return: for rebuilding
    /// stored extrema, as Python's `pickle` does. The parts are kept as
    /// they are, not checked or sorted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{CurveCurveExtrema, Line, Point, Tolerance, Vector, curve_curve_extrema};
    /// let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    /// let extrema = curve_curve_extrema(&l1, &l2, Tolerance::DEFAULT);
    /// let rebuilt =
    ///     CurveCurveExtrema::from_parts(extrema.minima().to_vec(), extrema.completeness().clone());
    /// assert_eq!(rebuilt, extrema);
    /// ```
    pub fn from_parts(
        minima: Vec<GenericCurveCurveExtremum<P>>,
        completeness: Completeness<CurveCurveRegion>,
    ) -> Self {
        GenericCurveCurveExtrema {
            minima,
            completeness,
        }
    }

    /// The same extrema with each point mapped by `point`: the extrema of
    /// two 2D curves read back from the search on them embedded in space.
    pub(crate) fn map_points<Q: TierPoint>(
        self,
        point: impl Fn(P) -> Q,
    ) -> GenericCurveCurveExtrema<Q> {
        let minima = self
            .minima
            .into_iter()
            .map(|minimum| GenericCurveCurveExtremum {
                first_parameter: minimum.first_parameter,
                second_parameter: minimum.second_parameter,
                first_point: point(minimum.first_point),
                second_point: point(minimum.second_point),
                distance: minimum.distance,
            })
            .collect();
        GenericCurveCurveExtrema {
            minima,
            completeness: self.completeness,
        }
    }

    /// Minima found from samples.
    pub(crate) fn sampled(minima: Vec<GenericCurveCurveExtremum<P>>) -> Self {
        Self::new(minima, Completeness::Sampled)
    }

    /// The minima, closest first.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Tolerance, Vector, curve_curve_extrema};
    /// let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    /// let minima = curve_curve_extrema(&l1, &l2, Tolerance::DEFAULT);
    /// assert_eq!(minima.minima()[0].distance, 1.0);
    /// ```
    pub fn minima(&self) -> &[GenericCurveCurveExtremum<P>] {
        &self.minima
    }

    /// The closest minimum, the first of [`minima`](Self::minima), or
    /// `None` when there is none.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Tolerance, Vector, curve_curve_extrema};
    /// let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    /// let closest = curve_curve_extrema(&l1, &l2, Tolerance::DEFAULT).closest().copied().unwrap();
    /// assert_eq!(closest.first_point, Point::ORIGIN);
    /// assert_eq!(closest.second_point, Point::new(0.0, 0.0, 1.0));
    /// assert_eq!(closest.distance, 1.0);
    /// ```
    pub fn closest(&self) -> Option<&GenericCurveCurveExtremum<P>> {
        self.minima.first()
    }

    /// The minima, closest first, without the completeness.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Tolerance, Vector, curve_curve_extrema};
    /// let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    /// let minima: Vec<_> = curve_curve_extrema(&l1, &l2, Tolerance::DEFAULT).into_minima();
    /// assert_eq!(minima.len(), 1);
    /// ```
    pub fn into_minima(self) -> Vec<GenericCurveCurveExtremum<P>> {
        self.minima
    }

    /// Whether the minima are all of them, and if not, why and where
    /// minima may be missing.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Completeness, Line, Point, Tolerance, Vector, curve_curve_extrema};
    /// let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    /// let minima = curve_curve_extrema(&l1, &l2, Tolerance::DEFAULT);
    /// assert_eq!(minima.completeness(), &Completeness::Sampled);
    /// ```
    pub fn completeness(&self) -> &Completeness<CurveCurveRegion> {
        &self.completeness
    }

    /// Whether the minima are all of them ([`Completeness::Complete`]).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Tolerance, Vector, curve_curve_extrema};
    /// let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    /// // Only a search with a B-spline on one side proves its minima complete.
    /// assert!(!curve_curve_extrema(&l1, &l2, Tolerance::DEFAULT).is_complete());
    /// ```
    pub fn is_complete(&self) -> bool {
        self.completeness.is_complete()
    }
}

impl<P: TierPoint> Deref for GenericCurveCurveExtrema<P> {
    type Target = [GenericCurveCurveExtremum<P>];

    fn deref(&self) -> &[GenericCurveCurveExtremum<P>] {
        &self.minima
    }
}

impl<P: TierPoint> IntoIterator for GenericCurveCurveExtrema<P> {
    type Item = GenericCurveCurveExtremum<P>;
    type IntoIter = std::vec::IntoIter<GenericCurveCurveExtremum<P>>;

    fn into_iter(self) -> std::vec::IntoIter<GenericCurveCurveExtremum<P>> {
        self.minima.into_iter()
    }
}

impl<'a, P: TierPoint> IntoIterator for &'a GenericCurveCurveExtrema<P> {
    type Item = &'a GenericCurveCurveExtremum<P>;
    type IntoIter = std::slice::Iter<'a, GenericCurveCurveExtremum<P>>;

    fn into_iter(self) -> std::slice::Iter<'a, GenericCurveCurveExtremum<P>> {
        self.minima.iter()
    }
}

/// Grid points per finite parameter domain in the 2D extrema search.
const SEED_GRID: usize = 32;
/// Samples per finite parameter domain in the 1D (line against curve)
/// extrema search.
const SEED_GRID_1D: usize = 256;
/// Most seeds taken from one sampled grid (or edge) of the extrema
/// search.
pub(crate) const MAX_SEEDS: usize = 64;
/// Most halvings of span pairs in the search for the points where two
/// curves may meet (see [`meeting_seeds`]).
const MAX_MEETING_HALVINGS: usize = 256;
/// Unbounded parameters are seeded until the curve has moved this many
/// configuration scales from its reference point.
pub(crate) const SEED_REACH: f64 = 64.0;

impl ParameterDomain {
    /// Seed parameters, ascending. A finite domain gets `count` uniform
    /// samples. An unbounded direction gets a geometric progression away
    /// from the reference parameter, starting at a step that moves the
    /// curve `scale / 16` and growing by a quarter per sample until the
    /// curve is [`SEED_REACH`] scales from its reference point: dense
    /// where the curves interact, sparse far out, at any absolute size.
    pub(crate) fn seeds<C: ParametricCurve>(self, curve: &C, scale: f64, count: usize) -> Vec<f64> {
        match self {
            ParameterDomain::Periodic { start, period } => (0..count)
                .map(|i| start + period * i as f64 / count as f64)
                .collect(),
            // The ends exactly, so a minimum held there is seeded on it.
            ParameterDomain::Interval(lo, hi) if self.is_finite() => (0..count)
                .map(|i| match i {
                    0 => lo,
                    _ if i == count - 1 => hi,
                    _ => lo + (hi - lo) * i as f64 / (count - 1) as f64,
                })
                .collect(),
            ParameterDomain::Interval(lo, hi) => {
                let u0 = self.reference();
                let (p0, d1) = curve.first_derivatives(u0);
                let speed = d1.magnitude();
                let first = scale / (16.0 * speed);
                let first = if first.is_finite() && first > 0.0 {
                    first
                } else {
                    1.0
                };
                let mut seeds: Vec<f64> = self.finite_ends().collect();
                seeds.push(u0);
                for side in [-1.0, 1.0] {
                    let mut step = first;
                    for _ in 0..96 {
                        let u = u0 + side * step;
                        if u < lo || u > hi {
                            break;
                        }
                        seeds.push(u);
                        if curve.point_at(u).distance(p0) >= SEED_REACH * scale {
                            break;
                        }
                        step *= 1.25;
                    }
                }
                seeds.sort_by(f64::total_cmp);
                seeds.dedup();
                seeds
            }
        }
    }
}

/// A straight curve with constant speed: `origin + (u - u0) * dir`.
#[derive(Debug, Clone, Copy)]
struct Straight {
    u0: f64,
    origin: Point,
    dir: Vector,
}

impl Straight {
    /// `curve` as a straight curve, if it is one over its whole domain
    /// (lines, straight segments): checked against the affine extension
    /// of its reference point and tangent at a spread of parameters.
    fn of<C: ParametricCurve>(curve: &C, dom: ParameterDomain) -> Option<Straight> {
        const SPREAD: [f64; 6] = [0.137, 0.382, 0.5, 0.731, 0.919, 1.0];
        const REACH: [f64; 4] = [0.7, 3.1, 17.3, 101.9];
        let probes: Vec<f64> = match dom {
            ParameterDomain::Periodic { .. } => return None,
            ParameterDomain::Interval(lo, hi) if dom.is_finite() => {
                SPREAD.iter().map(|f| lo + (hi - lo) * f).collect()
            }
            ParameterDomain::Interval(lo, _) if lo.is_finite() => {
                REACH.iter().map(|r| lo + r).collect()
            }
            ParameterDomain::Interval(_, hi) if hi.is_finite() => {
                REACH.iter().map(|r| hi - r).collect()
            }
            ParameterDomain::Interval(..) => REACH.iter().flat_map(|&r| [-r, r]).collect(),
        };
        let u0 = dom.reference();
        let (origin, dir) = curve.first_derivatives(u0);
        let line = Straight { u0, origin, dir };
        if line.dir.magnitude() <= 0.0 {
            return None;
        }
        let size = |p: Point| (p - Point::ORIGIN).magnitude();
        probes
            .into_iter()
            .all(|u| {
                let (p, q) = (curve.point_at(u), line.point(u));
                p.distance(q) <= 1e-12 * (size(p) + size(line.origin) + p.distance(line.origin))
            })
            .then_some(line)
    }

    fn point(&self, u: f64) -> Point {
        self.origin + self.dir * (u - self.u0)
    }

    /// Parameter of the foot of `p` on the (unbounded) line.
    fn foot(&self, p: Point) -> f64 {
        self.u0 + (p - self.origin).dot(self.dir) / self.dir.square_magnitude()
    }
}

/// The extremum record of the pair `(s, t)`.
fn curve_pair_extremum<A, B>(a: &A, b: &B, s: f64, t: f64) -> CurveCurveExtremum
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    CurveCurveExtremum::new(s, a.point_at(s), t, b.point_at(t))
}

/// Generic closest-point enumeration between two curves of one tier,
/// two 3D curves or two 2D curves: every local minimum of the distance
/// `|C1(s) - C2(t)|` over the two parameter domains, ordered by distance.
///
/// Any 3D [`ParametricCurve`] is taken, a curve of the caller's own
/// included. In the plane, the 2D curves of the crate are taken
/// ([`Curve2D`](crate::Curve2D), each concrete 2D curve and
/// [`TrimmedCurve2D`](crate::TrimmedCurve2D) over them); the search runs
/// on them embedded in space at `z = 0`, with the same parameters, so
/// everything below holds in the plane as in space.
///
/// Bounded (trimmed) domains are honoured exactly: a minimum on a domain
/// end or a corner of the parameter box is reported, and only true local
/// minima of the distance restricted to the box are. The squared
/// distance is minimised over the box by Newton steps with exact second
/// derivatives (parameters held at a bound by an outward gradient stay
/// there), and each result must pass the box optimality check: every
/// parameter stationary (tangential residual within `tol.confusion`)
/// with a positive semidefinite Hessian, or held at a bound.
///
/// Two straight curves (lines, segments) are solved in closed form (one
/// minimum). A full line against any curve reduces exactly to the minima
/// of that curve's distance to the line, along the curve alone. Otherwise
/// the minimiser starts from the discrete local minima of a parameter
/// grid, and of dense samples along each finite edge of the box, and from
/// the points where the curves may meet: the pairs of grid spans whose
/// chords come within the curves' estimated offset from them are halved
/// until each holds a single minimum, so that crossings closer together
/// than the grid spacing are each seeded. An
/// unbounded parameter is seeded geometrically from its reference
/// parameter, out to many times the configuration's own length scale
/// (the gap between the curves plus their sizes), and the minimiser may
/// leave the seeded range.
///
/// Configurations with a continuous family of equal minima (parallel
/// lines, concentric circles, overlapping curves, a curve against a piece
/// of itself) report the distance once, with a representative pair read
/// from the curves, the same wherever they are placed: the first curve at
/// its reference parameter (`0` on a line or circle, the start of a
/// bounded domain) with its closest point on the second, when that pair
/// belongs to the family, else the second curve at its reference
/// parameter with its closest point on the first (the start of a piece
/// the family runs along), else a member found by the search.
///
/// A B-spline or Bézier curve, or trimmed curve over one, against
/// another or against a segment (a line trimmed to a range), unless both
/// are straight, is
/// seeded without a grid, the segment converted to a degree-1 B-spline.
/// Every pair of their Bézier pieces is halved until the derivatives of
/// the squared distance provably keep a sign over it (no minimum inside),
/// or its Hessian is provably indefinite, or provably positive definite,
/// so that it holds at most one minimum, which Newton's method reaches;
/// the box's edges are searched the same way along one curve. When every
/// pair settles, the minima are
/// [`Completeness::Complete`](crate::Completeness::Complete). A pair that
/// never settles (a tangent contact or a straight overlap, where the
/// pieces run along one line within `tol.confusion`), or a search that
/// runs out of its budget of halvings (a valley of equal distance, an
/// overlap along a curved stretch), leaves them
/// [`Completeness::Incomplete`](crate::Completeness::Incomplete), with the
/// reason and the regions where minima may be missing; the grid above
/// then seeds the minimiser too. Any other pair is
/// [`Completeness::Sampled`](crate::Completeness::Sampled): two straight
/// curves, two elementary curves, and a B-spline against a full line or a
/// curved elementary curve, where the grid is faster than the search on
/// its exact rational B-spline.
///
/// # Examples
///
/// ```
/// use geomcore::{Line, Point, Tolerance, Vector, curve_curve_extrema};
/// let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
/// let ext = curve_curve_extrema(&l1, &l2, Tolerance::DEFAULT);
/// assert_eq!(ext.len(), 1);
/// assert_eq!(ext[0].distance, 1.0);
///
/// // Parallel lines: one representative pair at the distance.
/// let l3 = Line::new(Point::new(0.0, 2.0, 0.0), Vector::X).unwrap();
/// let ext = curve_curve_extrema(&l1, &l3, Tolerance::DEFAULT);
/// assert_eq!(ext.len(), 1);
/// assert_eq!((ext[0].first_parameter, ext[0].distance), (0.0, 2.0));
///
/// // Two 2D curves: a circle above a line.
/// use geomcore::{Circle2D, Line2D, Point2D, Vector2D};
/// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// let circle = Circle2D::new(Point2D::new(1.0, 3.0), 1.0).unwrap();
/// let ext = curve_curve_extrema(&line, &circle, Tolerance::DEFAULT);
/// assert_eq!(ext.len(), 1);
/// assert!((ext[0].distance - 2.0).abs() < 1e-12);
/// assert!(ext[0].first_point.distance(Point2D::new(1.0, 0.0)) < 1e-9);
/// ```
pub fn curve_curve_extrema<A, B>(a: &A, b: &B, tol: Tolerance) -> GenericCurveCurveExtrema<A::Point>
where
    A: CurveExtremaOperand,
    B: CurveExtremaOperand<Point = A::Point>,
{
    A::extrema_from_space(space_curve_curve_extrema(&a.in_space(), &b.in_space(), tol))
}

/// [`curve_curve_extrema`] on two 3D curves: the search itself.
pub(crate) fn space_curve_curve_extrema<A, B>(a: &A, b: &B, tol: Tolerance) -> CurveCurveExtrema
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let (da, db) = (curve_domain(a), curve_domain(b));
    let straight = (Straight::of(a, da), Straight::of(b, db));
    if !matches!(straight, (Some(_), Some(_)))
        && let Some((mut seeds, completeness)) = curve_curve_seeds(a, b, tol)
    {
        if !completeness.is_complete() {
            seeds.extend(box_seeds(a, da, b, db));
        }
        let minima = refine_minima(a, da, b, db, seeds, tol);
        return CurveCurveExtrema::new(minima, completeness);
    }
    sampled_minima(a, da, b, db, straight, tol)
}

/// [`curve_curve_extrema`] through the grid and the closed forms alone,
/// never the piece-pair search: the sampled fallback of the intersection,
/// and the reference the searches are tested against.
pub(crate) fn sampled_curve_curve_extrema<A, B>(a: &A, b: &B, tol: Tolerance) -> CurveCurveExtrema
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let (da, db) = (curve_domain(a), curve_domain(b));
    let straight = (Straight::of(a, da), Straight::of(b, db));
    sampled_minima(a, da, b, db, straight, tol)
}

/// The minima on the grid and the closed forms, for the domains `da`
/// and `db` and the straight curves among the two.
fn sampled_minima<A, B>(
    a: &A,
    da: ParameterDomain,
    b: &B,
    db: ParameterDomain,
    straight: (Option<Straight>, Option<Straight>),
    tol: Tolerance,
) -> CurveCurveExtrema
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let minima = match straight {
        (Some(la), Some(lb)) => vec![straight_pair_minimum(a, la, da, b, lb, db, tol)],
        (Some(la), None) if da.is_unbounded() => line_curve_minima(a, la, b, db, tol),
        (None, Some(lb)) if db.is_unbounded() => line_curve_minima(b, lb, a, da, tol)
            .into_iter()
            .map(CurveCurveExtremum::swapped)
            .collect(),
        _ => refine_minima(a, da, b, db, box_seeds(a, da, b, db), tol),
    };
    CurveCurveExtrema::sampled(minima)
}

/// The minimum distance between two straight curves. The squared
/// distance is a convex quadratic, so its minimum over the parameter box
/// is the unconstrained (skew) solution when that lies inside, and
/// otherwise lies on an edge of the box, where it is the clamped foot of
/// the edge's end point. A parallel pair has a valley of minima; its
/// representative is the first edge candidate (the first curve's start)
/// or, for two full lines, the first curve at `0` and its foot.
fn straight_pair_minimum<A, B>(
    a: &A,
    la: Straight,
    da: ParameterDomain,
    b: &B,
    lb: Straight,
    db: ParameterDomain,
    tol: Tolerance,
) -> CurveCurveExtremum
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let (d1, d2) = (la.dir, lb.dir);
    let w = la.origin - lb.origin;
    let (aa, ab, bb) = (d1.dot(d1), d1.dot(d2), d2.dot(d2));
    // `aa * bb - ab²` without its cancellation near parallel;
    // `denom / (aa * bb)` is the squared sine of the angle between them.
    let denom = d1.cross(d2).square_magnitude();
    if denom > tol.angular * tol.angular * aa * bb {
        let (e1, e2) = (d1.dot(w), d2.dot(w));
        let s = la.u0 + (ab * e2 - bb * e1) / denom;
        let t = lb.u0 + (aa * e2 - ab * e1) / denom;
        if da.wrap(s) == s && db.wrap(t) == t {
            return curve_pair_extremum(a, b, s, t);
        }
    }
    let edges = da
        .finite_ends()
        .map(|s| (s, db.wrap(lb.foot(la.point(s)))))
        .chain(db.finite_ends().map(|t| (da.wrap(la.foot(lb.point(t))), t)));
    let mut best: Option<CurveCurveExtremum> = None;
    for (s, t) in edges {
        let e = curve_pair_extremum(a, b, s, t);
        // Strict: on a parallel valley the first edge stays.
        if best.is_none_or(|m| e.distance < m.distance) {
            best = Some(e);
        }
    }
    best.unwrap_or_else(|| {
        let s = da.reference();
        curve_pair_extremum(a, b, s, lb.foot(la.point(s)))
    })
}

/// Local minima of the distance between a full line and a curve, as
/// `(line parameter, curve parameter)` extrema. The line's closest point
/// to each curve point is its foot, so the minima are exactly those of
/// the curve's distance to the line along the curve parameter alone:
/// seeded at the discrete minima of dense samples and refined by bounded
/// Newton with exact second derivatives.
fn line_curve_minima<L, C>(
    line_curve: &L,
    line: Straight,
    curve: &C,
    dom: ParameterDomain,
    tol: Tolerance,
) -> Vec<CurveCurveExtremum>
where
    L: ParametricCurve,
    C: ParametricCurve,
{
    let unit = Direction::from_unit_vector(line.dir * (1.0 / line.dir.magnitude()));
    let perp = |v: Vector| v.component_perpendicular_to(unit);
    // Half the squared distance to the line, its first two derivatives
    // along the curve, and the curve's speed.
    let local = |t: f64| {
        let (p, c1, c2) = curve.second_derivatives(t);
        let r = perp(p - line.origin);
        let c1p = perp(c1);
        (
            (0.5 * r.dot(r), r.dot(c1), c1p.dot(c1p) + r.dot(c2)),
            c1.magnitude(),
        )
    };
    let extremum = |t: f64| curve_pair_extremum(line_curve, curve, line.foot(curve.point_at(t)), t);
    let gap = perp(curve.point_at(dom.reference()) - line.origin).magnitude();
    let scale = configuration_scale(gap, 0.0, curve_size(curve, dom));
    let seeds = dom.seeds(curve, scale, SEED_GRID_1D);
    let points: Vec<Point> = seeds.iter().map(|&t| curve.point_at(t)).collect();
    let values: Vec<f64> = points
        .iter()
        .map(|&p| perp(p - line.origin).magnitude())
        .collect();
    let rounding = distance_rounding(coordinate_size(points.iter().chain([&line.origin])));
    let bounds = dom.bounds();
    let mut found = Vec::new();
    let shape = [seeds.len(), 1];
    let same_points = |p: usize, q: usize| points[p].distance(points[q]) <= rounding;
    for (k, _) in grid_minima(
        &values,
        shape,
        [dom.is_periodic(), false],
        rounding,
        same_points,
    ) {
        let t = dom.wrap(newton_min_1d(|t| local(t).0, seeds[k], bounds, 1e-14, 100));
        let ((_, g, h), speed) = local(t);
        let optimal = is_box_minimum(
            [t, 0.0],
            [g, 0.0],
            [[h, 0.0], [0.0, 0.0]],
            [speed, 0.0],
            [bounds, None],
            tol,
        );
        if t.is_finite() && optimal {
            found.push(extremum(t));
        }
    }
    let distance = |t: f64| {
        line_curve
            .point_at(line.foot(curve.point_at(t)))
            .distance(curve.point_at(t))
    };
    let same_valley = |p: &CurveCurveExtremum, q: &CurveCurveExtremum| {
        let dt = dom.step(p.second_parameter, q.second_parameter);
        flat_between(
            p.distance,
            q.distance,
            |f| distance(p.second_parameter + f * dt),
            tol,
        )
    };
    let representative = |e: &CurveCurveExtremum| {
        let r = extremum(dom.reference());
        same_valley(e, &r).then_some(r)
    };
    merge_valleys(found, same_valley, representative)
}

/// Seeds for the local minima of the distance between two curves over
/// their parameter box, for [`refine_minima`]: the discrete local minima
/// of a parameter grid and of dense samples along each finite edge, and
/// the points where the curves may meet (see [`meeting_seeds`]).
fn box_seeds<A, B>(a: &A, da: ParameterDomain, b: &B, db: ParameterDomain) -> Vec<[f64; 2]>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let gap = a
        .point_at(da.reference())
        .distance(b.point_at(db.reference()));
    let scale = configuration_scale(gap, curve_size(a, da), curve_size(b, db));
    let (ss, ts) = (da.seeds(a, scale, SEED_GRID), db.seeds(b, scale, SEED_GRID));
    let pa: Vec<Point> = ss.iter().map(|&s| a.point_at(s)).collect();
    let pb: Vec<Point> = ts.iter().map(|&t| b.point_at(t)).collect();
    let values: Vec<f64> = pa
        .iter()
        .flat_map(|p| pb.iter().map(move |q| p.distance(*q)))
        .collect();
    let shape = [ss.len(), ts.len()];
    let rounding = distance_rounding(coordinate_size(pa.iter().chain(&pb)));
    let mut seeds: Vec<[f64; 2]> = grid_minima(
        &values,
        shape,
        [da.is_periodic(), db.is_periodic()],
        rounding,
        |p, q| {
            let m = shape[1];
            pa[p / m].distance(pa[q / m]) <= rounding && pb[p % m].distance(pb[q % m]) <= rounding
        },
    )
    .into_iter()
    .map(|(i, j)| [ss[i], ts[j]])
    .collect();
    // The grid sees one minimum where two lie closer together than its
    // spacing (two crossings of a thin lens, a crossing next to a near
    // touch), and the minimiser reaches only one of them.
    seeds.extend(meeting_seeds(
        a,
        &Span::chain(a, da, &ss, &pa),
        b,
        &Span::chain(b, db, &ts, &pb),
    ));
    // A minimum on an edge of the box can have a basin narrower than the
    // grid: seed each finite edge at the discrete minima of dense samples.
    let edge = |p: Point, curve: &dyn Fn(f64) -> Point, us: &[f64], dom: ParameterDomain| {
        let points: Vec<Point> = us.iter().map(|&u| curve(u)).collect();
        let values: Vec<f64> = points.iter().map(|q| p.distance(*q)).collect();
        let rounding = distance_rounding(coordinate_size(points.iter().chain([&p])));
        let same_points = |p: usize, q: usize| points[p].distance(points[q]) <= rounding;
        grid_minima(
            &values,
            [us.len(), 1],
            [dom.is_periodic(), false],
            rounding,
            same_points,
        )
        .into_iter()
        .map(|(k, _)| us[k])
        .collect::<Vec<f64>>()
    };
    let (dense_s, dense_t) = (
        da.seeds(a, scale, SEED_GRID_1D),
        db.seeds(b, scale, SEED_GRID_1D),
    );
    for s in da.finite_ends() {
        let on_b = edge(a.point_at(s), &|t| b.point_at(t), &dense_t, db);
        seeds.extend(on_b.into_iter().map(|t| [s, t]));
    }
    for t in db.finite_ends() {
        let on_a = edge(b.point_at(t), &|s| a.point_at(s), &dense_s, da);
        seeds.extend(on_a.into_iter().map(|s| [s, t]));
    }
    // A span's closest point may sit on the end of a period.
    seeds
        .into_iter()
        .map(|[s, t]| [da.wrap(s), db.wrap(t)])
        .collect()
}

/// A stretch of a curve between two parameters: its chord, and an
/// estimate of how far the curve strays from the chord.
#[derive(Debug, Clone, Copy)]
struct Span {
    lo: f64,
    hi: f64,
    start: Point,
    end: Point,
    middle: Point,
    /// Twice the larger of the curve's offset from the chord at the
    /// middle parameter and the offset `|C''| h^2 / 8` its second
    /// derivative there predicts: exact for constant curvature, with a
    /// margin for curvature that varies along the span.
    deviation: f64,
    /// The chord's midpoint, and the radius about it of a ball that holds
    /// the curve over the span: half the chord plus the deviation.
    centre: Point,
    reach: f64,
}

impl Span {
    fn new<C: ParametricCurve>(curve: &C, lo: f64, hi: f64, start: Point, end: Point) -> Span {
        let step = hi - lo;
        let (middle, _, second) = curve.second_derivatives(lo + 0.5 * step);
        let centre = start + (end - start) * 0.5;
        let deviation = 2.0
            * middle
                .distance(centre)
                .max(second.magnitude() * step * step / 8.0);
        Span {
            lo,
            hi,
            start,
            end,
            middle,
            deviation,
            centre,
            reach: 0.5 * start.distance(end) + deviation,
        }
    }

    /// The spans between consecutive seed parameters of a curve; a
    /// periodic curve's last span closes the period.
    fn chain<C: ParametricCurve>(
        curve: &C,
        dom: ParameterDomain,
        parameters: &[f64],
        points: &[Point],
    ) -> Vec<Span> {
        let closing = match dom {
            ParameterDomain::Periodic { start, period } => Some((start + period, points[0])),
            ParameterDomain::Interval(..) => None,
        };
        let nodes: Vec<(f64, Point)> = parameters
            .iter()
            .copied()
            .zip(points.iter().copied())
            .chain(closing)
            .collect();
        nodes
            .windows(2)
            .map(|pair| Span::new(curve, pair[0].0, pair[1].0, pair[0].1, pair[1].1))
            .collect()
    }

    /// The two halves of the span.
    fn halves<C: ParametricCurve>(&self, curve: &C) -> [Span; 2] {
        let middle = 0.5 * (self.lo + self.hi);
        [
            Span::new(curve, self.lo, middle, self.start, self.middle),
            Span::new(curve, middle, self.hi, self.middle, self.end),
        ]
    }

    fn chord(&self) -> Vector {
        self.end - self.start
    }

    /// The closest points of this span's chord and `other`'s, as
    /// parameters interpolated linearly along the spans, and their
    /// distance.
    fn closest(&self, other: &Span) -> ([f64; 2], f64) {
        let (along_p, along_q) = (self.chord(), other.chord());
        let [f, g] = closest_on_segments(self.start, along_p, other.start, along_q);
        (
            [
                self.lo + (self.hi - self.lo) * f,
                other.lo + (other.hi - other.lo) * g,
            ],
            (self.start + along_p * f).distance(other.start + along_q * g),
        )
    }

    /// Whether the curves over the two spans may meet: their chords come
    /// within the sum of the spans' deviations.
    fn may_meet(&self, other: &Span, chord_distance: f64) -> bool {
        chord_distance <= self.deviation + other.deviation
    }

    /// A cheap necessary condition for [`Span::may_meet`]: the balls
    /// about the chords' midpoints holding the curves over the spans
    /// overlap.
    fn within_reach(&self, other: &Span) -> bool {
        let reach = self.reach + other.reach;
        (self.centre - other.centre).square_magnitude() <= reach * reach
    }

    /// Whether the distance between the curves over the two spans has at
    /// most one local minimum: the chords cross at an angle wider than
    /// the curves' tangents can turn within the spans (about
    /// `8 deviation / chord length` each), so the curves cross at most
    /// once and their tangents are never parallel.
    fn single_minimum(&self, other: &Span) -> bool {
        let (along_p, along_q) = (self.chord(), other.chord());
        let (length_p, length_q) = (along_p.magnitude(), along_q.magnitude());
        let sine = along_p.cross(along_q).magnitude() / (length_p * length_q);
        sine > 8.0 * (self.deviation / length_p + other.deviation / length_q)
    }
}

/// Seeds near every point where the curves may meet: the pairs of spans
/// whose chords come within their deviations are halved until the curves
/// over them part ([`Span::may_meet`] fails) or hold a single minimum
/// ([`Span::single_minimum`]), which is seeded at the chords' closest
/// points. A uniform grid resolves two minima only when they lie a few
/// spans apart; two crossings of a thin lens, or a crossing next to a near
/// touch, can lie closer, at the scale where the curves stray from their
/// chords by the gap between them, and halving reaches that scale. Pairs
/// still unsettled after [`MAX_MEETING_HALVINGS`] halvings in all (a
/// tangency, an overlap) are seeded as they stand, the closest first, at
/// most [`MAX_SEEDS`].
fn meeting_seeds<A, B>(a: &A, spans_a: &[Span], b: &B, spans_b: &[Span]) -> Vec<[f64; 2]>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let mut pending: Vec<(Span, Span)> = spans_a
        .iter()
        .flat_map(|p| spans_b.iter().map(move |q| (*p, *q)))
        .filter(|(p, q)| p.within_reach(q))
        .collect();
    let mut settled: Vec<([f64; 2], f64)> = Vec::new();
    let mut halvings = 0;
    while let Some((p, q)) = pending.pop() {
        let (closest, distance) = p.closest(&q);
        if !p.may_meet(&q, distance) {
            continue;
        }
        if halvings == MAX_MEETING_HALVINGS || p.single_minimum(&q) {
            settled.push((closest, distance));
            continue;
        }
        halvings += 1;
        let (halves_p, halves_q) = (p.halves(a), q.halves(b));
        for half_p in halves_p {
            pending.extend(halves_q.iter().map(|half_q| (half_p, *half_q)));
        }
    }
    settled.sort_by(|x, y| x.1.total_cmp(&y.1));
    settled.truncate(MAX_SEEDS);
    settled.into_iter().map(|(seed, _)| seed).collect()
}

/// The fractions `[f, g]` in `[0, 1]` of the closest points
/// `p + f * along_p` and `q + g * along_q` of two segments; a degenerate
/// segment is its start point, and parallel segments take the first
/// segment's start.
fn closest_on_segments(p: Point, along_p: Vector, q: Point, along_q: Vector) -> [f64; 2] {
    let r = p - q;
    let (aa, bb) = (along_p.dot(along_p), along_q.dot(along_q));
    let (ar, br) = (along_p.dot(r), along_q.dot(r));
    if aa <= 0.0 {
        let g = if bb > 0.0 {
            (br / bb).clamp(0.0, 1.0)
        } else {
            0.0
        };
        return [0.0, g];
    }
    if bb <= 0.0 {
        return [(-ar / aa).clamp(0.0, 1.0), 0.0];
    }
    let ab = along_p.dot(along_q);
    let denominator = along_p.cross(along_q).square_magnitude();
    let f = if denominator > 0.0 {
        ((ab * br - ar * bb) / denominator).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let g = (ab * f + br) / bb;
    if g < 0.0 {
        [(-ar / aa).clamp(0.0, 1.0), 0.0]
    } else if g > 1.0 {
        [((ab - ar) / aa).clamp(0.0, 1.0), 1.0]
    } else {
        [f, g]
    }
}

/// The local minima reached by Newton minimisation of the squared
/// distance over the box from `seeds`, kept when box-optimal, one per
/// valley of equal distance.
pub(crate) fn refine_minima<A, B>(
    a: &A,
    da: ParameterDomain,
    b: &B,
    db: ParameterDomain,
    seeds: Vec<[f64; 2]>,
    tol: Tolerance,
) -> Vec<CurveCurveExtremum>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    // Half the squared distance with its gradient and Hessian, and the
    // two curves' speeds.
    let local = |s: f64, t: f64| {
        let (pa, a1, a2) = a.second_derivatives(s);
        let (pb, b1, b2) = b.second_derivatives(t);
        let r = pa - pb;
        let cross = -a1.dot(b1);
        (
            (
                0.5 * r.dot(r),
                [r.dot(a1), -r.dot(b1)],
                [
                    [a1.dot(a1) + r.dot(a2), cross],
                    [cross, b1.dot(b1) - r.dot(b2)],
                ],
            ),
            [a1.magnitude(), b1.magnitude()],
        )
    };
    let bounds = [da.bounds(), db.bounds()];
    let mut found = Vec::new();
    for seed in seeds {
        let [s, t] = newton_min_2d(|s, t| local(s, t).0, seed, bounds, 1e-14, 100);
        let (s, t) = (da.wrap(s), db.wrap(t));
        let ((_, g, h), speed) = local(s, t);
        if s.is_finite() && t.is_finite() && is_box_minimum([s, t], g, h, speed, bounds, tol) {
            found.push(curve_pair_extremum(a, b, s, t));
        }
    }
    let same_valley =
        |p: &CurveCurveExtremum, q: &CurveCurveExtremum| one_valley(a, da, b, db, p, q, tol);
    let representative = |e: &CurveCurveExtremum| valley_representative(a, da, b, db, e, tol);
    merge_valleys(found, same_valley, representative)
}

/// Steps of the walk along a valley in [`one_valley`].
const VALLEY_WALK_STEPS: usize = 16;

/// Whether two minima lie on one valley of equal distance: their
/// distances agree within `tol.confusion`, and one curve, walked from one
/// minimum's parameter to the other's, keeps its closest point on the
/// other curve (each found from the last) at that distance, and arrives
/// on the other minimum's point. Either curve may lead: a valley along a
/// piece of a closed curve is walked along the piece, as the shorter way
/// round the closed curve may leave the piece; a valley of one curve over
/// a single point of the other (a circle about a line through its
/// centre) is walked along the circle. The walk follows the valley where
/// it bends, which the straight way in the parameters does not. Copies of
/// one minimum pass at once, and so do two minima joined by the straight
/// way in both parameters, the valley of two copies of one curve. Only
/// minima where the distance is flat in some direction (see
/// [`lies_flat`]) are walked.
fn one_valley<A, B>(
    a: &A,
    da: ParameterDomain,
    b: &B,
    db: ParameterDomain,
    p: &CurveCurveExtremum,
    q: &CurveCurveExtremum,
    tol: Tolerance,
) -> bool
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    if (p.distance - q.distance).abs() > tol.confusion {
        return false;
    }
    if !p.parameters_apart(q, tol) {
        return true;
    }
    let ds = da.step(p.first_parameter, q.first_parameter);
    let dt = db.step(p.second_parameter, q.second_parameter);
    let straight = |f: f64| {
        a.point_at(p.first_parameter + f * ds)
            .distance(b.point_at(p.second_parameter + f * dt))
    };
    let along_first = (p.first_parameter, ds, p.second_parameter, q.second_point);
    let along_second = (p.second_parameter, dt, p.first_parameter, q.first_point);
    if flat_between(p.distance, q.distance, straight, tol) {
        return true;
    }
    lies_flat(a, b, p, tol)
        && lies_flat(a, b, q, tol)
        && (walks_along_valley(a, b, db, along_first, p.distance, tol)
            || walks_along_valley(b, a, da, along_second, p.distance, tol))
}

/// The smallest eigenvalue of the Hessian of half the squared distance
/// at a point of a valley, read per unit length along each curve, is zero
/// to rounding; at an isolated minimum where the curves cross at the angle
/// `a` it is about `sin(a)^2`. Above this, a minimum is on no valley.
const VALLEY_FLATNESS: f64 = 1e-4;

/// Whether the distance is flat to second order in some direction at
/// the minimum `e`, read by [`flat_minimum`]. A minimum on a valley of
/// equal distance is; an isolated one, but for a tangent contact or a
/// shallow crossing, is not, and needs no walk nor representative.
fn lies_flat<A, B>(a: &A, b: &B, e: &CurveCurveExtremum, tol: Tolerance) -> bool
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let (pa, a1, a2) = a.second_derivatives(e.first_parameter);
    let (pb, b1, b2) = b.second_derivatives(e.second_parameter);
    let r = pa - pb;
    let mixed = -a1.dot(b1);
    flat_minimum(
        [r.dot(a1), -r.dot(b1)],
        [
            [a1.dot(a1) + r.dot(a2), mixed],
            [mixed, b1.dot(b1) - r.dot(b2)],
        ],
        [a1.magnitude(), b1.magnitude()],
        tol,
    )
}

/// Whether a box minimum of half the squared distance, with gradient `g`,
/// Hessian `h` and parameter speeds `speed`, is flat to second order
/// along some direction it is free to move in: the smallest eigenvalue of
/// the Hessian over its stationary parameters (those not held at a bound,
/// as [`is_box_minimum`] reads them), each read per unit length, is at
/// most [`VALLEY_FLATNESS`]. A minimum held in both parameters is a
/// corner of the box, on no valley. The test is dimensionless, so it does
/// not change with the size of the shapes.
pub(crate) fn flat_minimum(g: [f64; 2], h: [[f64; 2]; 2], speed: [f64; 2], tol: Tolerance) -> bool {
    if !(speed[0] > 0.0 && speed[1] > 0.0) {
        return true;
    }
    let free = [0, 1].map(|i| g[i].abs() <= tol.confusion * speed[i]);
    let first = h[0][0] / (speed[0] * speed[0]);
    let mixed = 0.5 * (h[0][1] + h[1][0]) / (speed[0] * speed[1]);
    let second = h[1][1] / (speed[1] * speed[1]);
    let smallest = match free {
        [true, true] => {
            0.5 * (first + second)
                - (0.25 * (first - second) * (first - second) + mixed * mixed).sqrt()
        }
        [true, false] => first,
        [false, true] => second,
        [false, false] => return false,
    };
    smallest.is_nan() || smallest <= VALLEY_FLATNESS
}

/// Whether `leader`, walked from `start` by `step`, keeps its closest
/// point on `follower` (from `warm`, then from the last one) within
/// `tol.confusion` of `distance`, and the last of them is `arrival`.
fn walks_along_valley<L, F>(
    leader: &L,
    follower: &F,
    follower_domain: ParameterDomain,
    (start, step, warm, arrival): (f64, f64, f64, Point),
    distance: f64,
    tol: Tolerance,
) -> bool
where
    L: ParametricCurve,
    F: ParametricCurve,
{
    let mut warm = warm;
    let mut last = None;
    for k in 1..=VALLEY_WALK_STEPS {
        let point = leader.point_at(start + step * k as f64 / VALLEY_WALK_STEPS as f64);
        let (u, foot) = closest_from(follower, follower_domain, point, warm);
        if (foot.distance(point) - distance).abs() > tol.confusion {
            return false;
        }
        warm = u;
        last = Some(foot);
    }
    last.is_some_and(|foot| foot.distance(arrival) <= tol.confusion)
}

/// The closest point of `curve` to `point` that Newton's method reaches
/// from the parameter `start`, with its parameter.
pub(crate) fn closest_from<C: ParametricCurve>(
    curve: &C,
    domain: ParameterDomain,
    point: Point,
    start: f64,
) -> (f64, Point) {
    let objective = |u: f64| {
        let (on_curve, first, second) = curve.second_derivatives(u);
        let r = point - on_curve;
        (
            0.5 * r.dot(r),
            -r.dot(first),
            first.dot(first) - r.dot(second),
        )
    };
    let u = domain.wrap(newton_min_1d(objective, start, domain.bounds(), 1e-14, 100));
    (u, curve.point_at(u))
}

/// The representative of the valley of equal distance through `e`, or
/// `None` when it lies on none (see [`lies_flat`]): the first curve at its reference
/// parameter with its closest point on the second, when that pair is on
/// the valley, else the second curve at its reference parameter with its
/// closest point on the first. Both are read from the curves alone, so the
/// pair does not depend on where the search entered the valley nor on the
/// placement (ADR 0051). One of them is on any valley whose ends are the
/// ends of the curves' domains (an overlap of two pieces, a whole curve
/// against a piece of itself, parallel segments, concentric arcs): the
/// valley starts where the later of the two pieces starts.
fn valley_representative<A, B>(
    a: &A,
    da: ParameterDomain,
    b: &B,
    db: ParameterDomain,
    e: &CurveCurveExtremum,
    tol: Tolerance,
) -> Option<CurveCurveExtremum>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    if !lies_flat(a, b, e, tol) {
        return None;
    }
    let s = da.reference();
    let (t, on_second) = closest_from(b, db, a.point_at(s), e.second_parameter);
    let from_first = CurveCurveExtremum::new(s, a.point_at(s), t, on_second);
    if one_valley(a, da, b, db, e, &from_first, tol) {
        return Some(from_first);
    }
    let t = db.reference();
    let (s, on_first) = closest_from(a, da, b.point_at(t), e.first_parameter);
    let from_second = CurveCurveExtremum::new(s, on_first, t, b.point_at(t));
    one_valley(a, da, b, db, e, &from_second, tol).then_some(from_second)
}

/// Length scale of a configuration: the gap between the curves'
/// reference points plus both curves' sizes (`1` if that degenerates).
pub(crate) fn configuration_scale(gap: f64, size_a: f64, size_b: f64) -> f64 {
    let scale = gap + size_a + size_b;
    if scale.is_finite() && scale > 0.0 {
        scale
    } else {
        1.0
    }
}

/// Size of a curve for [`configuration_scale`]: the reach of a finite
/// domain's grid points from its reference point, or for an unbounded
/// domain the radius of curvature at the reference (`0` when straight).
pub(crate) fn curve_size<C: ParametricCurve>(curve: &C, dom: ParameterDomain) -> f64 {
    let u0 = dom.reference();
    if dom.is_finite() {
        let p0 = curve.point_at(u0);
        dom.seeds(curve, 1.0, SEED_GRID)
            .into_iter()
            .map(|u| curve.point_at(u).distance(p0))
            .fold(0.0, f64::max)
    } else {
        let (_, d1, d2) = curve.second_derivatives(u0);
        let bend = d1.cross(d2).magnitude();
        if bend > 0.0 {
            d1.magnitude().powi(3) / bend
        } else {
            0.0
        }
    }
}

/// The largest coordinate, in absolute value, of `points`: the size about
/// the origin their rounding is measured against.
pub(crate) fn coordinate_size<'a>(points: impl IntoIterator<Item = &'a Point>) -> f64 {
    points
        .into_iter()
        .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
        .fold(0.0, f64::max)
}

/// How many times `EPSILON` times the coordinates a distance between two
/// sampled points is rounded within: each point carries the rounding of
/// its evaluation, and the difference adds both.
const DISTANCE_ROUNDING: f64 = 4.0;

/// The rounding of a distance between points whose coordinates are at
/// most `size` (see [`coordinate_size`]): about the operands, not their
/// own distance (ADR 0051, rule 2).
pub(crate) fn distance_rounding(size: f64) -> f64 {
    DISTANCE_ROUNDING * f64::EPSILON * size
}

/// Nodes of a row-major `shape[0] x shape[1]` grid of samples where the
/// values (distances, rounded within `rounding`) have a discrete local
/// minimum, one per plateau, lowest first, at most [`MAX_SEEDS`].
///
/// A node is a minimum when none of its 8 neighbours (wrapping across a
/// periodic direction) is lower by more than `rounding`. Minima that
/// neighbour each other within `rounding` are one plateau, and only its
/// lowest node is kept (the first in the grid on a tie). A valley of equal
/// distance along a grid direction ties every node along it with the next,
/// and rounding, so the placement of the shapes, decides which of them
/// come out lower: taken one by one they would fill the [`MAX_SEEDS`]
/// ahead of a farther minimum, and how many depended on the placement
/// (#448). The refinement from any node of the plateau reaches the valley.
///
/// Two nodes whose samples are the same points (`same_points`) never join
/// a plateau: their distances tie because the shapes do not move between
/// them, as along the pole row of a sphere, not because a valley runs
/// there. Each keeps its own seed, since the refinement from a node of a
/// degenerate row cannot turn towards a minimum beside it.
pub(crate) fn grid_minima(
    values: &[f64],
    shape: [usize; 2],
    wrap: [bool; 2],
    rounding: f64,
    same_points: impl Fn(usize, usize) -> bool,
) -> Vec<(usize, usize)> {
    let [n, m] = shape;
    let neighbour = |i: usize, di: isize, len: usize, wrap: bool| {
        let k = i as isize + di;
        if (0..len as isize).contains(&k) {
            Some(k as usize)
        } else if wrap && len > 1 {
            Some(k.rem_euclid(len as isize) as usize)
        } else {
            None
        }
    };
    // Whether some neighbour of node `(i, j)` passes `test`.
    let any_neighbour = |i: usize, j: usize, test: &mut dyn FnMut(usize) -> bool| {
        for di in -1..=1 {
            let Some(x) = neighbour(i, di, n, wrap[0]) else {
                continue;
            };
            for dj in -1..=1 {
                if let Some(y) = neighbour(j, dj, m, wrap[1])
                    && test(x * m + y)
                {
                    return true;
                }
            }
        }
        false
    };
    let mut minimum = vec![false; n * m];
    for i in 0..n {
        for j in 0..m {
            let v = values[i * m + j];
            minimum[i * m + j] =
                v.is_finite() && !any_neighbour(i, j, &mut |x| values[x] < v - rounding);
        }
    }
    let lower = |p: usize, q: usize| values[p] < values[q] || (values[p] == values[q] && p < q);
    let mut seen = vec![false; n * m];
    let mut plateau = Vec::new();
    let mut out = Vec::new();
    for k in 0..n * m {
        if !minimum[k] || seen[k] {
            continue;
        }
        seen[k] = true;
        plateau.push(k);
        let mut lowest = k;
        while let Some(p) = plateau.pop() {
            if lower(p, lowest) {
                lowest = p;
            }
            any_neighbour(p / m, p % m, &mut |x| {
                if minimum[x]
                    && !seen[x]
                    && (values[x] - values[p]).abs() <= rounding
                    && !same_points(x, p)
                {
                    seen[x] = true;
                    plateau.push(x);
                }
                false
            });
        }
        out.push(lowest);
    }
    out.sort_by(|&p, &q| values[p].total_cmp(&values[q]).then(p.cmp(&q)));
    out.truncate(MAX_SEEDS);
    out.into_iter().map(|k| (k / m, k % m)).collect()
}

/// Box optimality of a minimiser's result `x` for half the squared
/// distance: every parameter is either stationary (tangential residual
/// `|g_i| / speed_i` within `tol.confusion`) or held at a bound by an
/// outward gradient, and the Hessian over the stationary parameters is
/// positive semidefinite (rejecting saddles a seed happened to start on).
/// The Hessian is read per unit length along each parameter, as
/// [`flat_minimum`] reads it: an angle's curvature grows with the square
/// of the size and a length's does not, so a slack against the raw entries
/// let a saddle pass on small shapes that it failed on large ones
/// (ADR 0051, rule 4).
pub(crate) fn is_box_minimum(
    x: [f64; 2],
    g: [f64; 2],
    h: [[f64; 2]; 2],
    speed: [f64; 2],
    bounds: [Option<(f64, f64)>; 2],
    tol: Tolerance,
) -> bool {
    let mut stationary = [false; 2];
    for i in 0..2 {
        let held = match bounds[i] {
            Some((lo, hi)) => (x[i] <= lo && g[i] > 0.0) || (x[i] >= hi && g[i] < 0.0),
            None => false,
        };
        stationary[i] = g[i].abs() <= tol.confusion * speed[i];
        if !stationary[i] && !held {
            return false;
        }
    }
    // A parameter that does not move the point keeps its own units.
    let per_length = speed.map(|v| if v > 0.0 { 1.0 / v } else { 1.0 });
    let h: [[f64; 2]; 2] =
        std::array::from_fn(|i| std::array::from_fn(|j| h[i][j] * per_length[i] * per_length[j]));
    let slack = 1e-6 * (h[0][0].abs() + h[1][1].abs() + h[0][1].abs() + h[1][0].abs());
    match stationary {
        [true, true] => {
            let (p, q, r) = (h[0][0], 0.5 * (h[0][1] + h[1][0]), h[1][1]);
            0.5 * (p + r) - (0.25 * (p - r) * (p - r) + q * q).sqrt() >= -slack
        }
        [true, false] => h[0][0] >= -slack,
        [false, true] => h[1][1] >= -slack,
        [false, false] => true,
    }
}

/// Whether two minima at distances `d0` and `d1` lie on one valley: the
/// distances agree and so does the distance along the straight parameter
/// path between them (`along(f)`, `f` in `[0, 1]`), within
/// `tol.confusion`. Copies of one minimum pass trivially.
pub(crate) fn flat_between(d0: f64, d1: f64, along: impl Fn(f64) -> f64, tol: Tolerance) -> bool {
    (d0 - d1).abs() <= tol.confusion
        && (1..8).all(|k| (along(k as f64 / 8.0) - d0).abs() <= tol.confusion)
}

/// Collapses minima on one valley of equal distance (and copies of one
/// minimum) into a single entry, lowest distance first, each reported at
/// its `representative`: the valley's own pair, read from the shapes
/// rather than from the member the search reached, or `None` to keep the
/// minimum, which lies on no valley. A minimum the search reached alone on
/// its valley is replaced too, so that the pair does not depend on how
/// many of the valley's points the search met. A representative that
/// falls on the valley of another entry the merge kept apart is reported
/// once.
pub(crate) fn merge_valleys<E: Extremum>(
    mut found: Vec<E>,
    same_valley: impl Fn(&E, &E) -> bool,
    representative: impl Fn(&E) -> Option<E>,
) -> Vec<E> {
    found.sort_by(|x, y| x.distance().total_cmp(&y.distance()));
    let mut kept: Vec<E> = Vec::new();
    for e in found {
        if !kept.iter().any(|k| same_valley(k, &e)) {
            kept.push(e);
        }
    }
    let (mut represented, mut moved) = (Vec::with_capacity(kept.len()), Vec::new());
    for e in &kept {
        match representative(e) {
            Some(r) => moved.push(r),
            None => represented.push(*e),
        }
    }
    for r in moved {
        if !represented.iter().any(|k| same_valley(k, &r)) {
            represented.push(r);
        }
    }
    represented.sort_by(|x, y| x.distance().total_cmp(&y.distance()));
    represented
}

#[cfg(test)]
mod tests;
