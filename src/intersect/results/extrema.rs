//! The local minima of the distance between two curves
//! ([`CurveCurveExtrema`]).

use crate::intersect::{Completeness, CurveCurveExtremum, CurveCurveRegion};
use std::ops::Deref;

/// The local minima of the distance between two curves, closest first,
/// and whether they are all of them (see
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
#[derive(Debug, Clone, PartialEq)]
pub struct CurveCurveExtrema {
    minima: Vec<CurveCurveExtremum>,
    completeness: Completeness<CurveCurveRegion>,
}

impl CurveCurveExtrema {
    /// `minima`, put in order of distance, with their completeness.
    pub(crate) fn new(
        mut minima: Vec<CurveCurveExtremum>,
        completeness: Completeness<CurveCurveRegion>,
    ) -> Self {
        minima.sort_by(|x, y| x.distance.total_cmp(&y.distance));
        CurveCurveExtrema {
            minima,
            completeness,
        }
    }

    /// Minima found from samples.
    pub(crate) fn sampled(minima: Vec<CurveCurveExtremum>) -> Self {
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
    pub fn minima(&self) -> &[CurveCurveExtremum] {
        &self.minima
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
    pub fn into_minima(self) -> Vec<CurveCurveExtremum> {
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

impl Deref for CurveCurveExtrema {
    type Target = [CurveCurveExtremum];

    fn deref(&self) -> &[CurveCurveExtremum] {
        &self.minima
    }
}

impl IntoIterator for CurveCurveExtrema {
    type Item = CurveCurveExtremum;
    type IntoIter = std::vec::IntoIter<CurveCurveExtremum>;

    fn into_iter(self) -> std::vec::IntoIter<CurveCurveExtremum> {
        self.minima.into_iter()
    }
}

impl<'a> IntoIterator for &'a CurveCurveExtrema {
    type Item = &'a CurveCurveExtremum;
    type IntoIter = std::slice::Iter<'a, CurveCurveExtremum>;

    fn into_iter(self) -> std::slice::Iter<'a, CurveCurveExtremum> {
        self.minima.iter()
    }
}
