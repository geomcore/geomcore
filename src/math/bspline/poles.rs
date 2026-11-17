//! The poles of a B-spline, read in place through the tier's pole algebra
//! of [`base::pole`](crate::base::pole).

use crate::Point;
use crate::base::pole::{HomogeneousPole, PolePoint};

/// The poles of a B-spline, read in place: Cartesian points and, if
/// rational, one weight per point. A surface's poles are row-major, rows
/// along `u`.
#[derive(Clone, Copy)]
pub(crate) struct Poles<'a, P: PolePoint = Point> {
    pub(crate) points: &'a [P],
    pub(crate) weights: Option<&'a [f64]>,
}

impl<P: PolePoint> Poles<'_, P> {
    /// Homogeneous coordinates of pole `i` (`w = 1` if nonrational).
    #[inline]
    pub(crate) fn homogeneous(&self, i: usize) -> P::Homogeneous {
        let p = self.points[i];
        match self.weights {
            Some(weights) => HomogeneousPole::fold(p, weights[i]),
            None => HomogeneousPole::fold(p, 1.0),
        }
    }
}
