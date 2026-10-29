//! Exact intersections: closed-form solutions between elementary shapes,
//! always complete.
//!
//! The solvers and their result types live by family, in `curve_curve`,
//! `curve_surface` and `surface_surface`; `conic` and `conic_2d` hold the
//! shared conic machinery. This module holds the trait and the reversed
//! argument orders.
//!
//! [`IntersectExact`] is the zero-cost typed path, implemented once per
//! exact pair, in both argument orders, with the pair's own result type as
//! its output. Most results live on the stack; the ones that hold curves
//! or a variable number of hits allocate:
//! [`LineTorusIntersection`](crate::LineTorusIntersection),
//! [`TorusConeIntersection::Rings`](crate::TorusConeIntersection::Rings),
//! [`CircleSurfaceIntersection`](crate::CircleSurfaceIntersection) and
//! [`ConicConicIntersection`](crate::ConicConicIntersection) hold a `Vec`. The router reads each typed
//! result straight into its family result (`routing/exact_family.rs`).

pub(crate) mod conic;
pub(crate) mod conic_2d;
pub(crate) mod curve_curve;
pub(crate) mod curve_surface;
pub(crate) mod surface_surface;

use crate::{
    Circle, Circle2D, Cone, Cylinder, Ellipse, Hyperbola, Line, Line2D, Parabola, Plane, Sphere,
    Tolerance, Torus,
};

/// Exact intersection of two shapes: the typed, stack-only path.
///
/// Each implementation returns the result type of its pair, so the
/// compiler knows which outcomes a pair can produce, and dispatch is
/// static. On success a result holds the intersecting geometry; otherwise
/// it classifies the pair as tangent, disjoint, parallel or coincident.
/// Tangency and coincidence are decided with [`Tolerance`], so
/// near-degenerate configurations classify robustly instead of collapsing
/// to noise. A configuration with no closed form reports "not analytic";
/// this trait never runs a numeric intersection. To get an answer for
/// every pair, with the method that produced it, use [`Intersect`](crate::Intersect).
pub trait IntersectExact<Rhs: ?Sized = Self> {
    /// The result of intersecting `Self` with `Rhs`.
    type Output;

    /// Intersects `self` with `other`, deciding tangency and coincidence
    /// with `tol`.
    fn intersect_exact(&self, other: &Rhs, tol: Tolerance) -> Self::Output;
}

/// Implements the reversed argument order of a pair whose result does
/// not depend on the order, by swapping the arguments.
macro_rules! reversed {
    ($($a:ty => $b:ty),* $(,)?) => {$(
        impl IntersectExact<$b> for $a {
            type Output = <$b as IntersectExact<$a>>::Output;

            /// The same result as in the other argument order.
            fn intersect_exact(&self, other: &$b, tol: Tolerance) -> Self::Output {
                other.intersect_exact(self, tol)
            }
        }
    )*};
}

reversed!(
    Sphere => Plane,
    Cylinder => Plane,
    Cone => Plane,
    Plane => Torus,
    Cylinder => Sphere,
    Cone => Sphere,
    Sphere => Torus,
    Cylinder => Cone,
    Cylinder => Torus,
    Cone => Torus,
);

reversed!(
    Plane => Line,
    Sphere => Line,
    Cylinder => Line,
    Cone => Line,
    Torus => Line,
    Plane => Circle,
    Sphere => Circle,
    Cylinder => Circle,
    Cone => Circle,
    Plane => Ellipse,
    Sphere => Ellipse,
    Cylinder => Ellipse,
    Cone => Ellipse,
    Plane => Parabola,
    Sphere => Parabola,
    Cylinder => Parabola,
    Cone => Parabola,
    Plane => Hyperbola,
    Sphere => Hyperbola,
    Cylinder => Hyperbola,
    Cone => Hyperbola,
);

reversed!(
    Circle => Line,
    Ellipse => Line,
    Parabola => Line,
    Hyperbola => Line,
    Circle2D => Line2D,
);

#[cfg(test)]
mod tests;
