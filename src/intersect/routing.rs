//! The routed intersection entry point.
//!
//! [`Intersect`] gives the full answer for every supported pair: the exact
//! intersection when a closed form exists, the numeric one otherwise, and
//! the result says which method ran. [`IntersectExact`] is the typed exact
//! path, and never falls back to numeric work.
//!
//! Which pairs have a closed form, and which solver answers each, is the
//! closed-form pair table in `closed_form`; `exact_family` reads each
//! pair's typed exact result into its family result. This module holds
//! the helpers the routed implementations share: the operand traits
//! ([`CurveOperand`], [`SurfaceOperand`], [`CurveOperand2D`]), which read
//! an operand as the basis the table matches on and give the parameter
//! range restricting a curve, and the parameter-range clipping that turns
//! an exact intersection of the basis into the intersection of the
//! trimmed curve. The operand lists of the [`Intersect`] implementations
//! come from the shape enums' variant lists, so a new shape type edits no
//! list here.

use crate::curves::PointExtrema;
use crate::curves::PointExtrema2D;
use crate::curves::{
    BSplineCurve, BSplineCurve2D, BezierCurve, BezierCurve2D, Curve, Curve2D, ElementaryCurve,
    ParametricCurve, ParametricCurve2D, TrimmedCurve, TrimmedCurve2D,
};
use crate::intersect::numeric::embedding_2d::{Embed, intersect_numeric_2d};
use crate::intersect::numeric::marching::marching_intersection;
use crate::intersect::results::{CurveCurveOverlapList, SurfaceSurfaceIntersection};
use crate::intersect::{
    CurveCurveIntersection, CurveCurveIntersection2D, CurveCurveOverlap, CurveSurfaceIntersection,
    GenericCurveCurveIntersection, InlineList, intersect_curve_curve, intersect_curve_surface,
};
use crate::surfaces::{
    BSplineSurface, BezierSurface, ElementarySurface, ExtrusionSurface, ParametricSurface,
    RevolutionSurface, Surface,
};
use crate::tier::operations::VectorOperations;
use crate::{Point, Point2D, TierPoint, Tolerance};

// The list macros come before the submodules, which use `shapes!` too.

/// Calls `$then! { $arguments [shapes] }` with the shapes of a shape
/// enum's variant list (ADR 0017), so the routing lists follow the
/// variants.
macro_rules! shapes {
    ([$then:ident $($arguments:tt)*] $enum:ty; $($variant:path: $shape:ty),+) => {
        $then! { $($arguments)* [$($shape),+] }
    };
}

/// Calls `$then! { $arguments [operands] }` with the operands of a
/// family: the shapes of a variant list, then the enums named after them.
/// A shape is an operand of every pair once it is a variant.
macro_rules! operands {
    (
        [$then:ident [$($arguments:tt)*] $($enums:ty),+]
        $enum:ty;
        $($variant:path: $shape:ty),+
    ) => {
        $then! { $($arguments)* [$($shape,)+ $($enums),+] }
    };
}

/// Calls `$then! { $arguments [curves] }` with every curve operand.
macro_rules! curve_operands {
    ($then:ident $($arguments:tt)*) => {
        crate::curves::curve_variants!(
            Curve,
            operands! { $then [$($arguments)*] Curve, ElementaryCurve }
        );
    };
}

/// Calls `$then! { $arguments [surfaces] }` with every surface operand.
macro_rules! surface_operands {
    ($then:ident $($arguments:tt)*) => {
        crate::surfaces::surface_variants!(
            Surface,
            operands! { $then [$($arguments)*] Surface, ElementarySurface }
        );
    };
}

/// Calls `$then! { $arguments [curves] }` with every 2D curve operand.
macro_rules! curve_operands_2d {
    ($then:ident $($arguments:tt)*) => {
        crate::curves::curve_2d_variants!(operands! { $then [$($arguments)*] Curve2D });
    };
}

mod closed_form;
mod exact_family;
mod swept_plane;
pub(crate) use closed_form::ExactFamily;
use closed_form::{CurveKind, CurveKind2D, SurfaceKind};
use swept_plane::Swept;

/// Routed intersection: exact when a closed form exists, numeric otherwise.
///
/// Unlike [`IntersectExact`](crate::IntersectExact), the typed exact path that never runs numeric
/// work, `intersect` works for every supported pair, and its result says which
/// method ran (`method()`) and whether it is complete.
///
/// The operands may come in either order; a curve-surface result always
/// reads the curve parameter first.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{IntersectionMethod, Line, Point, Sphere, Tolerance, Vector};
///
/// // A line along z through a sphere of radius 2 crosses it at both poles.
/// let line = Line::new(Point::new(0.0, 0.0, -5.0), Vector::Z).unwrap();
/// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
/// let result = line.intersect(&sphere, Tolerance::DEFAULT);
/// assert_eq!(result.method(), IntersectionMethod::Exact);
/// assert!(result.is_complete());
/// let parameters: Vec<Option<f64>> = result.hits().iter().map(|h| h.curve_parameter).collect();
/// assert_eq!(parameters, [Some(3.0), Some(7.0)]);
/// assert_eq!(sphere.intersect(&line, Tolerance::DEFAULT), result);
/// ```
pub trait Intersect<Rhs: ?Sized = Self> {
    /// The result of intersecting `Self` with `Rhs`.
    type Output;

    /// Intersects `self` with `other`, deciding tangency and coincidence
    /// with `tol`.
    fn intersect(&self, other: &Rhs, tol: Tolerance) -> Self::Output;
}

/// A curve in the tier of `P` as the readers of `exact_family` and the
/// overlap mapping take it: a concrete elementary curve, or the tier's
/// elementary enum ([`ElementaryCurve`] or [`Curve2D`]).
///
/// [`parameter_of_point`](Self::parameter_of_point) gives the parameter of
/// a hit's point when the curve is restricted to a range and the solver
/// did not compute it: the range test needs it. A coincidence is mapped on
/// the tier's enum, so the mapping is compiled once per tier.
pub(crate) trait CurveBasis<P: CurveTier>: Clone {
    /// The curve as the tier's elementary enum.
    fn coincident(&self) -> P::Elementary;

    /// The parameter bounds.
    fn range(&self) -> (f64, f64);

    /// The period, or `None` when the curve is not periodic.
    fn cycle_length(&self) -> Option<f64>;

    /// The point at `parameter`.
    fn point_at(&self, parameter: f64) -> P;

    /// The first derivative at `parameter`.
    fn first_derivative(&self, parameter: f64) -> P::Vector;

    /// The parameter of the curve's point closest to `point`.
    fn parameter_of_point(&self, point: P, tol: Tolerance) -> f64;

    /// The length of the first derivative at `parameter`: the line
    /// direction's length, the circle's radius.
    #[inline]
    fn speed(&self, parameter: f64) -> f64 {
        self.first_derivative(parameter).magnitude()
    }
}

/// A tier as the curve router reads it: its point, and the enum of its
/// elementary curves, on which a coincidence is mapped.
pub(crate) trait CurveTier: TierPoint {
    /// [`ElementaryCurve`] in space, [`Curve2D`] in the plane.
    type Elementary: CurveBasis<Self>;
}

impl CurveTier for Point {
    type Elementary = ElementaryCurve;
}

impl CurveTier for Point2D {
    type Elementary = Curve2D;
}

/// Makes every curve of a tier that converts into the tier's elementary
/// enum a [`CurveBasis`], from the tier's curve traits.
macro_rules! tier_curve_basis {
    ($point:ty, $extrema:ident, $parametric:ident, $enum:ty) => {
        impl<C: $extrema + Clone + Into<$enum>> CurveBasis<$point> for C {
            #[inline]
            fn coincident(&self) -> $enum {
                self.clone().into()
            }

            #[inline]
            fn range(&self) -> (f64, f64) {
                $parametric::bounds(self)
            }

            #[inline]
            fn cycle_length(&self) -> Option<f64> {
                $parametric::period(self)
            }

            #[inline]
            fn point_at(&self, parameter: f64) -> $point {
                $parametric::eval_point(self, parameter)
            }

            #[inline]
            fn first_derivative(&self, parameter: f64) -> <$point as TierPoint>::Vector {
                $parametric::eval_derivative(self, parameter, 1)
            }

            #[inline]
            fn parameter_of_point(&self, point: $point, tol: Tolerance) -> f64 {
                $extrema::parameter_of(self, point, tol)
            }
        }
    };
}

tier_curve_basis!(Point, PointExtrema, ParametricCurve, ElementaryCurve);
tier_curve_basis!(Point2D, PointExtrema2D, ParametricCurve2D, Curve2D);

/// How a curve enters the router: its basis, the shape the closed-form
/// pair table matches it on (see `closed_form`), and the parameter range
/// it is restricted to (`None` = the basis's own bounds).
pub(crate) trait CurveOperand: ParametricCurve {
    /// The basis type: the curve's own type when it is a concrete
    /// elementary curve, so the pair is solved by its typed solver;
    /// [`ElementaryCurve`], or a [`CurveKind`], when its kind is only known
    /// at run time.
    type Basis: Copy;

    /// The basis and the restricting parameter range.
    fn basis_and_range(&self) -> Operand<Self::Basis>;
}

/// Makes each concrete elementary curve, and [`ElementaryCurve`], its own
/// basis.
macro_rules! own_curve_basis {
    ([$($curve:ty),*]) => {$(
        impl CurveOperand for $curve {
            type Basis = $curve;

            #[inline]
            fn basis_and_range(&self) -> Operand<$curve> {
                (*self, None)
            }
        }
    )*};
}

crate::curves::curve_variants!(ElementaryCurve, shapes! { own_curve_basis });
own_curve_basis!([ElementaryCurve]);

impl CurveOperand for BSplineCurve {
    type Basis = CurveKind;

    #[inline]
    fn basis_and_range(&self) -> Operand<CurveKind> {
        (CurveKind::Unlisted, None)
    }
}

impl CurveOperand for BezierCurve {
    type Basis = CurveKind;

    #[inline]
    fn basis_and_range(&self) -> Operand<CurveKind> {
        (CurveKind::Unlisted, None)
    }
}

impl CurveOperand for Curve {
    type Basis = CurveKind;

    #[inline]
    fn basis_and_range(&self) -> Operand<CurveKind> {
        (CurveKind::from(self), None)
    }
}

impl<C: CurveOperand> CurveOperand for TrimmedCurve<C> {
    type Basis = C::Basis;

    /// A trimmed curve of a trimmed curve keeps the inner basis and the
    /// outer range: the outer range lies within the inner one.
    #[inline]
    fn basis_and_range(&self) -> Operand<C::Basis> {
        let (basis, _) = self.basis().basis_and_range();
        (basis, Some((self.first(), self.last())))
    }
}

/// How a surface enters the router: its basis, the shape the closed-form
/// pair table matches it on (see `closed_form`).
pub(crate) trait SurfaceOperand: ParametricSurface {
    /// The basis type: the surface's own type when it is a concrete
    /// elementary surface, so a curve meeting it is solved by the pair's
    /// typed solver; [`ElementarySurface`], or a [`SurfaceKind`], when its
    /// kind is only known at run time. A surface pair is matched on the
    /// kinds the two bases convert into.
    type Basis<'a>: Copy + Into<SurfaceKind<'a>>
    where
        Self: 'a;

    /// The basis.
    fn surface_basis(&self) -> Self::Basis<'_>;
}

/// Makes each concrete elementary surface, and [`ElementarySurface`], its
/// own basis.
macro_rules! own_surface_basis {
    ([$($surface:ty),*]) => {$(
        impl SurfaceOperand for $surface {
            type Basis<'a> = $surface;

            #[inline]
            fn surface_basis(&self) -> $surface {
                *self
            }
        }
    )*};
}

crate::surfaces::surface_variants!(ElementarySurface, shapes! { own_surface_basis });
own_surface_basis!([ElementarySurface]);

impl SurfaceOperand for ExtrusionSurface {
    type Basis<'a> = SurfaceKind<'a>;

    #[inline]
    fn surface_basis(&self) -> SurfaceKind<'_> {
        SurfaceKind::Swept(Swept::Extrusion(self))
    }
}

impl SurfaceOperand for RevolutionSurface {
    type Basis<'a> = SurfaceKind<'a>;

    #[inline]
    fn surface_basis(&self) -> SurfaceKind<'_> {
        SurfaceKind::Swept(Swept::Revolution(self))
    }
}

impl SurfaceOperand for BSplineSurface {
    type Basis<'a> = SurfaceKind<'a>;

    #[inline]
    fn surface_basis(&self) -> SurfaceKind<'_> {
        SurfaceKind::Unlisted
    }
}

impl SurfaceOperand for BezierSurface {
    type Basis<'a> = SurfaceKind<'a>;

    #[inline]
    fn surface_basis(&self) -> SurfaceKind<'_> {
        SurfaceKind::Unlisted
    }
}

impl SurfaceOperand for Surface {
    type Basis<'a> = SurfaceKind<'a>;

    /// The kind of the surface the variant holds, built in one match.
    #[allow(
        clippy::useless_conversion,
        reason = "a swept or unlisted variant's basis is its kind already"
    )]
    #[inline]
    fn surface_basis(&self) -> SurfaceKind<'_> {
        crate::surfaces::surface_dispatch!(self, surface => surface.surface_basis().into())
    }
}

/// The representative in `range` of a point at `parameter`, or `None`
/// when the point lies farther than `tol.confusion` in space beyond the
/// range.
///
/// On a periodic curve the parameter is first shifted by whole periods
/// into `[range.0, range.0 + period)`. A parameter inside the range is
/// returned as it is. One beyond an end is measured in space: the
/// parameter excess times `speed` at that end (the line direction's
/// length, the circle's radius, the derivative's length on other curves).
/// Within `tol.confusion` the point is kept, snapped to that end. On a
/// periodic curve the end measured is the nearer one modulo the period,
/// so a range starting or ending at the seam keeps a point a hair on the
/// other side of it. `speed` is only called for a point outside the range.
pub(crate) fn in_range(
    parameter: f64,
    range: (f64, f64),
    period: Option<f64>,
    speed: impl FnOnce(f64) -> f64,
    tol: Tolerance,
) -> Option<f64> {
    let (before_start, past_end) = match period {
        None => {
            if parameter >= range.0 && parameter <= range.1 {
                return Some(parameter);
            }
            (range.0 - parameter, parameter - range.1)
        }
        Some(period) => {
            let shifted = parameter - period * ((parameter - range.0) / period).floor();
            if shifted <= range.1 {
                return Some(shifted);
            }
            (range.0 + period - shifted, shifted - range.1)
        }
    };
    let (end, excess) = if past_end > 0.0 && (before_start <= 0.0 || past_end <= before_start) {
        (range.1, past_end)
    } else {
        (range.0, before_start)
    };
    (excess * speed(end) <= tol.confusion).then_some(end)
}

/// The end of `range` nearest `parameter`, measured modulo `period` on a
/// periodic curve: the representative of a point that touches the range
/// from outside it. Clamping would pick the wrong end for a point just past
/// the seam (angle `5e-8` against a range ending at 2π).
fn nearest_end(parameter: f64, range: (f64, f64), period: Option<f64>) -> f64 {
    let Some(period) = period else {
        return parameter.clamp(range.0, range.1);
    };
    let distance = |end: f64| {
        let past = (parameter - end).rem_euclid(period);
        past.min(period - past)
    };
    if distance(range.0) <= distance(range.1) {
        range.0
    } else {
        range.1
    }
}

/// The parameter a hit keeps on one curve, or `None` when the hit lies
/// outside the curve's range.
///
/// On a curve restricted to `range` the range test needs the parameter:
/// the `stored` one, or the one `recover` computes from the point when the
/// solver did not, and the hit keeps its representative in the range (see
/// [`in_range`]). On an unrestricted curve the hit keeps what the solver
/// stored, computed or not.
#[inline]
fn kept_parameter(
    stored: Option<f64>,
    range: Option<(f64, f64)>,
    period: Option<f64>,
    recover: impl FnOnce() -> f64,
    speed: impl FnOnce(f64) -> f64,
    tol: Tolerance,
) -> Option<Option<f64>> {
    match range {
        Some(range) => {
            in_range(stored.unwrap_or_else(recover), range, period, speed, tol).map(Some)
        }
        None => Some(stored),
    }
}

/// Clips an overlap `(a, b)` to `range`, shifting by periods like
/// [`in_range`].
///
/// Returns the pieces of the overlap that lie in the range, ordered by
/// start, each expressed in the range's representation: none, one, or two
/// when the overlap meets the range on both sides of the seam of a periodic
/// curve. An overlap covering the whole period returns the range as one
/// piece. The three slots hold the three shifted copies the search tries,
/// so the list never fills up.
///
/// Lengths are measured in space: `speed` gives the length of the curve's
/// derivative at a parameter (the line direction's length, the circle's
/// radius). A piece whose length in space is at most `tol.confusion`,
/// including a gap that small between the overlap and the range, is kept
/// as a single point (see [`is_single_point`]): the overlap only touches
/// the range there. Such a point is dropped when it is, up to a period,
/// the end of another piece, so a touch is never reported twice.
pub(crate) fn clip_range(
    overlap: (f64, f64),
    range: (f64, f64),
    period: Option<f64>,
    speed: &impl Fn(f64) -> f64,
    tol: Tolerance,
) -> InlineList<(f64, f64), 3> {
    let clip = |(low, high): (f64, f64)| {
        let clipped = (low.max(range.0), high.min(range.1));
        (clipped.1 >= clipped.0 || is_single_point(clipped, speed, tol)).then_some(clipped)
    };
    let mut pieces = InlineList::new();
    let Some(period) = period else {
        if let Some(piece) = clip(overlap) {
            pieces.push(piece);
        }
        return pieces;
    };
    let missing = period - (overlap.1 - overlap.0);
    if missing <= 0.0 || missing * speed(overlap.1) <= tol.confusion {
        pieces.push(range);
        return pieces;
    }
    let first_copy = ((range.0 - overlap.1) / period).floor();
    let mut copies = InlineList::<(f64, f64), 3>::new();
    for step in 0..3 {
        let shift = (first_copy + f64::from(step)) * period;
        if let Some(piece) = clip((overlap.0 + shift, overlap.1 + shift)) {
            copies.push(piece);
        }
    }
    let same_point = |a: f64, b: f64| {
        let gap = (a - b).rem_euclid(period);
        gap.min(period - gap) * speed(a) <= tol.confusion
    };
    for (index, &piece) in copies.iter().enumerate() {
        let repeated = is_single_point(piece, speed, tol)
            && copies.iter().enumerate().any(|(other_index, &other)| {
                let earlier_point = is_single_point(other, speed, tol) && other_index < index;
                let longer = !is_single_point(other, speed, tol);
                (earlier_point || longer)
                    && (same_point(piece.0, other.0) || same_point(piece.0, other.1))
            });
        if !repeated {
            pieces.push(piece);
        }
    }
    pieces.sort_by(|a, b| a.0.total_cmp(&b.0));
    pieces
}

/// Whether a piece clipped by [`clip_range`] is at most `tol.confusion`
/// long in space, the parameter length times the curve's `speed` there: a
/// single point where the overlap touches the range. A piece with its ends
/// swapped is a gap, measured the same way.
pub(crate) fn is_single_point(
    (low, high): (f64, f64),
    speed: &impl Fn(f64) -> f64,
    tol: Tolerance,
) -> bool {
    let length = high - low;
    length.is_finite() && length.abs() * speed(0.5 * (low + high)) <= tol.confusion
}

/// Intersects a curve with a surface: in closed form when the table lists
/// the pair, by the numeric search otherwise.
#[inline]
pub(crate) fn route_curve_surface<'s, C: CurveOperand, S: SurfaceOperand>(
    curve: &C,
    surface: &'s S,
    tol: Tolerance,
) -> CurveSurfaceIntersection
where
    C::Basis: ExactFamily<S::Basis<'s>, Family = CurveSurfaceIntersection>,
{
    // A trimmed curve is searched over its own range only: no clipping.
    let numeric_search = || intersect_curve_surface(curve, surface, tol);
    let (basis, range) = curve.basis_and_range();
    basis.exact_family(&surface.surface_basis(), (range, None), tol, numeric_search)
}

/// Implements [`Intersect`] for every curve of the first list with every
/// surface of the second, in both orders. The surface-first order returns
/// the same result as the curve-first one.
macro_rules! curve_surface {
    (@pairs $curve:ty; [$($surface:ty),*]) => {$(
        impl Intersect<$surface> for $curve {
            type Output = CurveSurfaceIntersection;

            #[inline]
            fn intersect(&self, other: &$surface, tol: Tolerance) -> CurveSurfaceIntersection {
                route_curve_surface(self, other, tol)
            }
        }

        impl Intersect<$curve> for $surface {
            type Output = CurveSurfaceIntersection;

            #[inline]
            fn intersect(&self, other: &$curve, tol: Tolerance) -> CurveSurfaceIntersection {
                route_curve_surface(other, self, tol)
            }
        }
    )*};
    ([$($curve:ty),*] $surfaces:tt) => {$(
        curve_surface!(@pairs $curve; $surfaces);
    )*};
}

curve_operands!(surface_operands curve_surface);

/// Implements [`Intersect`] between a trimmed curve over any curve operand
/// and each surface, in both orders.
macro_rules! trimmed_curve_surface {
    ([$($surface:ty),*]) => {$(
        impl<X: CurveOperand> Intersect<$surface> for TrimmedCurve<X>
        where
            for<'a> X::Basis:
                ExactFamily<<$surface as SurfaceOperand>::Basis<'a>, Family = CurveSurfaceIntersection>,
        {
            type Output = CurveSurfaceIntersection;

            #[inline]
            fn intersect(&self, other: &$surface, tol: Tolerance) -> CurveSurfaceIntersection {
                route_curve_surface(self, other, tol)
            }
        }

        impl<X: CurveOperand> Intersect<TrimmedCurve<X>> for $surface
        where
            for<'a> X::Basis:
                ExactFamily<<$surface as SurfaceOperand>::Basis<'a>, Family = CurveSurfaceIntersection>,
        {
            type Output = CurveSurfaceIntersection;

            #[inline]
            fn intersect(&self, other: &TrimmedCurve<X>, tol: Tolerance) -> CurveSurfaceIntersection {
                route_curve_surface(other, self, tol)
            }
        }
    )*};
}

surface_operands!(trimmed_curve_surface);

/// Intersects two surfaces: in closed form when the table lists the
/// pair, by marching otherwise. The pair is matched on the two kinds: the
/// match costs nothing next to a closed-form solver or a march.
#[inline]
pub(crate) fn route_surface_surface<A: SurfaceOperand, B: SurfaceOperand>(
    a: &A,
    b: &B,
    tol: Tolerance,
) -> SurfaceSurfaceIntersection {
    let numeric_search = || marching_intersection(a, b, tol);
    let first: SurfaceKind<'_> = a.surface_basis().into();
    let second: SurfaceKind<'_> = b.surface_basis().into();
    first.exact_family(&second, (None, None), tol, numeric_search)
}

/// Implements [`Intersect`] for every surface of the list with every
/// surface of the list, in both orders (the full product).
macro_rules! surface_surface {
    (@pairs $first:ty; [$($second:ty),*]) => {$(
        impl Intersect<$second> for $first {
            type Output = SurfaceSurfaceIntersection;

            #[inline]
            fn intersect(&self, other: &$second, tol: Tolerance) -> SurfaceSurfaceIntersection {
                route_surface_surface(self, other, tol)
            }
        }
    )*};
    ([$($first:ty),*] $seconds:tt) => {$(
        surface_surface!(@pairs $first; $seconds);
    )*};
}

surface_operands!(surface_operands surface_surface);

/// Intersects two curves of one tier, each read as its basis and range:
/// in closed form when the table lists the pair, by the tier's numeric
/// search otherwise. A trimmed curve is searched over its own range only:
/// no clipping.
///
/// The operands are read here, not by the caller: a basis read from an
/// enum is as large as the enum, and moving two of them in measured 5 %
/// slower on a line and a circle read through `Curve`.
#[inline]
fn route_curve_curve<
    P: CurveTier,
    A: ExactFamily<B, Family = GenericCurveCurveIntersection<P>>,
    B,
>(
    read_first: impl FnOnce() -> Operand<A>,
    read_second: impl FnOnce() -> Operand<B>,
    tol: Tolerance,
    numeric_search: impl FnOnce() -> GenericCurveCurveIntersection<P>,
) -> GenericCurveCurveIntersection<P> {
    let (first, first_range) = read_first();
    let (second, second_range) = read_second();
    first.exact_family(&second, (first_range, second_range), tol, numeric_search)
}

/// A basis and the parameter range restricting it.
type Operand<C> = (C, Option<(f64, f64)>);

/// The overlaps of two coincident curves restricted to `first_range` and
/// `second_range`.
///
/// Each overlap holds a monotone range in each curve's own range: the
/// second range is mapped onto the first curve and clipped there (two
/// pieces when it crosses the first curve's seam), then each piece is
/// mapped back onto the second curve and split where it crosses the
/// second curve's seam.
///
/// A piece at most `tol.confusion` long in space, where the ranges only
/// touch (two consecutive pieces of one line or circle), is not an overlap: it
/// comes back in [`CoincidentPieces::touches`] as a parameter pair, each
/// parameter the point's representative in its range.
fn coincident_overlaps<P: CurveTier, C: CurveBasis<P>>(
    (first, first_range): (&C, (f64, f64)),
    (second, second_range): (&C, (f64, f64)),
    tol: Tolerance,
) -> CoincidentPieces {
    let same_direction = runs_along(second, second_range, first, tol);
    let on_first = mapped_range(second, second_range, first, same_direction, tol);
    let mut overlaps = CurveCurveOverlapList::new();
    let mut touches = InlineList::new();
    let first_speed = |parameter: f64| first.speed(parameter);
    let second_speed = |parameter: f64| second.speed(parameter);
    let first_pieces = clip_range(
        on_first,
        first_range,
        first.cycle_length(),
        &first_speed,
        tol,
    );
    for &piece in first_pieces.iter() {
        if is_single_point(piece, &first_speed, tol) {
            let first_parameter = (0.5 * (piece.0 + piece.1)).clamp(first_range.0, first_range.1);
            let found = parameter_of_image(second, first, first_parameter, tol);
            let period = second.cycle_length();
            let second_parameter = in_range(found, second_range, period, |t| second.speed(t), tol)
                .unwrap_or_else(|| nearest_end(found, second_range, period));
            touches.push((first_parameter, second_parameter));
            continue;
        }
        let on_second = mapped_range(first, piece, second, same_direction, tol);
        let pieces_on_second = clip_range(
            on_second,
            second_range,
            second.cycle_length(),
            &second_speed,
            tol,
        );
        if let [single] = pieces_on_second[..] {
            overlaps.push(CurveCurveOverlap {
                first_range: piece,
                second_range: single,
                same_direction,
            });
            continue;
        }
        // The piece crosses the second curve's seam: read the first curve's
        // part of each side back from its end points.
        for &part in pieces_on_second.iter() {
            let end = |parameter: f64| {
                let found = parameter_of_image(first, second, parameter, tol);
                let period = first.cycle_length();
                in_range(found, piece, period, |t| first.speed(t), tol)
                    .unwrap_or_else(|| nearest_end(found, piece, period))
            };
            let (start, finish) = (end(part.0), end(part.1));
            overlaps.push(CurveCurveOverlap {
                first_range: (start.min(finish), start.max(finish)),
                second_range: part,
                same_direction,
            });
        }
    }
    overlaps.sort_by(|p, q| p.first_range.0.total_cmp(&q.first_range.0));
    CoincidentPieces { overlaps, touches }
}

/// What two coincident curves restricted to parameter ranges share.
struct CoincidentPieces {
    /// The overlaps, ordered by start on the first curve.
    overlaps: CurveCurveOverlapList,
    /// The `(first, second)` parameters of the points where the ranges
    /// only touch, ordered by first parameter. At most one per clipped
    /// piece, so three slots never fill up.
    touches: InlineList<(f64, f64), 3>,
}

/// Whether the coincident curves `from` and `to` run the same way: the
/// sign of the tangents' dot product at a sample of `range`.
fn runs_along<P: CurveTier, C: CurveBasis<P>>(
    from: &C,
    range: (f64, f64),
    to: &C,
    tol: Tolerance,
) -> bool {
    let sample = match (range.0.is_finite(), range.1.is_finite()) {
        (true, true) => 0.5 * (range.0 + range.1),
        (true, false) => range.0,
        (false, true) => range.1,
        (false, false) => 0.0,
    };
    let image = parameter_of_image(to, from, sample, tol);
    from.first_derivative(sample)
        .dot(to.first_derivative(image))
        >= 0.0
}

/// The parameter, on `curve`, of the point `other` has at `parameter`;
/// `other` coincides with `curve`.
#[inline]
fn parameter_of_image<P: CurveTier, C: CurveBasis<P>>(
    curve: &C,
    other: &C,
    parameter: f64,
    tol: Tolerance,
) -> f64 {
    curve.parameter_of_point(other.point_at(parameter), tol)
}

/// The image of `range` of curve `from` on the coincident curve `to`, as a
/// monotone interval that may lie outside `to`'s own representation.
/// `along` says whether the curves run the same way.
///
/// An infinite end maps to the infinity its direction gives on `to`. On a
/// periodic curve the image starts at the image of the range's start in
/// the direction of travel, and a full period maps to a full period.
fn mapped_range<P: CurveTier, C: CurveBasis<P>>(
    from: &C,
    range: (f64, f64),
    to: &C,
    along: bool,
    tol: Tolerance,
) -> (f64, f64) {
    let image = |parameter: f64| parameter_of_image(to, from, parameter, tol);
    let Some(period) = to.cycle_length() else {
        let end = |parameter: f64| {
            if parameter.is_finite() {
                image(parameter)
            } else if along {
                parameter
            } else {
                -parameter
            }
        };
        let (start, finish) = (end(range.0), end(range.1));
        return (start.min(finish), start.max(finish));
    };
    let full_period = from.cycle_length().unwrap_or(period);
    if range.1 - range.0 >= full_period - tol.parametric {
        let start = image(range.0);
        return (start, start + period);
    }
    let (start, finish) = if along {
        (image(range.0), image(range.1))
    } else {
        (image(range.1), image(range.0))
    };
    (start, start + (finish - start).rem_euclid(period))
}

/// Implements [`Intersect`] between the curves of one tier, described as
/// `(operand trait, trimmed curve type, numeric search, result type)`:
/// every curve of the list with every curve of the list, the receiver read
/// first; a trimmed curve over any operand with each curve of the list, in
/// both orders; and two trimmed curves.
macro_rules! curve_curve {
    (@row $tier:tt $first:ty; [$($second:ty),*]) => {$(
        curve_curve!(@pair $tier [] $first, $second;);
    )*};
    (
        @pair ($operand:ident, $trimmed:ident, $search:ident, $output:ty)
        [$($generics:tt)*] $first:ty, $second:ty; $($bounds:tt)*
    ) => {
        impl<$($generics)*> Intersect<$second> for $first
        where
            $($bounds)*
        {
            type Output = $output;

            #[inline]
            fn intersect(&self, other: &$second, tol: Tolerance) -> $output {
                let numeric_search = || $search(self, other, tol);
                let (read_first, read_second) = (|| self.basis_and_range(), || other.basis_and_range());
                route_curve_curve(read_first, read_second, tol, numeric_search)
            }
        }
    };
    (($operand:ident, $trimmed:ident, $search:ident, $output:ty) [$($curve:ty),*] $seconds:tt) => {
        $(curve_curve!(@row ($operand, $trimmed, $search, $output) $curve; $seconds);)*
        $(
            curve_curve!(
                @pair ($operand, $trimmed, $search, $output) [X: $operand] $trimmed<X>, $curve;
                X::Basis: ExactFamily<<$curve as $operand>::Basis, Family = $output>,
            );
            curve_curve!(
                @pair ($operand, $trimmed, $search, $output) [X: $operand] $curve, $trimmed<X>;
                <$curve as $operand>::Basis: ExactFamily<X::Basis, Family = $output>,
            );
        )*
        curve_curve!(
            @pair ($operand, $trimmed, $search, $output)
            [X: $operand, Y: $operand] $trimmed<X>, $trimmed<Y>;
            X::Basis: ExactFamily<Y::Basis, Family = $output>,
        );
    };
}

curve_operands!(curve_operands curve_curve (
    CurveOperand, TrimmedCurve, intersect_curve_curve, CurveCurveIntersection
));

/// How a 2D curve enters the router: its basis, the shape the closed-form
/// pair table matches it on (see `closed_form`), and the parameter range
/// it is restricted to.
pub(crate) trait CurveOperand2D: Embed {
    /// The basis type: the curve's own type when it is a concrete
    /// elementary curve, so the pair is solved by its typed solver; a
    /// [`CurveKind2D`] when its kind is only known at run time.
    type Basis: Clone;

    /// The basis and the restricting parameter range.
    fn basis_and_range(&self) -> Operand<Self::Basis>;
}

/// Makes each concrete elementary 2D curve its own basis.
macro_rules! own_curve_basis_2d {
    ([$($curve:ty),*]) => {$(
        impl CurveOperand2D for $curve {
            type Basis = $curve;

            #[inline]
            fn basis_and_range(&self) -> Operand<$curve> {
                (self.clone(), None)
            }
        }
    )*};
}

crate::curves::curve_2d_variants!(Elementary, shapes! { own_curve_basis_2d });

impl CurveOperand2D for Curve2D {
    type Basis = CurveKind2D;

    #[inline]
    fn basis_and_range(&self) -> Operand<CurveKind2D> {
        (CurveKind2D::from(self), None)
    }
}

impl CurveOperand2D for BSplineCurve2D {
    type Basis = CurveKind2D;

    #[inline]
    fn basis_and_range(&self) -> Operand<CurveKind2D> {
        (CurveKind2D::Unlisted, None)
    }
}

impl CurveOperand2D for BezierCurve2D {
    type Basis = CurveKind2D;

    #[inline]
    fn basis_and_range(&self) -> Operand<CurveKind2D> {
        (CurveKind2D::Unlisted, None)
    }
}

impl<C: CurveOperand2D> CurveOperand2D for TrimmedCurve2D<C> {
    type Basis = C::Basis;

    /// A trimmed curve of a trimmed curve keeps the inner basis and the
    /// outer range: the outer range lies within the inner one.
    #[inline]
    fn basis_and_range(&self) -> Operand<C::Basis> {
        let (basis, _) = self.basis().basis_and_range();
        (basis, Some((self.first(), self.last())))
    }
}

curve_operands_2d!(curve_operands_2d curve_curve (
    CurveOperand2D, TrimmedCurve2D, intersect_numeric_2d, CurveCurveIntersection2D
));

#[cfg(test)]
mod tests;
