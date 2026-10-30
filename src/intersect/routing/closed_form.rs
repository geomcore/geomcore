//! The closed-form pair table: the pairs of shapes the router solves in
//! closed form, and the solver that answers each.
//!
//! Numeric is the default. The router reads every operand as its basis,
//! the shape this table matches on: a concrete elementary curve or surface
//! is its own basis, so a pair whose types are known at compile time
//! reaches its typed solver with no dispatch at run time, and any other
//! operand is read as its kind ([`CurveKind`], [`SurfaceKind`],
//! [`CurveKind2D`]), matched at run time. [`ExactFamily`] is implemented
//! for every pair of bases, from the entries below. A listed pair is
//! answered by its solver, which hands a configuration it does not cover
//! to the numeric search (two cylinders with skew axes). Every other pair
//! goes to the numeric search: a pair the table leaves out, and every pair
//! with a shape it lists nowhere, such as a B-spline or Bézier shape.
//! [`ExactFamily::has_closed_form`] reads the same entries.
//!
//! The concrete types of each list come from the shape enums' variant
//! lists (ADR 0017), so a new shape type edits no list here. Adding a
//! closed form is one solver and one entry: a typed solver under
//! `intersect/exact/`, read into its family result in `exact_family` and
//! listed in its family's `solved` entry, or a solver over kinds listed in
//! the [`SurfaceKind`] pair match, as the plane cuts of `swept_plane` are.

use super::swept_plane::{Swept, swept_plane_intersection};
use crate::Tolerance;
use crate::curves::{
    Circle, Circle2D, Curve, Curve2D, ElementaryCurve, Ellipse, Ellipse2D, Hyperbola, Hyperbola2D,
    Line, Line2D, Parabola, Parabola2D, ParametricCurve,
};
use crate::intersect::{
    CurveCurveIntersection, CurveCurveIntersection2D, CurveSurfaceIntersection, IntersectExact,
    SurfaceSurfaceIntersection,
};
use crate::surfaces::{Cone, Cylinder, ElementarySurface, Plane, Sphere, Torus};

/// The parameter ranges restricting the two operands of a pair, `None`
/// for an operand that is not restricted. A surface is never restricted.
pub(crate) type Ranges = (Option<(f64, f64)>, Option<(f64, f64)>);

/// A pair of bases, as the table answers it: in closed form when it lists
/// the pair, by the numeric search otherwise.
pub(crate) trait ExactFamily<B> {
    /// The family result: [`CurveCurveIntersection`],
    /// [`CurveCurveIntersection2D`], [`CurveSurfaceIntersection`] or
    /// [`SurfaceSurfaceIntersection`].
    type Family;

    /// Whether the table lists the pair: a closed-form solver answers it,
    /// though it may still hand a configuration it does not cover to the
    /// numeric search.
    fn has_closed_form(&self, other: &B) -> bool;

    /// Intersects `self` with `other` in closed form, each restricted to
    /// its range in `ranges`, deciding tangency and coincidence with
    /// `tol`; where the table does not list the pair, or the configuration
    /// has no closed form, the result of `numeric_search` instead.
    // A fallback rather than an `Option`: a curve family result holds its
    // hits inline, and moving it out of an `Option` costs a copy per call.
    fn exact_family(
        &self,
        other: &B,
        ranges: Ranges,
        tol: Tolerance,
        numeric_search: impl FnOnce() -> Self::Family,
    ) -> Self::Family;
}

/// Whether `$value` is one of the listed variants.
macro_rules! is_variant {
    ([$value:expr] $enum:ty; $($variant:path: $shape:ty),+) => {
        matches!($value, $($variant(_))|+)
    };
}

/// The match forwarding to the curve an [`ElementaryCurve`] holds; every
/// variant has its entries, so `$otherwise` is never needed.
macro_rules! each_elementary_curve {
    ($value:expr, $binding:ident => $body:expr; $otherwise:expr) => {
        crate::curves::elementary_curve_dispatch!($value, $binding => $body)
    };
}

/// The match forwarding to the surface an [`ElementarySurface`] holds;
/// every variant has its entries, so `$otherwise` is never needed.
macro_rules! each_elementary_surface {
    ($value:expr, $binding:ident => $body:expr; $otherwise:expr) => {
        crate::surfaces::elementary_surface_dispatch!($value, $binding => $body)
    };
}

/// The match forwarding to the elementary curve a [`Curve2D`] holds, and
/// `$otherwise` for a B-spline or Bézier curve.
macro_rules! each_elementary_curve_2d {
    ($value:expr, $binding:ident => $body:expr; $otherwise:expr) => {
        crate::curves::curve_2d_variants!(
            Elementary,
            forward_or! { $value, $binding => $body; $otherwise }
        )
    };
}

/// The match forwarding to the shape of each listed variant, and
/// `$otherwise` for the variants the list leaves out.
macro_rules! forward_or {
    (
        [$value:expr, $binding:ident => $body:expr; $otherwise:expr]
        $enum:ty;
        $($variant:path: $shape:ty),+
    ) => {
        match $value {
            $($variant($binding) => $body,)+
            _ => $otherwise,
        }
    };
}

/// A curve whose kind the router knows only at run time (a [`Curve`], a
/// B-spline or a Bézier curve), as the table matches it.
#[derive(Clone, Copy)]
pub(crate) enum CurveKind {
    /// An elementary curve, matched on its variant.
    Elementary(ElementaryCurve),
    /// A curve the table pairs with nothing.
    Unlisted,
}

impl From<&Curve> for CurveKind {
    #[inline]
    fn from(curve: &Curve) -> CurveKind {
        curve
            .as_elementary()
            .map_or(CurveKind::Unlisted, CurveKind::Elementary)
    }
}

/// A surface whose kind the router knows only at run time, as the table
/// matches it. Every surface pair is matched on kinds.
#[derive(Clone, Copy)]
pub(crate) enum SurfaceKind<'a> {
    /// An elementary surface, matched on its variant.
    Elementary(ElementarySurface),
    /// A swept surface, which a plane cuts in closed form.
    Swept(Swept<'a>),
    /// A surface the table pairs with nothing.
    Unlisted,
}

impl From<ElementarySurface> for SurfaceKind<'_> {
    #[inline]
    fn from(surface: ElementarySurface) -> Self {
        SurfaceKind::Elementary(surface)
    }
}

/// Reads each concrete elementary surface as its kind.
macro_rules! elementary_surface_kind {
    ([] $enum:ty; $($variant:path: $shape:ty),+) => {$(
        impl From<$shape> for SurfaceKind<'_> {
            #[inline]
            fn from(surface: $shape) -> Self {
                SurfaceKind::Elementary($variant(surface))
            }
        }
    )+};
}

crate::surfaces::surface_variants!(ElementarySurface, elementary_surface_kind! {});

/// A 2D curve whose kind the router knows only at run time (a
/// [`Curve2D`], a B-spline or a Bézier curve), as the table matches it.
#[derive(Clone)]
pub(crate) enum CurveKind2D {
    /// A [`Curve2D`] holding an elementary curve, matched on its variant.
    Elementary(Curve2D),
    /// A curve the table pairs with nothing.
    Unlisted,
}

impl From<&Curve2D> for CurveKind2D {
    #[inline]
    fn from(curve: &Curve2D) -> CurveKind2D {
        if crate::curves::curve_2d_variants!(Elementary, is_variant! { curve }) {
            CurveKind2D::Elementary(curve.clone())
        } else {
            CurveKind2D::Unlisted
        }
    }
}

/// Implements [`ExactFamily`] between each listed basis and an enum on the
/// second side, forwarding to the shape the enum holds.
macro_rules! enum_second {
    ($family:ty: $enum:ty, $each:ident; [$($first:ty),*]) => {$(
        impl ExactFamily<$enum> for $first {
            type Family = $family;

            #[inline]
            fn has_closed_form(&self, other: &$enum) -> bool {
                $each!(other, shape => self.has_closed_form(shape); false)
            }

            #[inline]
            fn exact_family(
                &self,
                other: &$enum,
                ranges: Ranges,
                tol: Tolerance,
                numeric_search: impl FnOnce() -> $family,
            ) -> $family {
                $each!(
                    other,
                    shape => self.exact_family(shape, ranges, tol, numeric_search);
                    numeric_search()
                )
            }
        }
    )*};
}

/// Implements [`ExactFamily`] between an enum on the first side and each
/// listed basis, forwarding from the shape the enum holds.
macro_rules! enum_first {
    ($family:ty: $enum:ty, $each:ident; [$($second:ty),*]) => {$(
        impl ExactFamily<$second> for $enum {
            type Family = $family;

            #[inline]
            fn has_closed_form(&self, other: &$second) -> bool {
                $each!(self, shape => shape.has_closed_form(other); false)
            }

            #[inline]
            fn exact_family(
                &self,
                other: &$second,
                ranges: Ranges,
                tol: Tolerance,
                numeric_search: impl FnOnce() -> $family,
            ) -> $family {
                $each!(
                    self,
                    shape => shape.exact_family(other, ranges, tol, numeric_search);
                    numeric_search()
                )
            }
        }
    )*};
}

/// Implements [`ExactFamily`] between each listed basis and a kind on the
/// second side: forwarded to the elementary shape the kind holds, numeric
/// for any other kind.
macro_rules! kind_second {
    ($family:ty: $kind:ty, $elementary:path; [$($first:ty),*]) => {$(
        impl ExactFamily<$kind> for $first {
            type Family = $family;

            #[inline]
            fn has_closed_form(&self, other: &$kind) -> bool {
                match other {
                    $elementary(shape) => self.has_closed_form(shape),
                    _ => false,
                }
            }

            #[inline]
            fn exact_family(
                &self,
                other: &$kind,
                ranges: Ranges,
                tol: Tolerance,
                numeric_search: impl FnOnce() -> $family,
            ) -> $family {
                match other {
                    $elementary(shape) => self.exact_family(shape, ranges, tol, numeric_search),
                    _ => numeric_search(),
                }
            }
        }
    )*};
}

// Surface–surface: every pair of elementary surfaces through its typed
// solver, and a plane with a swept surface through `swept_plane`. The
// router reads every surface pair as two kinds.

/// Implements [`ExactFamily`] for every pair of concrete elementary
/// surfaces, through its typed solver.
macro_rules! surface_surface {
    ([$($first:ty),*]) => {
        surface_surface!(@rows [$($first),*]; [$($first),*]);
    };
    (@rows [$($first:ty),*]; $seconds:tt) => {$(
        surface_surface!(@row $first; $seconds);
    )*};
    (@row $first:ty; [$($second:ty),*]) => {$(
        impl ExactFamily<$second> for $first {
            type Family = SurfaceSurfaceIntersection;

            #[inline]
            fn has_closed_form(&self, _other: &$second) -> bool {
                true
            }

            #[inline]
            fn exact_family(
                &self,
                other: &$second,
                _ranges: Ranges,
                tol: Tolerance,
                numeric_search: impl FnOnce() -> SurfaceSurfaceIntersection,
            ) -> SurfaceSurfaceIntersection {
                self.intersect_exact(other, tol).into_family().unwrap_or_else(numeric_search)
            }
        }
    )*};
}

crate::surfaces::surface_variants!(ElementarySurface, shapes! { surface_surface });

crate::surfaces::surface_variants!(
    ElementarySurface,
    shapes! { enum_second SurfaceSurfaceIntersection: ElementarySurface, each_elementary_surface; }
);
crate::surfaces::surface_variants!(
    ElementarySurface,
    shapes! { enum_first SurfaceSurfaceIntersection: ElementarySurface, each_elementary_surface; }
);
enum_first!(
    SurfaceSurfaceIntersection: ElementarySurface, each_elementary_surface; [ElementarySurface]
);

impl<'b> ExactFamily<SurfaceKind<'b>> for SurfaceKind<'_> {
    type Family = SurfaceSurfaceIntersection;

    #[inline]
    fn has_closed_form(&self, other: &SurfaceKind<'b>) -> bool {
        match (self, other) {
            (SurfaceKind::Elementary(first), SurfaceKind::Elementary(second)) => {
                first.has_closed_form(second)
            }
            _ => plane_and_swept(*self, *other).is_some(),
        }
    }

    #[inline]
    fn exact_family(
        &self,
        other: &SurfaceKind<'b>,
        ranges: Ranges,
        tol: Tolerance,
        numeric_search: impl FnOnce() -> SurfaceSurfaceIntersection,
    ) -> SurfaceSurfaceIntersection {
        match (self, other) {
            (SurfaceKind::Elementary(first), SurfaceKind::Elementary(second)) => {
                first.exact_family(second, ranges, tol, numeric_search)
            }
            _ => match plane_and_swept(*self, *other) {
                Some((plane, swept)) => {
                    swept_plane_intersection(&plane, swept, tol).unwrap_or_else(numeric_search)
                }
                None => numeric_search(),
            },
        }
    }
}

/// The plane and the swept surface of a pair holding one of each, in
/// either order: the cuts `swept_plane` solves.
#[inline]
fn plane_and_swept<'a>(
    first: SurfaceKind<'a>,
    second: SurfaceKind<'a>,
) -> Option<(Plane, Swept<'a>)> {
    match (first, second) {
        (SurfaceKind::Elementary(ElementarySurface::Plane(plane)), SurfaceKind::Swept(swept))
        | (SurfaceKind::Swept(swept), SurfaceKind::Elementary(ElementarySurface::Plane(plane))) => {
            Some((plane, swept))
        }
        _ => None,
    }
}

// Curve–surface: a line with every elementary surface, a conic with every
// quadric, each through its typed solver.

/// Implements [`ExactFamily`] for each pair of a concrete curve and a
/// concrete surface the table lists (`solved`) through its typed solver,
/// and for each pair of concrete elementary types it leaves out
/// (`unlisted`) through the numeric search: a typed operand needs an
/// answer for every typed partner.
macro_rules! curve_surface {
    (solved: $($curve:ty => [$($surface:ty),*]),*) => {$($(
        impl ExactFamily<$surface> for $curve {
            type Family = CurveSurfaceIntersection;

            #[inline]
            fn has_closed_form(&self, _surface: &$surface) -> bool {
                true
            }

            #[inline]
            fn exact_family(
                &self,
                surface: &$surface,
                (range, _): Ranges,
                tol: Tolerance,
                _numeric_search: impl FnOnce() -> CurveSurfaceIntersection,
            ) -> CurveSurfaceIntersection {
                self.intersect_exact(surface, tol).into_family(self, range, tol)
            }
        }
    )*)*};
    (unlisted: $($curve:ty => $surface:ty),*) => {$(
        impl ExactFamily<$surface> for $curve {
            type Family = CurveSurfaceIntersection;

            #[inline]
            fn has_closed_form(&self, _surface: &$surface) -> bool {
                false
            }

            /// A circle or conic against a torus has no closed form.
            #[inline]
            fn exact_family(
                &self,
                _surface: &$surface,
                _ranges: Ranges,
                _tol: Tolerance,
                numeric_search: impl FnOnce() -> CurveSurfaceIntersection,
            ) -> CurveSurfaceIntersection {
                numeric_search()
            }
        }
    )*};
}

curve_surface!(
    solved:
    Line => [Plane, Cylinder, Cone, Sphere, Torus],
    Circle => [Plane, Cylinder, Cone, Sphere],
    Ellipse => [Plane, Cylinder, Cone, Sphere],
    Parabola => [Plane, Cylinder, Cone, Sphere],
    Hyperbola => [Plane, Cylinder, Cone, Sphere]
);

curve_surface!(
    unlisted:
    Circle => Torus,
    Ellipse => Torus,
    Parabola => Torus,
    Hyperbola => Torus
);

crate::curves::curve_variants!(
    ElementaryCurve,
    shapes! { enum_second CurveSurfaceIntersection: ElementarySurface, each_elementary_surface; }
);
crate::surfaces::surface_variants!(
    ElementarySurface,
    shapes! { enum_first CurveSurfaceIntersection: ElementaryCurve, each_elementary_curve; }
);
enum_first!(CurveSurfaceIntersection: ElementaryCurve, each_elementary_curve; [ElementarySurface]);
crate::curves::curve_variants!(
    ElementaryCurve,
    shapes! { kind_second CurveSurfaceIntersection: SurfaceKind<'_>, SurfaceKind::Elementary; }
);
kind_second!(
    CurveSurfaceIntersection: SurfaceKind<'_>, SurfaceKind::Elementary; [ElementaryCurve]
);

// Curve–curve: every pair of elementary curves, through its typed solver.

/// Implements [`ExactFamily`] for each pair of concrete curves solved in
/// this order (`solved`) through its typed solver, and for each pair
/// solved in the other order (`reversed`) by swapping the operands and
/// the result.
macro_rules! curve_curve {
    ($family:ty; solved: $($first:ty => [$($second:ty),*]),*) => {$($(
        impl ExactFamily<$second> for $first {
            type Family = $family;

            #[inline]
            fn has_closed_form(&self, _other: &$second) -> bool {
                true
            }

            #[inline]
            fn exact_family(
                &self,
                other: &$second,
                ranges: Ranges,
                tol: Tolerance,
                _numeric_search: impl FnOnce() -> $family,
            ) -> $family {
                self.intersect_exact(other, tol).into_family(self, other, ranges, tol)
            }
        }
    )*)*};
    ($family:ty; reversed: $($first:ty => $second:ty),*) => {$(
        impl ExactFamily<$second> for $first {
            type Family = $family;

            #[inline]
            fn has_closed_form(&self, other: &$second) -> bool {
                other.has_closed_form(self)
            }

            #[inline]
            fn exact_family(
                &self,
                other: &$second,
                (first_range, second_range): Ranges,
                tol: Tolerance,
                numeric_search: impl FnOnce() -> $family,
            ) -> $family {
                let numeric_search = || numeric_search().swapped();
                other
                    .exact_family(self, (second_range, first_range), tol, numeric_search)
                    .swapped()
            }
        }
    )*};
}

curve_curve!(
    CurveCurveIntersection;
    solved:
    Line => [Line, Circle, Ellipse, Parabola, Hyperbola],
    Circle => [Circle, Ellipse, Parabola, Hyperbola],
    Ellipse => [Circle, Ellipse, Parabola, Hyperbola],
    Parabola => [Circle, Ellipse, Parabola, Hyperbola],
    Hyperbola => [Circle, Ellipse, Parabola, Hyperbola]
);

curve_curve!(
    CurveCurveIntersection;
    reversed:
    Circle => Line,
    Ellipse => Line,
    Parabola => Line,
    Hyperbola => Line
);

crate::curves::curve_variants!(
    ElementaryCurve,
    shapes! { enum_second CurveCurveIntersection: ElementaryCurve, each_elementary_curve; }
);
crate::curves::curve_variants!(
    ElementaryCurve,
    shapes! { enum_first CurveCurveIntersection: ElementaryCurve, each_elementary_curve; }
);
enum_first!(CurveCurveIntersection: ElementaryCurve, each_elementary_curve; [ElementaryCurve]);
crate::curves::curve_variants!(
    ElementaryCurve,
    shapes! { kind_second CurveCurveIntersection: CurveKind, CurveKind::Elementary; }
);
kind_second!(CurveCurveIntersection: CurveKind, CurveKind::Elementary; [ElementaryCurve]);

impl<B> ExactFamily<B> for CurveKind
where
    ElementaryCurve: ExactFamily<B>,
{
    type Family = <ElementaryCurve as ExactFamily<B>>::Family;

    #[inline]
    fn has_closed_form(&self, other: &B) -> bool {
        match self {
            CurveKind::Elementary(curve) => curve.has_closed_form(other),
            CurveKind::Unlisted => false,
        }
    }

    #[inline]
    fn exact_family(
        &self,
        other: &B,
        ranges: Ranges,
        tol: Tolerance,
        numeric_search: impl FnOnce() -> Self::Family,
    ) -> Self::Family {
        match self {
            CurveKind::Elementary(curve) => curve.exact_family(other, ranges, tol, numeric_search),
            CurveKind::Unlisted => numeric_search(),
        }
    }
}

// 2D curve–curve: every pair of elementary curves, through its typed
// solver.

curve_curve!(
    CurveCurveIntersection2D;
    solved:
    Line2D => [Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D],
    Circle2D => [Circle2D, Ellipse2D, Parabola2D, Hyperbola2D],
    Ellipse2D => [Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D],
    Parabola2D => [Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D],
    Hyperbola2D => [Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D]
);

curve_curve!(CurveCurveIntersection2D; reversed: Circle2D => Line2D);

crate::curves::curve_2d_variants!(
    Elementary,
    shapes! { enum_second CurveCurveIntersection2D: Curve2D, each_elementary_curve_2d; }
);
crate::curves::curve_2d_variants!(
    Elementary,
    shapes! { enum_first CurveCurveIntersection2D: Curve2D, each_elementary_curve_2d; }
);
enum_first!(CurveCurveIntersection2D: Curve2D, each_elementary_curve_2d; [Curve2D]);
crate::curves::curve_2d_variants!(
    Elementary,
    shapes! { kind_second CurveCurveIntersection2D: CurveKind2D, CurveKind2D::Elementary; }
);
kind_second!(CurveCurveIntersection2D: CurveKind2D, CurveKind2D::Elementary; [Curve2D]);

impl<B> ExactFamily<B> for CurveKind2D
where
    Curve2D: ExactFamily<B>,
{
    type Family = <Curve2D as ExactFamily<B>>::Family;

    #[inline]
    fn has_closed_form(&self, other: &B) -> bool {
        match self {
            CurveKind2D::Elementary(curve) => curve.has_closed_form(other),
            CurveKind2D::Unlisted => false,
        }
    }

    #[inline]
    fn exact_family(
        &self,
        other: &B,
        ranges: Ranges,
        tol: Tolerance,
        numeric_search: impl FnOnce() -> Self::Family,
    ) -> Self::Family {
        match self {
            CurveKind2D::Elementary(curve) => {
                curve.exact_family(other, ranges, tol, numeric_search)
            }
            CurveKind2D::Unlisted => numeric_search(),
        }
    }
}

#[cfg(test)]
mod tests;
