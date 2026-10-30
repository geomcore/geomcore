//! The expansions of a shape enum's variant list.
//!
//! Each shape enum ([`Curve`](crate::Curve) and
//! [`ElementaryCurve`](crate::ElementaryCurve),
//! [`Curve2D`](crate::Curve2D), [`Surface`](crate::Surface)
//! and [`ElementarySurface`](crate::ElementarySurface)) writes
//! its variant list once, in a `*_variants!` macro beside the enum, and
//! every match that forwards a call to the shape a variant holds, and every
//! `From` impl of a shape, is expanded from that list by the workers here.
//! A new variant is one line in the list; the compiler then points at the
//! matches that do real per-variant work.
//!
//! A list macro calls `worker! { [arguments] Enum; Enum::Variant: Shape, … }`
//! with every path qualified from the crate root, so the expansion does not
//! depend on what the caller imports.

/// Calls `worker!` on a list of variants written as `Variant(Shape)`,
/// qualifying the enum, each variant and each shape with
/// `$crate::$module`. The `*_variants!` list macros are built on it.
macro_rules! expand_variants {
    (
        $module:ident::$enum:ident,
        $($worker:ident)::+ ! { $($arguments:tt)* },
        $($variant:ident($shape:ident)),+ $(,)?
    ) => {
        $($worker)::+! {
            [$($arguments)*]
            $crate::$module::$enum;
            $($crate::$module::$enum::$variant: $crate::$module::$shape),+
        }
    };
}
pub(crate) use expand_variants;

/// The match that binds the shape of every variant to `$binding` and
/// evaluates `$body` with it: the body of a call forwarded to the shape.
macro_rules! dispatch {
    (
        [$value:expr, $binding:ident => $body:expr]
        $enum:ty;
        $($variant:path: $shape:ty),+
    ) => {
        match $value {
            $($variant($binding) => $body,)+
        }
    };
}
pub(crate) use dispatch;

/// `From<Shape> for Enum` for every variant, and with `[&]`
/// `From<&Shape> for Enum` by cloning the shape.
macro_rules! from_shapes {
    ([] $enum:ty; $($variant:path: $shape:ty),+) => {$(
        impl From<$shape> for $enum {
            /// Wraps the shape in its variant.
            #[inline]
            fn from(shape: $shape) -> $enum {
                $variant(shape)
            }
        }
    )+};
    ([&] $enum:ty; $($variant:path: $shape:ty),+) => {$(
        impl From<&$shape> for $enum {
            /// Wraps a clone of the shape in its variant: a copy for an
            /// elementary shape, a deep clone of the poles for a B-spline
            /// or Bézier one.
            #[inline]
            fn from(shape: &$shape) -> $enum {
                $variant(shape.clone())
            }
        }
    )+};
}
pub(crate) use from_shapes;
