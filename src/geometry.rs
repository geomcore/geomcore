//! The curve and surface types, and the macros that expand their shape
//! enums' variant lists.

/// The curve types: lines, conics, B-spline and Bézier curves.
pub(crate) mod curves;
/// Internal macros that expand the variant list of each shape enum.
pub(crate) mod enum_variants;
/// The surface types: planes, cylinders, cones, spheres, tori, B-spline
/// and Bézier surfaces.
pub(crate) mod surfaces;
