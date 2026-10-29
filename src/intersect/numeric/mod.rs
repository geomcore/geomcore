//! Numeric intersections: the searches `intersect` runs on a pair or
//! configuration without a closed form. Each result states its
//! completeness.

pub(crate) mod bezier_pieces;
pub(crate) mod curve_curve;
pub(crate) mod curve_curve_hull;
pub(crate) mod curve_surface;
pub(crate) mod curve_surface_hull;
pub(crate) mod embedding_2d;
pub(crate) mod marching;
pub(crate) mod overlap_runs;
pub(crate) mod pairs_with_bspline;
pub(crate) mod patch_pairs;
pub(crate) mod subdivision;
