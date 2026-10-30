//! UV curves, the curves of a surface's `(u, v)` parameter space, and their
//! mapping to and from 3D curves, which the surface owns as it owns
//! `eval_point` and `parameters_of` at point level (ADR 0030). Curves do
//! not depend on surfaces.

pub(crate) mod to_curve;

pub use to_curve::UvCurveTo3DError;
