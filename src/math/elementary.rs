//! Internal closed-form formulas of the elementary curves and surfaces:
//! evaluation, derivatives and parameter inversion for lines, conics,
//! planes, cylinders, cones, spheres and tori, in 3D and 2D, and their
//! least-squares fits to points.

pub(crate) mod curve;
pub(crate) mod fitting;
pub(crate) mod surface;
