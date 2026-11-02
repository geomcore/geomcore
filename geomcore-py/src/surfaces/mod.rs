//! The `geomcore.surfaces` classes, one file per class as in the core crate's
//! `surfaces/`.

mod bezier;
mod bspline;
mod cone;
mod cylinder;
mod extrusion;
mod plane;
mod revolution;
mod sphere;
mod torus;
mod uv_curve;

pub(crate) use bezier::PyBezierSurface;
pub(crate) use bspline::PyBSplineSurface;
pub(crate) use cone::PyCone;
pub(crate) use cylinder::PyCylinder;
pub(crate) use extrusion::PyExtrusionSurface;
pub(crate) use plane::PyPlane;
pub(crate) use revolution::PyRevolutionSurface;
pub(crate) use sphere::PySphere;
pub(crate) use torus::PyTorus;
pub(crate) use uv_curve::py_curve_3d_to_uv_numeric;
