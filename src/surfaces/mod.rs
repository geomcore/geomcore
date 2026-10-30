//! Public surface types: the five elementary analytic surfaces (plane,
//! cylinder, cone, sphere, torus), the swept surfaces (extrusion and
//! revolution of a profile curve), plus B-spline and Bézier surfaces, built
//! on top of the crate's placement type ([`crate::Frame`]).

mod bezier;
mod bspline;
mod cone;
mod conversion;
mod cylinder;
mod extrusion;
mod fitting;
mod iso;
pub(crate) mod parametric_surface;
mod plane;
mod revolution;
mod sphere;
mod surface;
mod swept;
mod torus;

pub use bezier::BezierSurface;
pub use bspline::BSplineSurface;
pub use cone::{Cone, ConeConstructionError};
pub use cylinder::{Cylinder, CylinderConstructionError};
pub use extrusion::{ExtrusionSurface, ExtrusionSurfaceConstructionError};
pub use iso::SurfaceIsoError;
pub use parametric_surface::{ParametricSurface, ReversibleSurface, SurfaceCurvatures};
pub use plane::{Plane, PlaneConstructionError};
pub use revolution::{RevolutionSurface, RevolutionSurfaceConstructionError};
pub use sphere::{Sphere, SphereConstructionError};
pub use surface::{ElementarySurface, Surface};
pub(crate) use surface::{elementary_surface_dispatch, surface_dispatch, surface_variants};
pub use torus::{Torus, TorusConstructionError};
