//! The domain a numeric search runs over on a curve or on a surface
//! direction, read from the shape's bounds and period.

use crate::geometry::curves::ParametricCurve;
use crate::geometry::surfaces::ParametricSurface;
use crate::math::parameter_domain::ParameterDomain;
use crate::math::parameter_range::ParameterRange;

/// The domain of `curve`.
pub(crate) fn curve_domain<C: ParametricCurve>(curve: &C) -> ParameterDomain {
    ParameterRange::new(curve.bounds(), curve.period()).search_domain()
}

/// The domain of the `u` direction of `surface`.
pub(crate) fn u_domain<S: ParametricSurface + ?Sized>(surface: &S) -> ParameterDomain {
    ParameterRange::new(surface.u_bounds(), surface.u_period()).search_domain()
}

/// The domain of the `v` direction of `surface`.
pub(crate) fn v_domain<S: ParametricSurface + ?Sized>(surface: &S) -> ParameterDomain {
    ParameterRange::new(surface.v_bounds(), surface.v_period()).search_domain()
}
