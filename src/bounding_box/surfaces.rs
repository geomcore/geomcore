//! Boxes of surfaces: closed form for spheres and tori, control points for
//! B-splines, the profile's box turned about the axis for revolution
//! surfaces; planes, cylinders, cones and extrusion surfaces are unbounded.

use super::{BoundingBox, HasBoundingBox, coordinates, point};
use crate::Point;
use crate::surfaces::{
    BSplineSurface, BezierSurface, Cone, Cylinder, ElementarySurface, ExtrusionSurface, Plane,
    RevolutionSurface, Sphere, Surface, Torus,
};
use crate::surfaces::{elementary_surface_dispatch, surface_dispatch};

/// The box `center ± half` on each axis.
fn centered(center: Point, half: [f64; 3]) -> BoundingBox {
    let c = coordinates(center);
    BoundingBox {
        min: point(std::array::from_fn(|k| c[k] - half[k])),
        max: point(std::array::from_fn(|k| c[k] + half[k])),
    }
}

impl HasBoundingBox for Plane {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for Cylinder {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for Cone {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for Sphere {
    fn bounding_box(&self) -> Option<BoundingBox> {
        Some(centered(self.center(), [self.radius(); 3]))
    }
}

impl HasBoundingBox for Torus {
    /// Along a unit vector `e`, the torus reaches `R h + r` from its
    /// center, `R` and `r` its major and minor radii and `h` the length of
    /// the projection of `e` onto the torus's plane: the center circle
    /// reaches `R h`, and the tube adds `r` around it.
    fn bounding_box(&self) -> Option<BoundingBox> {
        let frame = self.frame();
        let (x, y) = (frame.x_direction(), frame.y_direction());
        let (x, y) = ([x.x, x.y, x.z], [y.x, y.y, y.z]);
        let (major, minor) = (self.major_radius(), self.minor_radius());
        Some(centered(
            self.center(),
            std::array::from_fn(|k| major * x[k].hypot(y[k]) + minor),
        ))
    }
}

impl HasBoundingBox for BSplineSurface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        BoundingBox::from_points(self.poles().iter().copied())
    }
}

impl HasBoundingBox for BezierSurface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        BoundingBox::from_points(self.poles().iter().copied())
    }
}

impl HasBoundingBox for ExtrusionSurface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for RevolutionSurface {
    /// The profile's box holds the profile, so its corners reach at least
    /// as far from the axis and as far along it: the surface lies in the
    /// cylinder about the axis through the farthest corner, between the
    /// lowest and highest corner. `None` for an unbounded profile.
    fn bounding_box(&self) -> Option<BoundingBox> {
        let profile_box = self.profile().bounding_box()?;
        let (low, high) = (coordinates(profile_box.min), coordinates(profile_box.max));
        let (mut bottom, mut top, mut radius) = (f64::INFINITY, f64::NEG_INFINITY, 0.0_f64);
        for corner in 0..8 {
            let corner = point(std::array::from_fn(|k| {
                if corner & (1 << k) == 0 {
                    low[k]
                } else {
                    high[k]
                }
            }));
            let (height, radial) = self.split(corner);
            bottom = bottom.min(height);
            top = top.max(height);
            radius = radius.max(radial.magnitude());
        }
        let axis = self.axis();
        let a = coordinates(Point::ORIGIN + axis.direction());
        // A disk of radius `r` about the unit axis `a` reaches
        // `r * sqrt(1 - a_k²)` along coordinate `k`.
        let half: [f64; 3] = std::array::from_fn(|k| radius * a[(k + 1) % 3].hypot(a[(k + 2) % 3]));
        let disk = |height: f64| centered(axis.origin() + height * axis.direction(), half);
        Some(disk(bottom).union(disk(top)))
    }
}

impl HasBoundingBox for Surface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        surface_dispatch!(self, s => s.bounding_box())
    }
}

impl HasBoundingBox for ElementarySurface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        elementary_surface_dispatch!(self, s => s.bounding_box())
    }
}
