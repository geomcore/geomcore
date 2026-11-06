//! Boxes of surfaces: closed form for spheres and tori, control points for
//! B-splines, the profile's box turned about the axis for revolution
//! surfaces; planes, cylinders, cones and extrusion surfaces are unbounded.
//! The tight boxes of B-splines and revolution surfaces come from branch
//! and bound over Bézier pieces (`tight.rs`).

use super::tight::Nets;
use super::{BoundingBox, HasBoundingBox, coordinates, point};
use crate::geometry::curves::{ParametricCurve, RangeOperations};
use crate::geometry::surfaces::ParametricSurface;
use crate::geometry::surfaces::{
    BSplineSurface, BezierSurface, Cone, Cylinder, ElementarySurface, ExtrusionSurface, Plane,
    RevolutionSurface, Sphere, Surface, Torus,
};
use crate::geometry::surfaces::{elementary_surface_dispatch, surface_dispatch};
use crate::{Point, Tolerance};

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

    /// The extremes of each coordinate over the knot-span patches, by
    /// branch and bound.
    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox> {
        Some(Nets::of_surface(self).coordinate_box::<Point>(tol.confusion))
    }
}

impl HasBoundingBox for BezierSurface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        BoundingBox::from_points(self.poles().iter().copied())
    }

    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox> {
        self.as_bspline()?.tight_bounding_box(tol)
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

    /// Along coordinate `k`, the circle a profile point turns on reaches
    /// `o_k + h a_k ± r √(1 − a_k²)`, `h` its height along the unit axis
    /// `a` from the axis origin `o` and `r` its distance from the axis.
    /// Both, the lowest negated, are convex in the point, as `r` is: their
    /// extremes over the profile's Bézier pieces follow by branch and
    /// bound. `None` for an unbounded profile.
    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox> {
        let (first, last) = self.profile().bounds();
        if !(first.is_finite() && last.is_finite()) {
            return None;
        }
        let spline = self.profile().bspline_range(first, last).ok()?;
        let nets = Nets::of_curve(&spline, spline.bounds());
        let axis = self.axis();
        let origin = coordinates(axis.origin());
        let a = coordinates(Point::ORIGIN + axis.direction());
        let reach = |p: [f64; 3], k: usize, side: f64| {
            let d: [f64; 3] = std::array::from_fn(|c| p[c] - origin[c]);
            let height = d[0] * a[0] + d[1] * a[1] + d[2] * a[2];
            let radial = (d[0] - height * a[0])
                .hypot(d[1] - height * a[1])
                .hypot(d[2] - height * a[2]);
            let across = a[(k + 1) % 3].hypot(a[(k + 2) % 3]);
            origin[k] + height * a[k] + side * radial * across
        };
        let high = std::array::from_fn(|k| nets.largest(|p| reach(p, k, 1.0), tol.confusion));
        let low = std::array::from_fn(|k| -nets.largest(|p| -reach(p, k, -1.0), tol.confusion));
        Some(BoundingBox {
            min: point(low),
            max: point(high),
        })
    }
}

impl HasBoundingBox for Surface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        surface_dispatch!(self, s => s.bounding_box())
    }

    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox> {
        surface_dispatch!(self, s => s.tight_bounding_box(tol))
    }
}

impl HasBoundingBox for ElementarySurface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        elementary_surface_dispatch!(self, s => s.bounding_box())
    }
}
