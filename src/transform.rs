//! Rigid and affine transformations of points and vectors.

use crate::{Axis, Frame, FrameConstructionError, Point, Vector};

/// An affine transformation of 3D space: a linear map plus a translation.
///
/// Represented as a linear part `linear` (row-major 3x3 matrix, which may
/// encode rotation, scale, and/or reflection) and a translation vector
/// `translation`. A point `p` maps to `linear * p + translation`; a vector
/// `v` maps to `linear * v` (translation-invariant).
///
/// # Examples
///
/// ```
/// use geomcore::{Point, Transform, Vector};
/// let t = Transform::translation(Vector::new(1.0, 2.0, 3.0));
/// assert_eq!(t.apply_point(Point::ORIGIN), Point::new(1.0, 2.0, 3.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform {
    linear: [[f64; 3]; 3],
    translation: Vector,
}

const IDENTITY_LINEAR: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn mat_vec(m: [[f64; 3]; 3], v: Vector) -> Vector {
    Vector::new(
        m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
        m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
        m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
    )
}

fn mat_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    result
}

fn scalar_mul(s: f64, m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] = s * m[i][j];
        }
    }
    result
}

fn mat_sub(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] = a[i][j] - b[i][j];
        }
    }
    result
}

fn mat_add(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    result
}

/// Outer product `d ⊗ d` of a vector with itself, as a 3x3 matrix.
fn outer(d: Vector) -> [[f64; 3]; 3] {
    [
        [d.x * d.x, d.x * d.y, d.x * d.z],
        [d.y * d.x, d.y * d.y, d.y * d.z],
        [d.z * d.x, d.z * d.y, d.z * d.z],
    ]
}

/// Cross-product matrix `[d]×` such that `[d]× * v == d.cross(v)`.
fn cross_matrix(d: Vector) -> [[f64; 3]; 3] {
    [[0.0, -d.z, d.y], [d.z, 0.0, -d.x], [-d.y, d.x, 0.0]]
}

/// Translation `t = P - L·P` that keeps point `p_fixed` fixed under the
/// linear part `linear`.
fn translation_fixing(linear: [[f64; 3]; 3], p_fixed: Point) -> Vector {
    let p = Vector::new(p_fixed.x, p_fixed.y, p_fixed.z);
    p - mat_vec(linear, p)
}

impl Transform {
    /// The identity transformation: leaves every point and vector unchanged.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Transform};
    /// assert_eq!(Transform::IDENTITY.apply_point(Point::new(1.0, 2.0, 3.0)), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub const IDENTITY: Transform = Transform {
        linear: IDENTITY_LINEAR,
        translation: Vector::ZERO,
    };

    /// A pure translation by `offset`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Transform, Vector};
    /// let t = Transform::translation(Vector::new(1.0, 0.0, 0.0));
    /// assert_eq!(t.apply_point(Point::ORIGIN), Point::new(1.0, 0.0, 0.0));
    /// ```
    pub fn translation(offset: Vector) -> Transform {
        Transform {
            linear: IDENTITY_LINEAR,
            translation: offset,
        }
    }

    /// Rotation by `angle` radians about `axis` (right-hand rule).
    ///
    /// Uses the Rodrigues rotation formula
    /// `R = c*I + s*[d]x + (1-c)*d⊗d` for the unit axis direction `d`, then
    /// derives the translation that keeps the axis origin fixed.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Transform, Vector};
    /// use std::f64::consts::PI;
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let r = Transform::rotation(axis, PI / 2.0);
    /// let p = r.apply_point(Point::new(1.0, 0.0, 0.0));
    /// assert!((p.x - 0.0).abs() < 1e-10);
    /// assert!((p.y - 1.0).abs() < 1e-10);
    /// ```
    pub fn rotation(axis: Axis, angle: f64) -> Transform {
        let d = axis.direction();
        let c = angle.cos();
        let s = angle.sin();
        // Rodrigues' rotation formula: R = c*I + s*[d]x + (1-c)*(d⊗d).
        let linear = mat_add(
            mat_add(
                scalar_mul(c, IDENTITY_LINEAR),
                scalar_mul(s, cross_matrix(d)),
            ),
            scalar_mul(1.0 - c, outer(d)),
        );
        let translation = translation_fixing(linear, axis.origin());
        Transform {
            linear,
            translation,
        }
    }

    /// Uniform scaling by `factor` about `center` (`factor` may be negative
    /// or zero).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Transform};
    /// let s = Transform::scaling(Point::new(1.0, 1.0, 1.0), 2.0);
    /// let p = s.apply_point(Point::new(2.0, 2.0, 2.0));
    /// assert!((p.x - 3.0).abs() < 1e-10);
    /// ```
    pub fn scaling(center: Point, factor: f64) -> Transform {
        let linear = scalar_mul(factor, IDENTITY_LINEAR);
        let translation = translation_fixing(linear, center);
        Transform {
            linear,
            translation,
        }
    }

    /// Point reflection (inversion) through `center`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Transform};
    /// let m = Transform::mirror_point(Point::new(1.0, 1.0, 1.0));
    /// let p = m.apply_point(Point::ORIGIN);
    /// assert!((p.x - 2.0).abs() < 1e-10);
    /// assert!((p.y - 2.0).abs() < 1e-10);
    /// assert!((p.z - 2.0).abs() < 1e-10);
    /// ```
    pub fn mirror_point(center: Point) -> Transform {
        let linear = scalar_mul(-1.0, IDENTITY_LINEAR);
        let translation = 2.0 * Vector::new(center.x, center.y, center.z);
        Transform {
            linear,
            translation,
        }
    }

    /// Mirror across `axis` (reflection that fixes the axis line).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Transform, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::X).unwrap();
    /// let m = Transform::mirror_axis(axis);
    /// let p = m.apply_point(Point::new(0.0, 1.0, 0.0));
    /// assert!((p.y - (-1.0)).abs() < 1e-10);
    /// ```
    pub fn mirror_axis(axis: Axis) -> Transform {
        let d = axis.direction();
        let linear = mat_sub(scalar_mul(2.0, outer(d)), IDENTITY_LINEAR);
        let translation = translation_fixing(linear, axis.origin());
        Transform {
            linear,
            translation,
        }
    }

    /// Mirror across the plane through `frame.origin()` with normal
    /// `frame.z_direction()`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Transform};
    /// let m = Transform::mirror_plane(Frame::WORLD);
    /// let p = m.apply_point(Point::new(0.0, 0.0, 5.0));
    /// assert!((p.z - (-5.0)).abs() < 1e-10);
    /// ```
    pub fn mirror_plane(frame: Frame) -> Transform {
        let n = frame.z_direction();
        let linear = mat_sub(IDENTITY_LINEAR, scalar_mul(2.0, outer(n)));
        let translation = translation_fixing(linear, frame.origin());
        Transform {
            linear,
            translation,
        }
    }

    /// Composes `self` with `next`: `self` is applied first, then `next`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Transform, Vector};
    /// use std::f64::consts::PI;
    /// let translate = Transform::translation(Vector::new(1.0, 0.0, 0.0));
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let rotate = Transform::rotation(axis, PI / 2.0);
    /// let composed = translate.then(rotate);
    /// let p = composed.apply_point(Point::ORIGIN);
    /// assert!((p.x - 0.0).abs() < 1e-10);
    /// assert!((p.y - 1.0).abs() < 1e-10);
    /// ```
    pub fn then(self, next: Transform) -> Transform {
        let linear = mat_mul(next.linear, self.linear);
        let translation = mat_vec(next.linear, self.translation) + next.translation;
        Transform {
            linear,
            translation,
        }
    }

    /// Applies this transformation to a point: `linear * p + translation`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Transform, Vector};
    /// let t = Transform::translation(Vector::new(1.0, 2.0, 3.0));
    /// assert_eq!(t.apply_point(Point::ORIGIN), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn apply_point(self, p: Point) -> Point {
        let v = mat_vec(self.linear, Vector::new(p.x, p.y, p.z)) + self.translation;
        Point::new(v.x, v.y, v.z)
    }

    /// Applies this transformation to a vector (translation-invariant):
    /// `linear * v`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Transform, Vector};
    /// let t = Transform::translation(Vector::new(1.0, 2.0, 3.0));
    /// let v = Vector::new(4.0, 5.0, 6.0);
    /// assert_eq!(t.apply_vector(v), v);
    /// ```
    pub fn apply_vector(self, v: Vector) -> Vector {
        mat_vec(self.linear, v)
    }

    /// Applies this transformation to a placement frame, returning the
    /// image frame and the uniform scale factor applied to lengths.
    ///
    /// The linear part of every `Transform` is a uniform scale composed
    /// with an orthogonal map, so the three transformed directions stay
    /// mutually perpendicular and share one scale factor (the magnitude
    /// of the transformed x direction). Each axis is mapped as-is and
    /// divided by that factor, so the image frame is the transform of the
    /// frame: anything placed by the frame (a point at local coordinates
    /// `(a, b, c)`) maps to the same local coordinates in the image. A
    /// mirror (negative determinant) therefore changes the handedness
    /// (see [`Frame::is_direct`]); a rotation keeps it. Fails only when
    /// the transform collapses directions (a zero scale factor).
    pub(crate) fn apply_frame(self, frame: Frame) -> Result<(Frame, f64), FrameConstructionError> {
        let x = self.apply_vector(frame.x_direction());
        let scale = x.magnitude();
        let mapped = Frame::from_axes(
            self.apply_point(frame.origin()),
            x,
            self.apply_vector(frame.y_direction()),
            self.apply_vector(frame.z_direction()),
        )?;
        Ok((mapped, scale))
    }
}

#[cfg(test)]
mod tests;
