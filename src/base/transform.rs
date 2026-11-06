//! Rigid and affine transformations of points and vectors, and the
//! [`Transformable`] shapes they move.

use crate::{Axis, Frame, FrameConstructionError, Point, Tolerance, Vector};
use std::fmt;

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

    /// Composes `self` with `other`: `self` is applied first, then `other`.
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
    pub fn then(self, other: Transform) -> Transform {
        let linear = mat_mul(other.linear, self.linear);
        let translation = mat_vec(other.linear, self.translation) + other.translation;
        Transform {
            linear,
            translation,
        }
    }

    /// The transformation of a row-major 4×4 homogeneous `matrix`: a point
    /// `p` maps to the upper 3×3 block times `p` plus the last column.
    /// [`to_matrix`](Self::to_matrix) gives it back unchanged.
    ///
    /// # Errors
    ///
    /// [`TransformMatrixError::NonFinite`] for a NaN or infinite entry,
    /// [`TransformMatrixError::NotAffine`] when the last row is not
    /// `[0, 0, 0, 1]`, and [`TransformMatrixError::NotSimilarity`] when the
    /// 3×3 block is not a uniform scale of a rotation or a mirror: its
    /// columns must have one length and be orthogonal within the default
    /// angular tolerance, as [`Frame::from_axes`] asks of its directions.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Transform, TransformMatrixError};
    /// let quarter_turn = [
    ///     [0.0, -1.0, 0.0, 5.0],
    ///     [1.0, 0.0, 0.0, 0.0],
    ///     [0.0, 0.0, 1.0, 0.0],
    ///     [0.0, 0.0, 0.0, 1.0],
    /// ];
    /// let pose = Transform::from_matrix(quarter_turn).unwrap();
    /// assert_eq!(pose.apply_point(Point::new(1.0, 0.0, 0.0)), Point::new(5.0, 1.0, 0.0));
    /// assert_eq!(pose.to_matrix(), quarter_turn);
    ///
    /// let mut sheared = quarter_turn;
    /// sheared[0][2] = 0.5;
    /// assert_eq!(Transform::from_matrix(sheared), Err(TransformMatrixError::NotSimilarity));
    /// ```
    pub fn from_matrix(matrix: [[f64; 4]; 4]) -> Result<Transform, TransformMatrixError> {
        if !matrix.iter().flatten().all(|entry| entry.is_finite()) {
            return Err(TransformMatrixError::NonFinite);
        }
        if matrix[3] != [0.0, 0.0, 0.0, 1.0] {
            return Err(TransformMatrixError::NotAffine);
        }
        let linear = [0, 1, 2].map(|i| [matrix[i][0], matrix[i][1], matrix[i][2]]);
        let column = |j: usize| Vector::new(linear[0][j], linear[1][j], linear[2][j]);
        let columns = [column(0), column(1), column(2)];
        if !is_similarity(columns) {
            return Err(TransformMatrixError::NotSimilarity);
        }
        Ok(Transform {
            linear,
            translation: Vector::new(matrix[0][3], matrix[1][3], matrix[2][3]),
        })
    }

    /// The row-major 4×4 homogeneous matrix of this transformation: the
    /// linear part in the upper 3×3 block, the translation in the last
    /// column, and `[0, 0, 0, 1]` below.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Transform, Vector};
    /// let shift = Transform::translation(Vector::new(1.0, 2.0, 3.0));
    /// assert_eq!(shift.to_matrix()[1], [0.0, 1.0, 0.0, 2.0]);
    /// ```
    pub fn to_matrix(self) -> [[f64; 4]; 4] {
        let [x, y, z] = self.linear;
        let t = self.translation;
        [
            [x[0], x[1], x[2], t.x],
            [y[0], y[1], y[2], t.y],
            [z[0], z[1], z[2], t.z],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }

    /// The transformation that undoes this one: `self.then(inverse)` maps
    /// every point back to itself, within rounding. The linear part is a
    /// uniform scale `s` of an orthogonal map, so its inverse is its
    /// transpose over `s²`.
    ///
    /// # Errors
    ///
    /// [`TransformError::Collapsed`] when the scale factor is zero, so the
    /// transformation collapses space onto a point, and
    /// [`TransformError::NonFinite`] when the inverse overflows `f64`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Transform, TransformError, Vector};
    /// let axis = Axis::new(Point::new(1.0, 0.0, 0.0), Vector::Z).unwrap();
    /// let turn = Transform::rotation(axis, 0.3).then(Transform::scaling(Point::ORIGIN, 2.0));
    /// let back = turn.inverse().unwrap();
    /// let p = Point::new(3.0, -1.0, 2.0);
    /// assert!(back.apply_point(turn.apply_point(p)).distance(p) < 1e-14);
    /// assert_eq!(
    ///     Transform::scaling(Point::ORIGIN, 0.0).inverse(),
    ///     Err(TransformError::Collapsed)
    /// );
    /// ```
    pub fn inverse(self) -> Result<Transform, TransformError> {
        let square_scale = self.apply_vector(Vector::X).square_magnitude();
        if square_scale == 0.0 {
            return Err(TransformError::Collapsed);
        }
        let l = self.linear;
        let linear = [0, 1, 2].map(|i| [0, 1, 2].map(|j| l[j][i] / square_scale));
        let translation = -mat_vec(linear, self.translation);
        let inverse = Transform {
            linear,
            translation,
        };
        let finite =
            linear.iter().flatten().all(|entry| entry.is_finite()) && translation.is_finite();
        if finite {
            Ok(inverse)
        } else {
            Err(TransformError::NonFinite)
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

    /// The rigid motion that takes the `source` placement to the `target`
    /// one: the point at local coordinates `(u, v, w)` in `source` maps to
    /// the point at the same local coordinates in `target`, and
    /// `source.transformed(t)` is `target`. When the two frames have
    /// opposite handedness the motion includes a mirror.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Transform, Vector};
    /// let source = Frame::from_z(Point::new(1.0, 0.0, 0.0), Vector::Z).unwrap();
    /// let target = Frame::new(Point::new(0.0, 2.0, 0.0), Vector::X, Vector::Y).unwrap();
    /// let t = Transform::between_frames(source, target);
    /// let p = t.apply_point(source.point_at(1.0, 2.0, 3.0));
    /// assert!(p.distance(target.point_at(1.0, 2.0, 3.0)) < 1e-12);
    /// ```
    pub fn between_frames(source: Frame, target: Frame) -> Transform {
        let (s, t) = (
            [
                source.x_direction(),
                source.y_direction(),
                source.z_direction(),
            ],
            [
                target.x_direction(),
                target.y_direction(),
                target.z_direction(),
            ],
        );
        // The source directions are orthonormal, so the inverse of the
        // matrix whose columns they are is its transpose: the linear part
        // is the sum of the outer products t_k ⊗ s_k.
        let mut linear = [[0.0; 3]; 3];
        for (sk, tk) in s.iter().zip(t.iter()) {
            let (sk, tk) = ([sk.x, sk.y, sk.z], [tk.x, tk.y, tk.z]);
            for (i, row) in linear.iter_mut().enumerate() {
                for (j, entry) in row.iter_mut().enumerate() {
                    *entry += tk[i] * sk[j];
                }
            }
        }
        let origin = source.origin();
        let translation = target.origin() - Point::ORIGIN - mat_vec(linear, origin - Point::ORIGIN);
        Transform {
            linear,
            translation,
        }
    }

    /// The rigid motion taking `frame` onto the world frame: the image of a
    /// point has, as its coordinates, the point's local coordinates in
    /// `frame` (see [`Frame::local_coordinates`]). Lengths, and so the
    /// parameters of curves, are kept.
    pub(crate) fn world_to_local(frame: Frame) -> Transform {
        let rows = [
            frame.x_direction(),
            frame.y_direction(),
            frame.z_direction(),
        ];
        let linear = rows.map(|r| [r.x, r.y, r.z]);
        let o = frame.origin();
        let translation = -mat_vec(linear, Vector::new(o.x, o.y, o.z));
        Transform {
            linear,
            translation,
        }
    }

    /// The uniform factor by which this transformation scales lengths:
    /// the magnitude of the transformed x direction. One for a rotation, a
    /// translation or a mirror; zero for a transformation that collapses
    /// space onto a point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Transform};
    /// assert_eq!(Transform::scaling(Point::ORIGIN, -2.0).scale_factor(), 2.0);
    /// ```
    pub fn scale_factor(self) -> f64 {
        self.apply_vector(Vector::X).magnitude()
    }

    /// Whether the transformation keeps the handedness of space: true for
    /// rotations, translations and positive scalings, false for mirrors
    /// (the determinant of the linear part is negative). The image of a
    /// surface under a transformation that does not preserve orientation
    /// has its normal on the opposite side of the mirrored normal (see
    /// [`Transformable`]). False as well for a collapsing transformation.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Transform};
    /// assert!(Transform::scaling(Point::ORIGIN, 2.0).preserves_orientation());
    /// assert!(!Transform::mirror_plane(Frame::WORLD).preserves_orientation());
    /// ```
    pub fn preserves_orientation(self) -> bool {
        let m = self.linear;
        let determinant = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        determinant > 0.0
    }

    /// Whether the transformation leaves space unchanged within `tol`:
    /// every entry of the linear part is within `tol.angular` of the
    /// identity's (a turn or a scale change that small) and the translation
    /// is shorter than `tol.confusion`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Tolerance, Transform, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let full_turn = Transform::rotation(axis, std::f64::consts::TAU);
    /// assert!(full_turn.is_identity(Tolerance::DEFAULT));
    /// assert!(!Transform::rotation(axis, 1e-3).is_identity(Tolerance::DEFAULT));
    /// ```
    pub fn is_identity(self, tol: Tolerance) -> bool {
        let linear_is_identity = self
            .linear
            .iter()
            .flatten()
            .zip(IDENTITY_LINEAR.iter().flatten())
            .all(|(entry, identity)| (entry - identity).abs() <= tol.angular);
        linear_is_identity && self.translation.magnitude() <= tol.confusion
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
    /// (see [`Frame::is_direct`]); a rotation keeps it. Fails with
    /// [`TransformError::Collapsed`] when the transform collapses
    /// directions (a zero scale factor) and [`TransformError::NonFinite`]
    /// when the image overflows.
    pub(crate) fn apply_frame(self, frame: Frame) -> Result<(Frame, f64), TransformError> {
        let x = self.apply_vector(frame.x_direction());
        let scale = x.magnitude();
        let mapped = Frame::from_axes(
            self.apply_point(frame.origin()),
            x,
            self.apply_vector(frame.y_direction()),
            self.apply_vector(frame.z_direction()),
        )
        .map_err(frame_image_error)?;
        Ok((mapped, scale))
    }
}

/// Whether three columns are those of a uniform scale of an orthogonal
/// map: all zero, or of one length and orthogonal within the default
/// angular tolerance, as [`Frame::from_axes`] checks its directions.
fn is_similarity(columns: [Vector; 3]) -> bool {
    let lengths = columns.map(Vector::magnitude);
    let longest = lengths[0].max(lengths[1]).max(lengths[2]);
    if longest == 0.0 {
        return true;
    }
    let angular = crate::Tolerance::DEFAULT.angular;
    let one_length = lengths
        .iter()
        .all(|&length| (longest - length) <= angular * longest);
    let orthogonal = [(0, 1), (0, 2), (1, 2)]
        .iter()
        .all(|&(i, j)| columns[i].dot(columns[j]).abs() <= angular * longest * longest);
    one_length && orthogonal
}

/// Error returned by [`Transform::from_matrix`] for a matrix that is not
/// a transformation this crate represents.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransformMatrixError {
    /// An entry is NaN or infinite.
    NonFinite,
    /// The last row is not `[0, 0, 0, 1]`, so the matrix is projective
    /// rather than affine.
    NotAffine,
    /// The 3×3 block is not a uniform scale of a rotation or a mirror: it
    /// shears, or scales the directions by different factors.
    NotSimilarity,
}

impl fmt::Display for TransformMatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            TransformMatrixError::NonFinite => "a matrix entry is NaN or infinite",
            TransformMatrixError::NotAffine => "the last row of the matrix is not [0, 0, 0, 1]",
            TransformMatrixError::NotSimilarity => {
                "the matrix is not a uniform scale of a rotation or a mirror"
            }
        })
    }
}

impl std::error::Error for TransformMatrixError {}

/// Error returned when a shape cannot be moved by a [`Transform`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransformError {
    /// The transform has a zero scale factor, so it collapses the shape
    /// onto a point.
    Collapsed,
    /// The image overflows `f64`: a coordinate or a length is no longer
    /// finite.
    NonFinite,
}

impl fmt::Display for TransformError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            TransformError::Collapsed => "the transform collapses the shape onto a point",
            TransformError::NonFinite => "the transformed shape overflows f64",
        })
    }
}

impl std::error::Error for TransformError {}

/// The error of an image frame or axis that could not be built: an
/// overflow, or directions collapsed by a zero scale factor.
pub(crate) fn frame_image_error(error: FrameConstructionError) -> TransformError {
    match error {
        FrameConstructionError::NonFinite => TransformError::NonFinite,
        _ => TransformError::Collapsed,
    }
}

/// The error of an image that its type's constructor rejected: an
/// overflow when the constructor reports `non_finite`, a collapse (a
/// length scaled down to zero) otherwise.
pub(crate) fn image_error<E: PartialEq>(non_finite: E) -> impl FnOnce(E) -> TransformError {
    move |error| {
        if error == non_finite {
            TransformError::NonFinite
        } else {
            TransformError::Collapsed
        }
    }
}

/// A shape that a [`Transform`] moves: every 3D curve and surface, and
/// their enums.
///
/// The image is parametrized by the image of the shape's parametrization:
/// under a rotation, a translation or a mirror,
/// `shape.transformed(t)?.point_at(u)` is `t.apply_point(shape.point_at(u))`
/// at every parameter. Under a scale `s`, a parameter that is a length
/// (along a line, a parabola, the generators of a plane, cylinder or cone,
/// the sweep of an extrusion) scales by `s` too; angles and B-spline
/// parameters do not. Every frame axis maps through the transform, so a
/// mirror turns a direct frame indirect and the image's normal is the
/// opposite of the mirrored normal; reversing the image restores the
/// mirrored side (see [`ReversibleSurface`](crate::ReversibleSurface)).
///
/// # Errors
///
/// Every shape fails the same way: with [`TransformError::Collapsed`] when
/// the transform has a zero scale factor, and with
/// [`TransformError::NonFinite`] when the image overflows `f64`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Curve, Point, Transform, TransformError, Vector};
///
/// let circle = Curve::from(Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap());
/// let shift = Transform::translation(Vector::new(0.0, 0.0, 5.0));
/// let moved = circle.transformed(shift).unwrap();
/// assert_eq!(moved.point_at(0.0), Point::new(2.0, 0.0, 5.0));
///
/// let collapse = Transform::scaling(Point::ORIGIN, 0.0);
/// assert_eq!(circle.transformed(collapse), Err(TransformError::Collapsed));
/// ```
pub trait Transformable: Sized {
    /// Returns the image of the shape under `transform`.
    ///
    /// # Errors
    ///
    /// [`TransformError::Collapsed`] for a zero scale factor and
    /// [`TransformError::NonFinite`] when the image overflows `f64`.
    fn transformed(&self, transform: Transform) -> Result<Self, TransformError>;
}

#[cfg(test)]
mod tests;
