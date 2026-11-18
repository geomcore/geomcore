use crate::base::direction::UNIT_ROUNDING;
use crate::{Direction, Direction2D, Transform, Transform2D};

/// Two-dimensional vector.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector2D {
    /// x-component
    pub x: f64,
    /// y-component
    pub y: f64,
}

impl Vector2D {
    /// Zero vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// assert_eq!(Vector2D::ZERO, Vector2D::new(0.0, 0.0));
    /// ```
    pub const ZERO: Vector2D = Vector2D { x: 0.0, y: 0.0 };

    /// Unit vector along x-axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// assert_eq!(Vector2D::X, Vector2D::new(1.0, 0.0));
    /// ```
    pub const X: Vector2D = Vector2D { x: 1.0, y: 0.0 };

    /// Unit vector along y-axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// assert_eq!(Vector2D::Y, Vector2D::new(0.0, 1.0));
    /// ```
    pub const Y: Vector2D = Vector2D { x: 0.0, y: 1.0 };

    /// Create a new 2D vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// let v = Vector2D::new(3.0, 4.0);
    /// assert_eq!(v.x, 3.0);
    /// assert_eq!(v.y, 4.0);
    /// ```
    pub const fn new(x: f64, y: f64) -> Vector2D {
        Vector2D { x, y }
    }

    /// Dot product.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// assert_eq!(Vector2D::X.dot(Vector2D::Y), 0.0);
    /// assert_eq!(Vector2D::X.dot(Vector2D::X), 1.0);
    /// ```
    pub fn dot(self, other: impl Into<Vector2D>) -> f64 {
        let other: Vector2D = other.into();
        self.x * other.x + self.y * other.y
    }

    /// Cross product (scalar z-component).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// assert_eq!(Vector2D::X.cross(Vector2D::Y), 1.0);
    /// assert_eq!(Vector2D::Y.cross(Vector2D::X), -1.0);
    /// ```
    pub fn cross(self, other: impl Into<Vector2D>) -> f64 {
        let other: Vector2D = other.into();
        self.x * other.y - self.y * other.x
    }

    /// Perpendicular vector (counterclockwise 90 degree rotation).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// assert_eq!(Vector2D::X.perpendicular(), Vector2D::Y);
    /// assert_eq!(Vector2D::Y.perpendicular(), Vector2D::new(-1.0, 0.0));
    /// ```
    pub fn perpendicular(self) -> Vector2D {
        Vector2D {
            x: -self.y,
            y: self.x,
        }
    }

    /// Magnitude (Euclidean length).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// assert_eq!(Vector2D::new(3.0, 4.0).magnitude(), 5.0);
    /// ```
    pub fn magnitude(self) -> f64 {
        self.square_magnitude().sqrt()
    }

    /// Squared magnitude.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// assert_eq!(Vector2D::new(3.0, 4.0).square_magnitude(), 25.0);
    /// ```
    pub fn square_magnitude(self) -> f64 {
        self.x * self.x + self.y * self.y
    }

    /// Whether every coordinate is finite (neither NaN nor infinite).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// assert!(Vector2D::new(1.0, 2.0).is_finite());
    /// assert!(!Vector2D::new(f64::NAN, 2.0).is_finite());
    /// ```
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }

    /// Normalized vector.
    ///
    /// Returns `None` if the squared magnitude is below `f64::MIN_POSITIVE`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// let v = Vector2D::new(3.0, 4.0);
    /// let normalized = v.normalized().unwrap();
    /// assert!((normalized.x - 0.6).abs() < 1e-10);
    /// assert!((normalized.y - 0.8).abs() < 1e-10);
    /// ```
    pub fn normalized(self) -> Option<Vector2D> {
        let sq_mag = self.square_magnitude();
        if sq_mag <= f64::MIN_POSITIVE {
            None
        } else {
            let mag = sq_mag.sqrt();
            Some(Vector2D {
                x: self.x / mag,
                y: self.y / mag,
            })
        }
    }

    /// The unsigned angle between the two vectors, in `[0, π]`, from
    /// `atan2(|a × b|, a · b)`, accurate for nearly parallel vectors too.
    /// It is `0` when either vector is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use geomcore::Vector2D;
    /// assert!((Vector2D::X.angle(Vector2D::new(-1.0, -1.0)) - 3.0 * PI / 4.0).abs() < 1e-15);
    /// ```
    pub fn angle(self, other: impl Into<Vector2D>) -> f64 {
        let other: Vector2D = other.into();
        self.cross(other).abs().atan2(self.dot(other))
    }

    /// The signed angle from this vector to `other`, in `(-π, π]`:
    /// positive when `other` is counterclockwise from this one. It is `0`
    /// when either vector is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use geomcore::Vector2D;
    /// assert!((Vector2D::X.signed_angle(Vector2D::new(0.0, -1.0)) + PI / 2.0).abs() < 1e-15);
    /// ```
    pub fn signed_angle(self, other: Vector2D) -> f64 {
        self.cross(other).atan2(self.dot(other))
    }

    /// The vector turned counterclockwise by `angle` radians.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use geomcore::Vector2D;
    /// let turned = Vector2D::new(2.0, 0.0).rotated(PI / 2.0);
    /// assert!((turned - Vector2D::new(0.0, 2.0)).magnitude() < 1e-15);
    /// ```
    pub fn rotated(self, angle: f64) -> Vector2D {
        let (sin, cos) = angle.sin_cos();
        Vector2D::new(cos * self.x - sin * self.y, sin * self.x + cos * self.y)
    }

    /// Whether the vectors are parallel, in either orientation: the sine of
    /// the angle between them, `|a × b| / (|a| |b|)`, is at most
    /// `angular`. Pass `tol.angular`, or `tol.angular_for_extent(size)`
    /// where a direction error matters as a distance over `size`. A zero
    /// vector is parallel to nothing. Between two directions, which need
    /// no normalizing, ask [`Direction2D::is_parallel`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Tolerance, Vector2D};
    /// let angular = Tolerance::DEFAULT.angular;
    /// assert!(Vector2D::new(2.0, 0.0).is_parallel(Vector2D::new(-5.0, 0.0), angular));
    /// assert!(!Vector2D::X.is_parallel(Vector2D::new(1.0, 1e-6), angular));
    /// assert!(!Vector2D::ZERO.is_parallel(Vector2D::X, angular));
    /// ```
    #[inline]
    pub fn is_parallel(self, other: impl Into<Vector2D>, angular: f64) -> bool {
        let other: Vector2D = other.into();
        let scale = self.magnitude() * other.magnitude();
        scale > 0.0 && self.cross(other).abs() <= angular * scale
    }

    /// Whether the vectors are perpendicular: the cosine of the angle
    /// between them, `|a · b| / (|a| |b|)`, is at most `angular`. A zero
    /// vector is perpendicular to nothing. Between two directions, ask
    /// [`Direction2D::is_perpendicular`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Tolerance, Vector2D};
    /// let angular = Tolerance::DEFAULT.angular;
    /// assert!(Vector2D::new(0.0, 3.0).is_perpendicular(Vector2D::new(2.0, 0.0), angular));
    /// assert!(!Vector2D::ZERO.is_perpendicular(Vector2D::X, angular));
    /// ```
    #[inline]
    pub fn is_perpendicular(self, other: impl Into<Vector2D>, angular: f64) -> bool {
        let other: Vector2D = other.into();
        let scale = self.magnitude() * other.magnitude();
        scale > 0.0 && self.dot(other).abs() <= angular * scale
    }

    /// The part of this vector along `direction`, `d (v · d)`. For a
    /// vector `w` of any length, pass `w.direction()?`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction2D, Vector2D};
    /// let v = Vector2D::new(3.0, 4.0);
    /// assert_eq!(v.component_along(Direction2D::X), Vector2D::new(3.0, 0.0));
    /// ```
    #[inline]
    pub fn component_along(self, direction: Direction2D) -> Vector2D {
        direction * self.dot(direction)
    }

    /// The part of this vector across `direction`: the vector less its
    /// [`component_along`](Self::component_along) `direction`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction2D, Vector2D};
    /// let v = Vector2D::new(3.0, 4.0);
    /// assert_eq!(v.component_perpendicular_to(Direction2D::X), Vector2D::new(0.0, 4.0));
    /// ```
    #[inline]
    pub fn component_perpendicular_to(self, direction: Direction2D) -> Vector2D {
        self - direction * self.dot(direction)
    }

    /// The direction of the vector, a unit vector: the vector itself when it
    /// is already unit within rounding, so that a unit vector costs no square
    /// root or division and comes back bit for bit (a frame or axis rebuilt
    /// from its own directions is the same frame or axis), else the vector
    /// [`normalized`](Self::normalized). `None` for a zero vector, which has no
    /// direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction2D, Vector2D};
    /// assert_eq!(Vector2D::new(3.0, 4.0).direction().unwrap(), Vector2D::new(0.6, 0.8));
    /// assert_eq!(Vector2D::X.direction(), Some(Direction2D::X));
    /// assert_eq!(Vector2D::ZERO.direction(), None);
    /// ```
    #[inline(always)]
    pub fn direction(self) -> Option<Direction2D> {
        if (self.square_magnitude() - 1.0).abs() <= UNIT_ROUNDING {
            Some(Direction2D::from_unit_vector(self))
        } else {
            self.normalized().map(Direction2D::from_unit_vector)
        }
    }

    /// The image of this free vector under `transform`, the same as
    /// [`Transform2D::apply_vector`]: turned, mirrored and scaled, never
    /// translated. A vector may shrink to zero, so this never fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point2D, Transform2D, Vector2D};
    /// let double = Transform2D::scaling(Point2D::new(5.0, 5.0), 2.0);
    /// assert_eq!(Vector2D::new(1.0, -1.0).transformed(double), Vector2D::new(2.0, -2.0));
    /// ```
    pub fn transformed(self, transform: Transform2D) -> Vector2D {
        transform.apply_vector(self)
    }
}

/// Three-dimensional vector.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector {
    /// x-component
    pub x: f64,
    /// y-component
    pub y: f64,
    /// z-component
    pub z: f64,
}

impl Vector {
    /// Zero vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// assert_eq!(Vector::ZERO, Vector::new(0.0, 0.0, 0.0));
    /// ```
    pub const ZERO: Vector = Vector {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// Unit vector along x-axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// assert_eq!(Vector::X, Vector::new(1.0, 0.0, 0.0));
    /// ```
    pub const X: Vector = Vector {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };

    /// Unit vector along y-axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// assert_eq!(Vector::Y, Vector::new(0.0, 1.0, 0.0));
    /// ```
    pub const Y: Vector = Vector {
        x: 0.0,
        y: 1.0,
        z: 0.0,
    };

    /// Unit vector along z-axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// assert_eq!(Vector::Z, Vector::new(0.0, 0.0, 1.0));
    /// ```
    pub const Z: Vector = Vector {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    };

    /// Create a new 3D vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// let v = Vector::new(1.0, 2.0, 3.0);
    /// assert_eq!(v.x, 1.0);
    /// assert_eq!(v.y, 2.0);
    /// assert_eq!(v.z, 3.0);
    /// ```
    pub const fn new(x: f64, y: f64, z: f64) -> Vector {
        Vector { x, y, z }
    }

    /// Dot product.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// assert_eq!(Vector::X.dot(Vector::Y), 0.0);
    /// assert_eq!(Vector::X.dot(Vector::X), 1.0);
    /// ```
    pub fn dot(self, other: impl Into<Vector>) -> f64 {
        let other: Vector = other.into();
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// assert_eq!(Vector::X.cross(Vector::Y), Vector::Z);
    /// assert_eq!(Vector::Y.cross(Vector::Z), Vector::X);
    /// assert_eq!(Vector::Z.cross(Vector::X), Vector::Y);
    /// ```
    pub fn cross(self, other: impl Into<Vector>) -> Vector {
        let other: Vector = other.into();
        Vector {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Magnitude (Euclidean length).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// assert!((Vector::new(1.0, 2.0, 2.0).magnitude() - 3.0).abs() < 1e-10);
    /// ```
    pub fn magnitude(self) -> f64 {
        self.square_magnitude().sqrt()
    }

    /// Squared magnitude.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// assert_eq!(Vector::new(1.0, 2.0, 2.0).square_magnitude(), 9.0);
    /// ```
    pub fn square_magnitude(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Whether every coordinate is finite (neither NaN nor infinite).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// assert!(Vector::new(1.0, 2.0, 3.0).is_finite());
    /// assert!(!Vector::new(f64::NAN, 2.0, 3.0).is_finite());
    /// ```
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Normalized vector.
    ///
    /// Returns `None` if the squared magnitude is below `f64::MIN_POSITIVE`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// let v = Vector::new(3.0, 4.0, 0.0);
    /// let normalized = v.normalized().unwrap();
    /// assert!((normalized.x - 0.6).abs() < 1e-10);
    /// assert!((normalized.y - 0.8).abs() < 1e-10);
    /// assert_eq!(normalized.z, 0.0);
    /// ```
    pub fn normalized(self) -> Option<Vector> {
        let sq_mag = self.square_magnitude();
        if sq_mag <= f64::MIN_POSITIVE {
            None
        } else {
            let mag = sq_mag.sqrt();
            Some(Vector {
                x: self.x / mag,
                y: self.y / mag,
                z: self.z / mag,
            })
        }
    }

    /// The direction of the vector, a unit vector: the vector itself when it
    /// is already unit within rounding, so that a unit vector costs no square
    /// root or division and comes back bit for bit (a frame or axis rebuilt
    /// from its own directions is the same frame or axis), else the vector
    /// [`normalized`](Self::normalized). `None` for a zero vector, which has no
    /// direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction, Vector};
    /// assert_eq!(Vector::new(0.0, 3.0, 4.0).direction().unwrap(), Vector::new(0.0, 0.6, 0.8));
    /// assert_eq!(Vector::Z.direction(), Some(Direction::Z));
    /// assert_eq!(Vector::ZERO.direction(), None);
    /// ```
    #[inline(always)]
    pub fn direction(self) -> Option<Direction> {
        if (self.square_magnitude() - 1.0).abs() <= UNIT_ROUNDING {
            Some(Direction::from_unit_vector(self))
        } else {
            self.normalized().map(Direction::from_unit_vector)
        }
    }

    /// Signed angle from self to other around reference (right-hand rule), in (-π, π].
    ///
    /// Uses `atan2(reference.normalized · (self × other), self · other)`.
    ///
    /// # Preconditions
    ///
    /// The `reference` vector must be non-zero. If `reference` is a zero vector (squared magnitude
    /// ≤ `f64::MIN_POSITIVE`), the function cannot determine the sign and returns the unsigned angle
    /// in [0, π]: `atan2(0.0, self · other)`. This collapses to 0 for `self · other > 0` and π for
    /// `self · other < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use geomcore::Vector;
    /// let angle = Vector::X.angle_with_reference(Vector::Y, Vector::Z);
    /// assert!((angle - PI / 2.0).abs() < 1e-10);
    /// ```
    pub fn angle_with_reference(
        self,
        other: impl Into<Vector>,
        reference: impl Into<Vector>,
    ) -> f64 {
        let other = other.into();
        let reference = reference.into();
        let cross_product = self.cross(other);
        let dot_product = self.dot(other);

        match reference.normalized() {
            Some(unit_reference) => {
                let signed_magnitude = unit_reference.dot(cross_product);
                signed_magnitude.atan2(dot_product)
            }
            None => {
                // Degenerate reference: sign is unknowable; collapse to the unsigned 0/π case
                f64::atan2(0.0, dot_product)
            }
        }
    }

    /// The unsigned angle between the two vectors, in `[0, π]`, from
    /// `atan2(|a × b|, a · b)`, accurate for nearly parallel vectors too.
    /// It is `0` when either vector is zero. For a signed angle about an
    /// axis, see [`angle_with_reference`](Self::angle_with_reference).
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use geomcore::Vector;
    /// assert!((Vector::X.angle(Vector::new(0.0, -2.0, 0.0)) - PI / 2.0).abs() < 1e-15);
    /// ```
    pub fn angle(self, other: impl Into<Vector>) -> f64 {
        let other: Vector = other.into();
        self.cross(other).magnitude().atan2(self.dot(other))
    }

    /// Whether the vectors are parallel, in either orientation: the sine of
    /// the angle between them, `|a × b| / (|a| |b|)`, is at most
    /// `angular`. Pass `tol.angular`, or `tol.angular_for_extent(size)`
    /// where a direction error matters as a distance over `size`. A zero
    /// vector is parallel to nothing. Between two directions, which need
    /// no normalizing, ask [`Direction::is_parallel`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Tolerance, Vector};
    /// let angular = Tolerance::DEFAULT.angular;
    /// assert!(Vector::new(0.0, 0.0, 2.0).is_parallel(-Vector::Z, angular));
    /// assert!(!Vector::Z.is_parallel(Vector::new(1e-6, 0.0, 1.0), angular));
    /// assert!(!Vector::ZERO.is_parallel(Vector::Z, angular));
    /// ```
    #[inline]
    pub fn is_parallel(self, other: impl Into<Vector>, angular: f64) -> bool {
        let other: Vector = other.into();
        let scale = self.magnitude() * other.magnitude();
        scale > 0.0 && self.cross(other).magnitude() <= angular * scale
    }

    /// Whether the vectors are perpendicular: the cosine of the angle
    /// between them, `|a · b| / (|a| |b|)`, is at most `angular`. A zero
    /// vector is perpendicular to nothing. Between two directions, ask
    /// [`Direction::is_perpendicular`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Tolerance, Vector};
    /// let angular = Tolerance::DEFAULT.angular;
    /// assert!(Vector::new(3.0, 0.0, 0.0).is_perpendicular(Vector::Z, angular));
    /// assert!(!Vector::ZERO.is_perpendicular(Vector::Z, angular));
    /// ```
    #[inline]
    pub fn is_perpendicular(self, other: impl Into<Vector>, angular: f64) -> bool {
        let other: Vector = other.into();
        let scale = self.magnitude() * other.magnitude();
        scale > 0.0 && self.dot(other).abs() <= angular * scale
    }

    /// The part of this vector along `direction`, `d (v · d)`. For a
    /// vector `w` of any length, pass `w.direction()?`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction, Vector};
    /// let v = Vector::new(1.0, 2.0, 3.0);
    /// assert_eq!(v.component_along(Direction::Z), Vector::new(0.0, 0.0, 3.0));
    /// ```
    #[inline]
    pub fn component_along(self, direction: Direction) -> Vector {
        direction * self.dot(direction)
    }

    /// The part of this vector across `direction`: the vector less its
    /// [`component_along`](Self::component_along) `direction`. Its length
    /// is the distance from the tip of the vector to the line along
    /// `direction`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction, Vector};
    /// let v = Vector::new(1.0, 2.0, 3.0);
    /// assert_eq!(v.component_perpendicular_to(Direction::Z), Vector::new(1.0, 2.0, 0.0));
    /// ```
    #[inline]
    pub fn component_perpendicular_to(self, direction: Direction) -> Vector {
        self - direction * self.dot(direction)
    }

    /// The image of this free vector under `transform`, the same as
    /// [`Transform::apply_vector`]: turned, mirrored and scaled, never
    /// translated. A vector may shrink to zero, so this never fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Transform, Vector};
    /// let double = Transform::scaling(Point::new(5.0, 5.0, 5.0), 2.0);
    /// assert_eq!(Vector::X.transformed(double), Vector::new(2.0, 0.0, 0.0));
    /// ```
    pub fn transformed(self, transform: Transform) -> Vector {
        transform.apply_vector(self)
    }
}

// Implement arithmetic operations for Vector
impl std::ops::Add for Vector {
    type Output = Vector;

    fn add(self, other: Vector) -> Vector {
        Vector {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl std::ops::Sub for Vector {
    type Output = Vector;

    fn sub(self, other: Vector) -> Vector {
        Vector {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl std::ops::Neg for Vector {
    type Output = Vector;

    fn neg(self) -> Vector {
        Vector {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl std::ops::Mul<f64> for Vector {
    type Output = Vector;

    fn mul(self, scalar: f64) -> Vector {
        Vector {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
        }
    }
}

impl std::ops::Mul<Vector> for f64 {
    type Output = Vector;

    fn mul(self, vector: Vector) -> Vector {
        Vector {
            x: self * vector.x,
            y: self * vector.y,
            z: self * vector.z,
        }
    }
}

/// Each component divided by `scalar`, as IEEE division does it: a zero
/// `scalar` gives infinite or NaN components.
///
/// ```
/// use geomcore::Vector;
/// assert_eq!(Vector::new(2.0, 4.0, 6.0) / 2.0, Vector::new(1.0, 2.0, 3.0));
/// ```
impl std::ops::Div<f64> for Vector {
    type Output = Vector;

    fn div(self, scalar: f64) -> Vector {
        Vector {
            x: self.x / scalar,
            y: self.y / scalar,
            z: self.z / scalar,
        }
    }
}

// Implement arithmetic operations for Vector2D
impl std::ops::Add for Vector2D {
    type Output = Vector2D;

    fn add(self, other: Vector2D) -> Vector2D {
        Vector2D {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl std::ops::Sub for Vector2D {
    type Output = Vector2D;

    fn sub(self, other: Vector2D) -> Vector2D {
        Vector2D {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl std::ops::Neg for Vector2D {
    type Output = Vector2D;

    fn neg(self) -> Vector2D {
        Vector2D {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl std::ops::Mul<f64> for Vector2D {
    type Output = Vector2D;

    fn mul(self, scalar: f64) -> Vector2D {
        Vector2D {
            x: self.x * scalar,
            y: self.y * scalar,
        }
    }
}

impl std::ops::Mul<Vector2D> for f64 {
    type Output = Vector2D;

    fn mul(self, vector: Vector2D) -> Vector2D {
        Vector2D {
            x: self * vector.x,
            y: self * vector.y,
        }
    }
}

#[cfg(test)]
mod tests;

/// Each component divided by `scalar`, as IEEE division does it.
///
/// ```
/// use geomcore::Vector2D;
/// assert_eq!(Vector2D::new(2.0, 4.0) / 2.0, Vector2D::new(1.0, 2.0));
/// ```
impl std::ops::Div<f64> for Vector2D {
    type Output = Vector2D;

    fn div(self, scalar: f64) -> Vector2D {
        Vector2D {
            x: self.x / scalar,
            y: self.y / scalar,
        }
    }
}
