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
    pub fn dot(self, other: Vector2D) -> f64 {
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
    pub fn cross(self, other: Vector2D) -> f64 {
        self.x * other.y - self.y * other.x
    }

    /// Perpendicular vector (counterclockwise 90 degree rotation).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// assert_eq!(Vector2D::X.perp(), Vector2D::Y);
    /// assert_eq!(Vector2D::Y.perp(), Vector2D::new(-1.0, 0.0));
    /// ```
    pub fn perp(self) -> Vector2D {
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
    pub fn angle(self, other: Vector2D) -> f64 {
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
    /// vector is parallel to nothing.
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
    pub fn is_parallel(self, other: Vector2D, angular: f64) -> bool {
        let scale = self.magnitude() * other.magnitude();
        scale > 0.0 && self.cross(other).abs() <= angular * scale
    }

    /// Whether the vectors are perpendicular: the cosine of the angle
    /// between them, `|a · b| / (|a| |b|)`, is at most `angular`. A zero
    /// vector is perpendicular to nothing.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Tolerance, Vector2D};
    /// let angular = Tolerance::DEFAULT.angular;
    /// assert!(Vector2D::new(0.0, 3.0).is_perpendicular(Vector2D::new(2.0, 0.0), angular));
    /// assert!(!Vector2D::ZERO.is_perpendicular(Vector2D::X, angular));
    /// ```
    pub fn is_perpendicular(self, other: Vector2D, angular: f64) -> bool {
        let scale = self.magnitude() * other.magnitude();
        scale > 0.0 && self.dot(other).abs() <= angular * scale
    }

    /// The part of this vector along `direction`, `d (v · d) / |d|²`, of
    /// any length; the zero vector when `direction` is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// let v = Vector2D::new(3.0, 4.0);
    /// assert_eq!(v.component_along(Vector2D::new(2.0, 0.0)), Vector2D::new(3.0, 0.0));
    /// ```
    pub fn component_along(self, direction: Vector2D) -> Vector2D {
        let square = direction.square_magnitude();
        if square <= f64::MIN_POSITIVE {
            return Vector2D::ZERO;
        }
        direction * (self.dot(direction) / square)
    }

    /// The part of this vector across `direction`: the vector less its
    /// [`component_along`](Self::component_along) `direction`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector2D;
    /// let v = Vector2D::new(3.0, 4.0);
    /// assert_eq!(v.component_perpendicular_to(Vector2D::X), Vector2D::new(0.0, 4.0));
    /// ```
    pub fn component_perpendicular_to(self, direction: Vector2D) -> Vector2D {
        self - self.component_along(direction)
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
    pub fn dot(self, other: Vector) -> f64 {
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
    pub fn cross(self, other: Vector) -> Vector {
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
    pub fn angle_with_reference(self, other: Vector, reference: Vector) -> f64 {
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
    pub fn angle(self, other: Vector) -> f64 {
        self.cross(other).magnitude().atan2(self.dot(other))
    }

    /// Whether the vectors are parallel, in either orientation: the sine of
    /// the angle between them, `|a × b| / (|a| |b|)`, is at most
    /// `angular`. Pass `tol.angular`, or `tol.angular_for_extent(size)`
    /// where a direction error matters as a distance over `size`. A zero
    /// vector is parallel to nothing.
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
    pub fn is_parallel(self, other: Vector, angular: f64) -> bool {
        let scale = self.magnitude() * other.magnitude();
        scale > 0.0 && self.cross(other).magnitude() <= angular * scale
    }

    /// Whether the vectors are perpendicular: the cosine of the angle
    /// between them, `|a · b| / (|a| |b|)`, is at most `angular`. A zero
    /// vector is perpendicular to nothing.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Tolerance, Vector};
    /// let angular = Tolerance::DEFAULT.angular;
    /// assert!(Vector::new(3.0, 0.0, 0.0).is_perpendicular(Vector::Z, angular));
    /// assert!(!Vector::ZERO.is_perpendicular(Vector::Z, angular));
    /// ```
    pub fn is_perpendicular(self, other: Vector, angular: f64) -> bool {
        let scale = self.magnitude() * other.magnitude();
        scale > 0.0 && self.dot(other).abs() <= angular * scale
    }

    /// The part of this vector along `direction`, `d (v · d) / |d|²`, of
    /// any length; the zero vector when `direction` is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// let v = Vector::new(1.0, 2.0, 3.0);
    /// assert_eq!(v.component_along(Vector::new(0.0, 0.0, 5.0)), Vector::new(0.0, 0.0, 3.0));
    /// ```
    pub fn component_along(self, direction: Vector) -> Vector {
        let square = direction.square_magnitude();
        if square <= f64::MIN_POSITIVE {
            return Vector::ZERO;
        }
        direction * (self.dot(direction) / square)
    }

    /// The part of this vector across `direction`: the vector less its
    /// [`component_along`](Self::component_along) `direction`. Its length
    /// is the distance from the tip of the vector to the line along
    /// `direction`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Vector;
    /// let v = Vector::new(1.0, 2.0, 3.0);
    /// assert_eq!(v.component_perpendicular_to(Vector::Z), Vector::new(1.0, 2.0, 0.0));
    /// ```
    pub fn component_perpendicular_to(self, direction: Vector) -> Vector {
        self - self.component_along(direction)
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
