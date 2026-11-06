use crate::base::tolerance::Tolerance;
use crate::base::transform::Transform;
use crate::base::vector::{Vector, Vector2D};

/// Two-dimensional point.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point2D {
    /// x-coordinate
    pub x: f64,
    /// y-coordinate
    pub y: f64,
}

impl Point2D {
    /// Origin (0, 0).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Point2D;
    /// assert_eq!(Point2D::ORIGIN, Point2D::new(0.0, 0.0));
    /// ```
    pub const ORIGIN: Point2D = Point2D { x: 0.0, y: 0.0 };

    /// Create a new 2D point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Point2D;
    /// let p = Point2D::new(1.0, 2.0);
    /// assert_eq!(p.x, 1.0);
    /// assert_eq!(p.y, 2.0);
    /// ```
    pub const fn new(x: f64, y: f64) -> Point2D {
        Point2D { x, y }
    }

    /// Distance to another point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Point2D;
    /// assert_eq!(Point2D::new(0.0, 0.0).distance(Point2D::new(3.0, 4.0)), 5.0);
    /// ```
    pub fn distance(self, other: Point2D) -> f64 {
        let dx = other.x - self.x;
        let dy = other.y - self.y;
        (dx * dx + dy * dy).sqrt()
    }

    /// Whether `other` is within `tol.confusion` of this point. Unlike `==`,
    /// which compares coordinates exactly, this is not transitive.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point2D, Tolerance};
    /// let p = Point2D::new(1.0, 2.0);
    /// let q = Point2D::new(1.0 + 1e-9, 2.0);
    /// assert!(p.coincides(q, Tolerance::DEFAULT));
    /// assert_ne!(p, q);
    /// ```
    pub fn coincides(self, other: Point2D, tol: Tolerance) -> bool {
        self.distance(other) <= tol.confusion
    }

    /// The point halfway between this point and `other`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Point2D;
    /// let m = Point2D::new(0.0, 0.0).midpoint(Point2D::new(2.0, 4.0));
    /// assert_eq!(m, Point2D::new(1.0, 2.0));
    /// ```
    pub fn midpoint(self, other: Point2D) -> Point2D {
        Point2D::new(0.5 * (self.x + other.x), 0.5 * (self.y + other.y))
    }

    /// Whether every coordinate is finite (neither NaN nor infinite).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Point2D;
    /// assert!(Point2D::new(1.0, 2.0).is_finite());
    /// assert!(!Point2D::new(f64::NAN, 2.0).is_finite());
    /// ```
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }
}

/// Three-dimensional point.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point {
    /// x-coordinate
    pub x: f64,
    /// y-coordinate
    pub y: f64,
    /// z-coordinate
    pub z: f64,
}

impl Point {
    /// Origin (0, 0, 0).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Point;
    /// assert_eq!(Point::ORIGIN, Point::new(0.0, 0.0, 0.0));
    /// ```
    pub const ORIGIN: Point = Point {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// Create a new 3D point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Point;
    /// let p = Point::new(1.0, 2.0, 3.0);
    /// assert_eq!(p.x, 1.0);
    /// assert_eq!(p.y, 2.0);
    /// assert_eq!(p.z, 3.0);
    /// ```
    pub const fn new(x: f64, y: f64, z: f64) -> Point {
        Point { x, y, z }
    }

    /// Distance to another point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Point;
    /// assert_eq!(Point::new(1.0, 2.0, 3.0).distance(Point::new(4.0, 6.0, 3.0)), 5.0);
    /// ```
    pub fn distance(self, other: Point) -> f64 {
        let dx = other.x - self.x;
        let dy = other.y - self.y;
        let dz = other.z - self.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Whether `other` is within `tol.confusion` of this point. Unlike `==`,
    /// which compares coordinates exactly, this is not transitive.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Tolerance};
    /// let p = Point::new(1.0, 2.0, 3.0);
    /// let q = Point::new(1.0, 2.0, 3.0 + 1e-9);
    /// assert!(p.coincides(q, Tolerance::DEFAULT));
    /// assert_ne!(p, q);
    /// ```
    pub fn coincides(self, other: Point, tol: Tolerance) -> bool {
        self.distance(other) <= tol.confusion
    }

    /// The point halfway between this point and `other`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Point;
    /// let m = Point::new(0.0, 0.0, 0.0).midpoint(Point::new(2.0, 4.0, -6.0));
    /// assert_eq!(m, Point::new(1.0, 2.0, -3.0));
    /// ```
    pub fn midpoint(self, other: Point) -> Point {
        Point::new(
            0.5 * (self.x + other.x),
            0.5 * (self.y + other.y),
            0.5 * (self.z + other.z),
        )
    }

    /// The image of this point under `transform`, the same as
    /// [`Transform::apply_point`] read in the order of the shapes'
    /// `transformed`. Unlike a shape, a point cannot collapse, so this
    /// never fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Transform, Vector};
    /// let shift = Transform::translation(Vector::new(1.0, 0.0, 0.0));
    /// assert_eq!(Point::ORIGIN.transformed(shift), Point::new(1.0, 0.0, 0.0));
    /// ```
    pub fn transformed(self, transform: Transform) -> Point {
        transform.apply_point(self)
    }

    /// Whether every coordinate is finite (neither NaN nor infinite).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Point;
    /// assert!(Point::new(1.0, 2.0, 3.0).is_finite());
    /// assert!(!Point::new(f64::NAN, 2.0, 3.0).is_finite());
    /// ```
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }
}

// Point + Vector -> Point
impl std::ops::Add<Vector> for Point {
    type Output = Point;

    fn add(self, vector: Vector) -> Point {
        Point {
            x: self.x + vector.x,
            y: self.y + vector.y,
            z: self.z + vector.z,
        }
    }
}

// Point - Vector -> Point
impl std::ops::Sub<Vector> for Point {
    type Output = Point;

    fn sub(self, vector: Vector) -> Point {
        Point {
            x: self.x - vector.x,
            y: self.y - vector.y,
            z: self.z - vector.z,
        }
    }
}

// Point - Point -> Vector
impl std::ops::Sub for Point {
    type Output = Vector;

    fn sub(self, other: Point) -> Vector {
        Vector {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

// Point2D + Vector2D -> Point2D
impl std::ops::Add<Vector2D> for Point2D {
    type Output = Point2D;

    fn add(self, vector: Vector2D) -> Point2D {
        Point2D {
            x: self.x + vector.x,
            y: self.y + vector.y,
        }
    }
}

// Point2D - Vector2D -> Point2D
impl std::ops::Sub<Vector2D> for Point2D {
    type Output = Point2D;

    fn sub(self, vector: Vector2D) -> Point2D {
        Point2D {
            x: self.x - vector.x,
            y: self.y - vector.y,
        }
    }
}

// Point2D - Point2D -> Vector2D
impl std::ops::Sub for Point2D {
    type Output = Vector2D;

    fn sub(self, other: Point2D) -> Vector2D {
        Vector2D {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

#[cfg(test)]
mod tests;
