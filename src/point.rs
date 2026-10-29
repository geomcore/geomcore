use crate::tolerance::Tolerance;
use crate::vector::{Vector, Vector2D};

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
mod tests {
    use super::*;

    #[test]
    fn test_point3_origin() {
        assert_eq!(Point::ORIGIN, Point::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_point3_new() {
        let p = Point::new(1.0, 2.0, 3.0);
        assert_eq!(p.x, 1.0);
        assert_eq!(p.y, 2.0);
        assert_eq!(p.z, 3.0);
    }

    #[test]
    fn test_point3_distance() {
        let p1 = Point::new(1.0, 2.0, 3.0);
        let p2 = Point::new(4.0, 6.0, 3.0);
        assert_eq!(p1.distance(p2), 5.0);
    }

    #[test]
    fn test_point3_coincides_within_confusion() {
        let p = Point::new(1.0, 2.0, 3.0);
        let tol = Tolerance::DEFAULT;
        assert!(p.coincides(Point::new(1.0, 2.0, 3.0 + 5e-8), tol));
        assert!(!p.coincides(Point::new(1.0, 2.0, 3.0 + 5e-7), tol));
        let loose = Tolerance::new(1e-3, tol.angular, tol.parametric);
        assert!(p.coincides(Point::new(1.0, 2.0005, 3.0), loose));
        assert!(!p.coincides(Point::new(f64::NAN, 2.0, 3.0), loose));
    }

    #[test]
    fn test_point3_add_vector() {
        let p = Point::new(1.0, 2.0, 3.0);
        let v = Vector::new(4.0, 5.0, 6.0);
        assert_eq!(p + v, Point::new(5.0, 7.0, 9.0));
    }

    #[test]
    fn test_point3_sub_vector() {
        let p = Point::new(5.0, 7.0, 9.0);
        let v = Vector::new(4.0, 5.0, 6.0);
        assert_eq!(p - v, Point::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_point3_sub_point() {
        let p1 = Point::new(4.0, 6.0, 8.0);
        let p2 = Point::new(1.0, 2.0, 3.0);
        assert_eq!(p1 - p2, Vector::new(3.0, 4.0, 5.0));
    }

    #[test]
    fn test_point3_op_identity() {
        let p = Point::new(1.0, 2.0, 3.0);
        let q = Point::new(4.0, 5.0, 6.0);
        let diff = q - p;
        assert_eq!(p + diff, q);
    }

    #[test]
    fn test_point2_origin() {
        assert_eq!(Point2D::ORIGIN, Point2D::new(0.0, 0.0));
    }

    #[test]
    fn test_point2_new() {
        let p = Point2D::new(1.0, 2.0);
        assert_eq!(p.x, 1.0);
        assert_eq!(p.y, 2.0);
    }

    #[test]
    fn test_point2_distance() {
        let p1 = Point2D::new(0.0, 0.0);
        let p2 = Point2D::new(3.0, 4.0);
        assert_eq!(p1.distance(p2), 5.0);
    }

    #[test]
    fn test_point2_coincides_within_confusion() {
        let p = Point2D::new(1.0, 2.0);
        let tol = Tolerance::DEFAULT;
        assert!(p.coincides(Point2D::new(1.0 + 5e-8, 2.0), tol));
        assert!(!p.coincides(Point2D::new(1.0 + 5e-7, 2.0), tol));
        assert!(!p.coincides(Point2D::new(1.0, f64::INFINITY), tol));
    }

    #[test]
    fn test_point2_add_vector() {
        let p = Point2D::new(1.0, 2.0);
        let v = Vector2D::new(3.0, 4.0);
        assert_eq!(p + v, Point2D::new(4.0, 6.0));
    }

    #[test]
    fn test_point2_sub_vector() {
        let p = Point2D::new(4.0, 6.0);
        let v = Vector2D::new(1.0, 2.0);
        assert_eq!(p - v, Point2D::new(3.0, 4.0));
    }

    #[test]
    fn test_point2_sub_point() {
        let p1 = Point2D::new(4.0, 6.0);
        let p2 = Point2D::new(1.0, 2.0);
        assert_eq!(p1 - p2, Vector2D::new(3.0, 4.0));
    }

    #[test]
    fn test_point2_op_identity() {
        let p = Point2D::new(1.0, 2.0);
        let q = Point2D::new(4.0, 6.0);
        let diff = q - p;
        assert_eq!(p + diff, q);
    }
}
