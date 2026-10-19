//! The poles a B-spline blends: one algebra for the 3D and the 2D tier.

use crate::{Point, Point2D, Vector, Vector2D};

/// A pole an algorithm combines linearly (knot insertion, de Casteljau,
/// degree elevation): a point of a nonrational shape, or a homogeneous
/// point `(w x, w y, (w z,) w)` of a rational one, the weight in the last
/// coordinate.
pub(crate) trait Pole: Copy + Default {
    /// The tier's Cartesian points, a homogeneous pole dehomogenizing to.
    type Cartesian: Copy;

    /// `a self + b other`.
    fn combine(self, a: f64, other: Self, b: f64) -> Self;

    /// `a self`.
    fn scale(self, a: f64) -> Self;

    /// The points of `poles` and their weights, if rational.
    fn parts(poles: &[Self]) -> (Vec<Self::Cartesian>, Option<Vec<f64>>);
}

impl Pole for Point {
    type Cartesian = Point;

    #[inline]
    fn combine(self, a: f64, other: Self, b: f64) -> Self {
        Point::new(
            a * self.x + b * other.x,
            a * self.y + b * other.y,
            a * self.z + b * other.z,
        )
    }

    #[inline]
    fn scale(self, a: f64) -> Self {
        Point::new(a * self.x, a * self.y, a * self.z)
    }

    fn parts(poles: &[Self]) -> (Vec<Point>, Option<Vec<f64>>) {
        (poles.to_vec(), None)
    }
}

impl Pole for Point2D {
    type Cartesian = Point2D;

    #[inline]
    fn combine(self, a: f64, other: Self, b: f64) -> Self {
        Point2D::new(a * self.x + b * other.x, a * self.y + b * other.y)
    }

    #[inline]
    fn scale(self, a: f64) -> Self {
        Point2D::new(a * self.x, a * self.y)
    }

    fn parts(poles: &[Self]) -> (Vec<Point2D>, Option<Vec<f64>>) {
        (poles.to_vec(), None)
    }
}

impl Pole for [f64; 4] {
    type Cartesian = Point;

    #[inline]
    fn combine(self, a: f64, other: Self, b: f64) -> Self {
        std::array::from_fn(|c| a * self[c] + b * other[c])
    }

    #[inline]
    fn scale(self, a: f64) -> Self {
        self.map(|c| a * c)
    }

    fn parts(poles: &[Self]) -> (Vec<Point>, Option<Vec<f64>>) {
        (
            poles.iter().map(|&h| h.cartesian()).collect(),
            Some(poles.iter().map(|&h| h.weight()).collect()),
        )
    }
}

impl Pole for [f64; 3] {
    type Cartesian = Point2D;

    #[inline]
    fn combine(self, a: f64, other: Self, b: f64) -> Self {
        std::array::from_fn(|c| a * self[c] + b * other[c])
    }

    #[inline]
    fn scale(self, a: f64) -> Self {
        self.map(|c| a * c)
    }

    fn parts(poles: &[Self]) -> (Vec<Point2D>, Option<Vec<f64>>) {
        (
            poles.iter().map(|&h| h.cartesian()).collect(),
            Some(poles.iter().map(|&h| h.weight()).collect()),
        )
    }
}

/// The tier's homogeneous pole: blended linearly as any [`Pole`], and read
/// back as a Cartesian point and a weight. As four-wide lane arrays, the
/// coordinates sit in `[0..DIM]`, the weight in the last lane, and the
/// spare coordinate of the 2D tier is zero.
pub(crate) trait HomogeneousPole: Pole<Cartesian: PolePoint<Homogeneous = Self>> {
    /// The homogeneous pole of `p` with the weight `w` folded in.
    fn fold(p: Self::Cartesian, w: f64) -> Self;

    /// The Cartesian point: the coordinates over the weight.
    fn cartesian(self) -> Self::Cartesian;

    /// The weight, the pole's last coordinate.
    fn weight(self) -> f64;

    /// The homogeneous pole as four-wide lanes, coordinates in `[0..DIM]`
    /// and the weight in the last, the spare coordinate zero.
    fn lanes(self) -> [f64; 4];

    /// The homogeneous pole of four-wide lanes as [`Self::lanes`] gives
    /// them, the spare coordinate dropped.
    fn of_lanes(s: [f64; 4]) -> Self;
}

impl HomogeneousPole for [f64; 4] {
    #[inline]
    fn fold(p: Point, w: f64) -> Self {
        [w * p.x, w * p.y, w * p.z, w]
    }

    #[inline]
    fn cartesian(self) -> Point {
        let [x, y, z, w] = self;
        let inverse = 1.0 / w;
        Point::new(x * inverse, y * inverse, z * inverse)
    }

    #[inline]
    fn weight(self) -> f64 {
        self[3]
    }

    #[inline]
    fn lanes(self) -> [f64; 4] {
        self
    }

    #[inline]
    fn of_lanes(s: [f64; 4]) -> Self {
        s
    }
}

impl HomogeneousPole for [f64; 3] {
    #[inline]
    fn fold(p: Point2D, w: f64) -> Self {
        [w * p.x, w * p.y, w]
    }

    #[inline]
    fn cartesian(self) -> Point2D {
        let [x, y, w] = self;
        let inverse = 1.0 / w;
        Point2D::new(x * inverse, y * inverse)
    }

    #[inline]
    fn weight(self) -> f64 {
        self[2]
    }

    #[inline]
    fn lanes(self) -> [f64; 4] {
        [self[0], self[1], 0.0, self[2]]
    }

    #[inline]
    fn of_lanes(s: [f64; 4]) -> Self {
        [s[0], s[1], s[3]]
    }
}

/// The points a B-spline's poles are, and the geometry its evaluation
/// produces: [`Point`] in space, [`Point2D`] in the plane. One
/// implementation per tier, monomorphized, so every loop over coordinates
/// unrolls and no evaluation allocates.
pub(crate) trait PolePoint: Copy + Pole<Cartesian = Self> {
    /// Euclidean coordinates per pole: 3 in space, 2 in the plane.
    const DIM: usize;

    /// The tier's homogeneous pole, the weight folded in.
    type Homogeneous: HomogeneousPole<Cartesian = Self>;

    /// A difference or derivative: the tier's vector type.
    type Vector: Copy
        + std::ops::Sub<Output = Self::Vector>
        + std::ops::Mul<f64, Output = Self::Vector>;

    /// Coordinate `c` of the pole, `c < DIM`.
    fn coord(self, c: usize) -> f64;

    /// The dot product of two of the tier's vectors.
    fn dot(a: Self::Vector, b: Self::Vector) -> f64;

    /// The vector from `Self::ORIGIN` to `self`, the point as a position.
    fn to_vector(self) -> Self::Vector;

    /// The point of homogeneous sums `s`: the coordinates in `[0..DIM]`
    /// and the weight in the last of a four-wide array, the spare
    /// coordinate zero and ignored.
    fn point_of(s: [f64; 4]) -> Self;

    /// The vector of the same sums, seen from the origin: a derivative
    /// before the quotient rule.
    fn vector_of(s: [f64; 4]) -> Self::Vector;
}

impl PolePoint for Point {
    const DIM: usize = 3;

    type Homogeneous = [f64; 4];

    type Vector = Vector;

    #[inline]
    fn coord(self, c: usize) -> f64 {
        match c {
            0 => self.x,
            1 => self.y,
            _ => self.z,
        }
    }

    #[inline]
    fn dot(a: Self::Vector, b: Self::Vector) -> f64 {
        a.dot(b)
    }

    #[inline]
    fn to_vector(self) -> Self::Vector {
        Vector::new(self.x, self.y, self.z)
    }

    #[inline]
    fn point_of(s: [f64; 4]) -> Self {
        Point::new(s[0], s[1], s[2])
    }

    #[inline]
    fn vector_of(s: [f64; 4]) -> Self::Vector {
        Vector::new(s[0], s[1], s[2])
    }
}

impl PolePoint for Point2D {
    const DIM: usize = 2;

    type Homogeneous = [f64; 3];

    type Vector = Vector2D;

    #[inline]
    fn coord(self, c: usize) -> f64 {
        match c {
            0 => self.x,
            _ => self.y,
        }
    }

    #[inline]
    fn dot(a: Self::Vector, b: Self::Vector) -> f64 {
        a.dot(b)
    }

    #[inline]
    fn to_vector(self) -> Self::Vector {
        Vector2D::new(self.x, self.y)
    }

    #[inline]
    fn point_of(s: [f64; 4]) -> Self {
        Point2D::new(s[0], s[1])
    }

    #[inline]
    fn vector_of(s: [f64; 4]) -> Self::Vector {
        Vector2D::new(s[0], s[1])
    }
}

/// The poles of a B-spline, read in place: Cartesian points and, if
/// rational, one weight per point. A surface's poles are row-major, rows
/// along `u`.
#[derive(Clone, Copy)]
pub(crate) struct Poles<'a, P: PolePoint = Point> {
    pub(crate) points: &'a [P],
    pub(crate) weights: Option<&'a [f64]>,
}

impl<P: PolePoint> Poles<'_, P> {
    /// Homogeneous coordinates of pole `i` (`w = 1` if nonrational).
    #[inline]
    pub(crate) fn homogeneous(&self, i: usize) -> P::Homogeneous {
        let p = self.points[i];
        match self.weights {
            Some(weights) => HomogeneousPole::fold(p, weights[i]),
            None => HomogeneousPole::fold(p, 1.0),
        }
    }
}
