//! Least-squares fits of planes, circles, spheres and cylinders to
//! points: [`Plane::fit`], [`Circle::fit`], [`Circle2D::fit`],
//! [`Sphere::fit`] and [`Cylinder::fit`], the [`Fit`] they return and the
//! [`FitError`] they refuse with. The numbers are in
//! `math/elementary/fitting.rs`.

use crate::base::tier::TierPoint;
use crate::math::elementary::fitting;
use crate::{
    Circle, Circle2D, Cylinder, Direction, Frame, Frame2D, Plane, Point, Point2D, Sphere,
    Tolerance, Vector, Vector2D,
};
use std::fmt;

/// Why a least-squares fit refused its points.
///
/// The degenerate configurations are judged at the confusion tolerance:
/// the points are collinear when every one of them lies within
/// `tol.confusion` of their least-squares line, and so on, so a set of
/// points the caller's model cannot tell from a line is refused rather
/// than fitted with an arbitrary shape (ADR 0055).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FitError {
    /// Fewer points than the shape needs: as many as its degrees of
    /// freedom, counting two equations per point for a circle in space
    /// (across its plane and within it). A plane, a circle and a circle in
    /// the plane need 3 points, a sphere 4 and a cylinder 5.
    TooFewPoints {
        /// The fewest points the shape needs.
        required: usize,
        /// The points given.
        given: usize,
    },
    /// A point has a NaN or infinite coordinate, or the fitted shape
    /// overflows.
    NonFinite,
    /// Every point lies within `tol.confusion` of their centroid.
    ConfusedPoints,
    /// Every point lies within `tol.confusion` of their least-squares
    /// line, which leaves a plane's normal or a circle's centre
    /// undetermined.
    CollinearPoints,
    /// Every point lies within `tol.confusion` of their least-squares
    /// plane, which leaves a sphere or a cylinder undetermined: a circle
    /// lies on a whole family of spheres, and an ellipse on two cylinders.
    CoplanarPoints,
}

impl fmt::Display for FitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FitError::TooFewPoints { required, given } => write!(
                f,
                "the fit needs at least {required} points, {given} were given"
            ),
            FitError::NonFinite => f.write_str("a point or the fitted shape is NaN or infinite"),
            FitError::ConfusedPoints => f.write_str("the points are confused"),
            FitError::CollinearPoints => f.write_str("the points are collinear"),
            FitError::CoplanarPoints => f.write_str("the points are coplanar"),
        }
    }
}

impl std::error::Error for FitError {}

/// A shape fitted to points by least squares, with how far the points lie
/// from it: the root mean square and the largest of their distances to
/// the shape, its residuals. A metrology caller accepts or rejects the
/// fit on them. The residuals are lengths: they scale with the points.
///
/// [`PlaneFit`], [`CircleFit`], [`CircleFit2D`], [`SphereFit`] and
/// [`CylinderFit`] name its instances.
///
/// # Examples
///
/// ```
/// use geomcore::{Plane, PlaneFit, Point, Tolerance};
/// let points = [
///     Point::new(0.0, 0.0, 0.1),
///     Point::new(1.0, 0.0, -0.1),
///     Point::new(1.0, 1.0, 0.1),
///     Point::new(0.0, 1.0, -0.1),
/// ];
/// let fit: PlaneFit = Plane::fit(&points, Tolerance::DEFAULT).unwrap();
/// assert!((fit.root_mean_square_residual() - 0.1).abs() < 1e-12);
/// assert!((fit.maximum_residual() - 0.1).abs() < 1e-12);
/// let plane: Plane = fit.into_shape();
/// assert!((plane.signed_distance(Point::new(0.5, 0.5, 2.0)).abs() - 2.0).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fit<S> {
    shape: S,
    root_mean_square_residual: f64,
    maximum_residual: f64,
}

/// A plane fitted to points (see [`Fit`] and [`Plane::fit`]).
pub type PlaneFit = Fit<Plane>;

/// A circle in space fitted to points (see [`Fit`] and [`Circle::fit`]).
pub type CircleFit = Fit<Circle>;

/// A circle in the plane fitted to points (see [`Fit`] and
/// [`Circle2D::fit`]).
pub type CircleFit2D = Fit<Circle2D>;

/// A sphere fitted to points (see [`Fit`] and [`Sphere::fit`]).
pub type SphereFit = Fit<Sphere>;

/// A cylinder fitted to points (see [`Fit`] and [`Cylinder::fit`]).
pub type CylinderFit = Fit<Cylinder>;

impl<S> Fit<S> {
    /// The fit made of the parts its getters return: for rebuilding a
    /// stored one, as Python's `pickle` does. The parts are kept as they
    /// are, not checked.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, PlaneFit, Point, Tolerance};
    /// let points = [Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)];
    /// let fit = Plane::fit(&points, Tolerance::DEFAULT).unwrap();
    /// let rebuilt = PlaneFit::from_parts(
    ///     *fit.shape(),
    ///     fit.root_mean_square_residual(),
    ///     fit.maximum_residual(),
    /// );
    /// assert_eq!(rebuilt, fit);
    /// ```
    pub fn from_parts(shape: S, root_mean_square_residual: f64, maximum_residual: f64) -> Self {
        Fit {
            shape,
            root_mean_square_residual,
            maximum_residual,
        }
    }

    /// The fitted shape.
    pub fn shape(&self) -> &S {
        &self.shape
    }

    /// The fitted shape, taken out of the result.
    pub fn into_shape(self) -> S {
        self.shape
    }

    /// The root mean square of the distances from the points to the
    /// shape: `sqrt(sum(d²) / n)`.
    pub fn root_mean_square_residual(&self) -> f64 {
        self.root_mean_square_residual
    }

    /// The largest distance from a point to the shape.
    pub fn maximum_residual(&self) -> f64 {
        self.maximum_residual
    }
}

/// `shape`, with the residuals of the `distances` of the points to it.
fn fitted<S>(shape: S, distances: impl Iterator<Item = f64>) -> Fit<S> {
    let (mut count, mut squares, mut largest) = (0usize, 0.0f64, 0.0f64);
    for distance in distances {
        count += 1;
        squares += distance * distance;
        largest = largest.max(distance);
    }
    Fit::from_parts(shape, (squares / count as f64).sqrt(), largest)
}

/// Refuses fewer than `required` points, or a point that is not finite.
fn check_points<P: TierPoint>(points: &[P], required: usize) -> Result<(), FitError> {
    if points.len() < required {
        return Err(FitError::TooFewPoints {
            required,
            given: points.len(),
        });
    }
    if !points.iter().all(|p| p.is_finite()) {
        return Err(FitError::NonFinite);
    }
    Ok(())
}

/// The largest distance from a point of `points` to the line through
/// `through` along the unit `direction`.
fn largest_distance_to_line(points: &[Point], through: Point, direction: Direction) -> f64 {
    points
        .iter()
        .map(|&p| {
            (p - through)
                .component_perpendicular_to(direction)
                .magnitude()
        })
        .fold(0.0, f64::max)
}

/// Refuses points within `tol.confusion` of their centroid, then of their
/// least-squares line, and, when `refuse_coplanar`, of their
/// least-squares plane.
fn check_spread(
    points: &[Point],
    spread: &fitting::Spread,
    tol: Tolerance,
    refuse_coplanar: bool,
) -> Result<(), FitError> {
    let largest_from_centroid = points
        .iter()
        .map(|&p| p.distance(spread.centroid))
        .fold(0.0, f64::max);
    if largest_from_centroid <= tol.confusion {
        return Err(FitError::ConfusedPoints);
    }
    let line = Direction::from_unit_vector(spread.axes[2]);
    if largest_distance_to_line(points, spread.centroid, line) <= tol.confusion {
        return Err(FitError::CollinearPoints);
    }
    if refuse_coplanar {
        let normal = spread.axes[0];
        let largest_from_plane = points
            .iter()
            .map(|&p| (p - spread.centroid).dot(normal).abs())
            .fold(0.0, f64::max);
        if largest_from_plane <= tol.confusion {
            return Err(FitError::CoplanarPoints);
        }
    }
    Ok(())
}

/// The sum of `(p_i − c) × (p_{i+1} − c)` over the points in their order,
/// the last back to the first, with `c` the `centroid`: twice the area
/// vector of the polygon through them, along which they turn
/// counterclockwise.
fn winding(points: &[Point], centroid: Point) -> Vector {
    let mut previous = points[points.len() - 1] - centroid;
    let mut sum = Vector::ZERO;
    for &p in points {
        let current = p - centroid;
        sum = sum + previous.cross(current);
        previous = current;
    }
    sum
}

/// `direction`, turned the way the points turn counterclockwise about.
fn along_winding(direction: Direction, points: &[Point], centroid: Point) -> Direction {
    if winding(points, centroid).dot(direction) < 0.0 {
        -direction
    } else {
        direction
    }
}

/// The frame at `origin` with z along `z_direction` and x towards
/// `first_point`, or any x when the first point lies on the z axis.
fn frame_towards(origin: Point, z_direction: Direction, first_point: Point) -> Frame {
    Frame::new(origin, z_direction, first_point - origin)
        .or_else(|_| Frame::from_z(origin, z_direction))
        .unwrap_or(Frame::WORLD)
}

impl Plane {
    /// The least-squares plane of `points`, with its residuals: the plane
    /// through their centroid across their axis of least spread, which
    /// minimizes the sum of the squared distances from the points (the
    /// orthogonal regression, closed form).
    ///
    /// The plane's origin is the centroid and its x direction points
    /// towards the first point. Its normal is the one the points, in
    /// their order, turn counterclockwise about, so a loop of points
    /// gives the normal of the right-hand rule.
    ///
    /// # Errors
    ///
    /// [`FitError::TooFewPoints`] below 3 points,
    /// [`FitError::NonFinite`] for a NaN or infinite point, and
    /// [`FitError::ConfusedPoints`] or [`FitError::CollinearPoints`] when
    /// every point lies within `tol.confusion` of their centroid or of
    /// their least-squares line.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Tolerance, Vector};
    /// let points = [
    ///     Point::new(0.0, 0.0, 1.0),
    ///     Point::new(2.0, 0.0, 1.0),
    ///     Point::new(2.0, 2.0, 1.0),
    ///     Point::new(0.0, 2.0, 1.0),
    /// ];
    /// let fit = Plane::fit(&points, Tolerance::DEFAULT).unwrap();
    /// assert_eq!(fit.maximum_residual(), 0.0);
    /// let plane = fit.shape();
    /// assert_eq!(plane.frame().origin(), Point::new(1.0, 1.0, 1.0));
    /// // The points turn counterclockwise about +z.
    /// assert_eq!(plane.frame().z_direction(), Vector::Z);
    /// ```
    pub fn fit(points: &[Point], tol: Tolerance) -> Result<PlaneFit, FitError> {
        check_points(points, 3)?;
        let spread = fitting::spread(points);
        check_spread(points, &spread, tol, false)?;
        let normal = Direction::from_unit_vector(spread.axes[0]);
        let normal = along_winding(normal, points, spread.centroid);
        let plane = Plane::from_frame(frame_towards(spread.centroid, normal, points[0]));
        Ok(fitted(
            plane,
            points
                .iter()
                .map(|&p| (p - spread.centroid).dot(normal).abs()),
        ))
    }
}

impl Circle {
    /// The least-squares circle in space of `points`, with its residuals:
    /// the circle that minimizes the sum of the squared distances from
    /// the points to it.
    ///
    /// It starts from the least-squares plane and the algebraic circle of
    /// the points in that plane, then refines the centre, the normal and
    /// the radius together on the distance itself by Levenberg–Marquardt
    /// steps, so the points need not be coplanar.
    ///
    /// The circle runs through the points in their order: its normal is
    /// the one they turn counterclockwise about, and its parameter 0 is
    /// the first point's direction from the centre.
    ///
    /// # Errors
    ///
    /// [`FitError::TooFewPoints`] below 3 points,
    /// [`FitError::NonFinite`] for a NaN or infinite point or a circle
    /// that overflows, and [`FitError::ConfusedPoints`] or
    /// [`FitError::CollinearPoints`] when every point lies within
    /// `tol.confusion` of their centroid or of their least-squares line.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Tolerance, Vector};
    /// use std::f64::consts::TAU;
    /// let points: Vec<Point> = (0..8)
    ///     .map(|k| {
    ///         let angle = k as f64 * TAU / 8.0;
    ///         Point::new(1.0 + 2.0 * angle.cos(), 2.0 * angle.sin(), 3.0)
    ///     })
    ///     .collect();
    /// let fit = Circle::fit(&points, Tolerance::DEFAULT).unwrap();
    /// assert!(fit.maximum_residual() < 1e-12);
    /// let circle = fit.shape();
    /// assert!((circle.radius() - 2.0).abs() < 1e-12);
    /// assert!(circle.center().distance(Point::new(1.0, 0.0, 3.0)) < 1e-12);
    /// assert!((circle.normal().dot(Vector::Z) - 1.0).abs() < 1e-12);
    /// ```
    pub fn fit(points: &[Point], tol: Tolerance) -> Result<CircleFit, FitError> {
        check_points(points, 3)?;
        let spread = fitting::spread(points);
        check_spread(points, &spread, tol, false)?;
        let (center, normal, radius) =
            fitting::circle(points, &spread).ok_or(FitError::CollinearPoints)?;
        let normal = along_winding(normal, points, spread.centroid);
        let frame = frame_towards(center, normal, points[0]);
        let circle = Circle::from_frame(frame, radius).map_err(|_| FitError::NonFinite)?;
        Ok(fitted(
            circle,
            points.iter().map(|&p| {
                let offset = p - center;
                let height = offset.dot(normal);
                let across = offset.component_perpendicular_to(normal).magnitude() - radius;
                height.hypot(across)
            }),
        ))
    }
}

impl Circle2D {
    /// The least-squares circle in the plane of `points`, with its
    /// residuals: the circle that minimizes the sum of the squared
    /// distances from the points to it.
    ///
    /// It starts from the algebraic circle, whose equation
    /// `x² + y² + D x + E y + F = 0` is linear in its unknowns, then
    /// refines the centre and the radius on the distance itself by
    /// Levenberg–Marquardt steps.
    ///
    /// The circle runs through the points in their order: counterclockwise
    /// (a direct frame) unless the points turn clockwise, and its
    /// parameter 0 is the first point's direction from the centre.
    ///
    /// # Errors
    ///
    /// [`FitError::TooFewPoints`] below 3 points,
    /// [`FitError::NonFinite`] for a NaN or infinite point or a circle
    /// that overflows, and [`FitError::ConfusedPoints`] or
    /// [`FitError::CollinearPoints`] when every point lies within
    /// `tol.confusion` of their centroid or of their least-squares line.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D, Tolerance};
    /// let points = [
    ///     Point2D::new(3.0, 0.0),
    ///     Point2D::new(0.0, 3.0),
    ///     Point2D::new(-3.0, 0.0),
    ///     Point2D::new(0.0, -3.0),
    /// ];
    /// let fit = Circle2D::fit(&points, Tolerance::DEFAULT).unwrap();
    /// let circle = fit.shape();
    /// assert!((circle.radius() - 3.0).abs() < 1e-12);
    /// assert!(circle.center().distance(Point2D::ORIGIN) < 1e-12);
    /// assert!(circle.is_counterclockwise());
    /// ```
    pub fn fit(points: &[Point2D], tol: Tolerance) -> Result<CircleFit2D, FitError> {
        check_points(points, 3)?;
        let spread = fitting::spread_2d(points);
        let largest_from_centroid = points
            .iter()
            .map(|&p| p.distance(spread.centroid))
            .fold(0.0, f64::max);
        if largest_from_centroid <= tol.confusion {
            return Err(FitError::ConfusedPoints);
        }
        let across_line = spread.axes[0];
        let largest_from_line = points
            .iter()
            .map(|&p| (p - spread.centroid).dot(across_line).abs())
            .fold(0.0, f64::max);
        if largest_from_line <= tol.confusion {
            return Err(FitError::CollinearPoints);
        }
        let (center, radius) =
            fitting::circle_2d(points, &spread).ok_or(FitError::CollinearPoints)?;
        let mut previous = points[points.len() - 1] - spread.centroid;
        let mut turning = 0.0;
        for &p in points {
            let current = p - spread.centroid;
            turning += previous.cross(current);
            previous = current;
        }
        let x_direction = (points[0] - center).normalized().unwrap_or(Vector2D::X);
        let y_direction = if turning < 0.0 {
            -x_direction.perpendicular()
        } else {
            x_direction.perpendicular()
        };
        let frame =
            Frame2D::new(center, x_direction, y_direction).map_err(|_| FitError::NonFinite)?;
        let circle = Circle2D::from_frame(frame, radius).map_err(|_| FitError::NonFinite)?;
        Ok(fitted(
            circle,
            points.iter().map(|&p| (p.distance(center) - radius).abs()),
        ))
    }
}

impl Sphere {
    /// The least-squares sphere of `points`, with its residuals: the
    /// sphere that minimizes the sum of the squared distances from the
    /// points to it.
    ///
    /// It starts from the algebraic sphere, whose equation
    /// `x² + y² + z² + D x + E y + F z + G = 0` is linear in its unknowns,
    /// then refines the centre and the radius on the distance itself by
    /// Levenberg–Marquardt steps. The sphere's frame is world-aligned, as
    /// [`Sphere::new`] places it.
    ///
    /// # Errors
    ///
    /// [`FitError::TooFewPoints`] below 4 points,
    /// [`FitError::NonFinite`] for a NaN or infinite point or a sphere
    /// that overflows, and [`FitError::ConfusedPoints`],
    /// [`FitError::CollinearPoints`] or [`FitError::CoplanarPoints`] when
    /// every point lies within `tol.confusion` of their centroid, of their
    /// least-squares line or of their least-squares plane.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere, Tolerance};
    /// let points = [
    ///     Point::new(1.0, 0.0, 0.0),
    ///     Point::new(-1.0, 0.0, 0.0),
    ///     Point::new(0.0, 1.0, 0.0),
    ///     Point::new(0.0, 0.0, 1.0),
    ///     Point::new(0.0, 0.0, -1.0),
    /// ];
    /// let fit = Sphere::fit(&points, Tolerance::DEFAULT).unwrap();
    /// assert!((fit.shape().radius() - 1.0).abs() < 1e-12);
    /// assert!(fit.shape().center().distance(Point::ORIGIN) < 1e-12);
    /// ```
    pub fn fit(points: &[Point], tol: Tolerance) -> Result<SphereFit, FitError> {
        check_points(points, 4)?;
        let spread = fitting::spread(points);
        check_spread(points, &spread, tol, true)?;
        let (center, radius) = fitting::sphere(points, &spread).ok_or(FitError::CoplanarPoints)?;
        let sphere = Sphere::new(center, radius).map_err(|_| FitError::NonFinite)?;
        Ok(fitted(
            sphere,
            points.iter().map(|&p| (p.distance(center) - radius).abs()),
        ))
    }
}

impl Cylinder {
    /// The least-squares cylinder of `points`, with its residuals: the
    /// cylinder that minimizes the sum of the squared distances from the
    /// points to it.
    ///
    /// The axis direction is first searched over a hemisphere of
    /// directions, one every 15 degrees, each scored by how well the
    /// points seen along it fit a circle; the best direction, with that
    /// circle's centre and radius, is then refined together on the
    /// distance itself by Levenberg–Marquardt steps; the three best
    /// directions at least 30 degrees apart are refined, and the closest
    /// fit kept. Points spread over an arc of the cylinder, not only
    /// around it, fit too. A narrow arc sampled by a few points may still
    /// leave the refinement in a local minimum, which its residuals show.
    ///
    /// The cylinder's origin is the point of its axis nearest the
    /// points' centroid, its x direction points towards the first point,
    /// and its axis is the direction the points, in their order, turn
    /// counterclockwise about.
    ///
    /// # Errors
    ///
    /// [`FitError::TooFewPoints`] below 5 points,
    /// [`FitError::NonFinite`] for a NaN or infinite point or a cylinder
    /// that overflows, and [`FitError::ConfusedPoints`],
    /// [`FitError::CollinearPoints`] or [`FitError::CoplanarPoints`] when
    /// every point lies within `tol.confusion` of their centroid, of their
    /// least-squares line or of their least-squares plane.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Tolerance, Vector};
    /// use std::f64::consts::TAU;
    /// // Three rings of six points around the z axis, radius 2.
    /// let points: Vec<Point> = (0..18)
    ///     .map(|k| {
    ///         let angle = k as f64 * TAU / 6.0;
    ///         Point::new(2.0 * angle.cos(), 2.0 * angle.sin(), (k / 6) as f64)
    ///     })
    ///     .collect();
    /// let fit = Cylinder::fit(&points, Tolerance::DEFAULT).unwrap();
    /// assert!(fit.maximum_residual() < 1e-12);
    /// let cylinder = fit.shape();
    /// assert!((cylinder.radius() - 2.0).abs() < 1e-12);
    /// assert!((cylinder.axis().direction().dot(Vector::Z) - 1.0).abs() < 1e-12);
    /// assert!(cylinder.frame().origin().distance(Point::new(0.0, 0.0, 1.0)) < 1e-12);
    /// ```
    pub fn fit(points: &[Point], tol: Tolerance) -> Result<CylinderFit, FitError> {
        check_points(points, 5)?;
        let spread = fitting::spread(points);
        check_spread(points, &spread, tol, true)?;
        let (origin, direction, radius) =
            fitting::cylinder(points, &spread).ok_or(FitError::CoplanarPoints)?;
        let direction = along_winding(direction, points, spread.centroid);
        let frame = frame_towards(origin, direction, points[0]);
        let cylinder = Cylinder::from_frame(frame, radius).map_err(|_| FitError::NonFinite)?;
        Ok(fitted(
            cylinder,
            points.iter().map(|&p| {
                ((p - origin)
                    .component_perpendicular_to(direction)
                    .magnitude()
                    - radius)
                    .abs()
            }),
        ))
    }
}
