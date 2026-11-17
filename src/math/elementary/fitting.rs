//! Least-squares fits of the elementary shapes to points: the numbers
//! under the plane, circle, sphere and cylinder fits.
//!
//! Every fit works on the points moved to their centroid and scaled to a
//! unit root mean square distance from it, so its numbers are near unit
//! size wherever the points sit and whatever their size (ADR 0051), and
//! maps its answer back. The plane is the covariance's axis of least
//! spread, which is already the geometric least-squares fit. A circle or
//! a sphere starts from the algebraic fit, linear in its unknowns, and is
//! refined on the geometric distance by Levenberg–Marquardt steps. A
//! cylinder's axis direction is first searched over a hemisphere of
//! directions, each scored by the algebraic circle fit of the points seen
//! along it, and the best one is refined the same way.

use crate::math::numeric::eigen::symmetric_eigen;
use crate::math::numeric::least_squares::levenberg_marquardt;
use crate::math::numeric::solve::solve_square;
use crate::{Direction, Point, Point2D, Vector, Vector2D};

/// Where points in space sit and how they spread.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Spread {
    /// The centroid.
    pub(crate) centroid: Point,
    /// The root mean square of the distances from the centroid.
    pub(crate) scale: f64,
    /// The unit axes of the covariance, from the least spread to the most:
    /// `axes[0]` is the normal of the least-squares plane, `axes[2]` the
    /// direction of the least-squares line.
    pub(crate) axes: [Vector; 3],
}

/// Where points in the plane sit and how they spread.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Spread2D {
    /// The centroid.
    pub(crate) centroid: Point2D,
    /// The root mean square of the distances from the centroid.
    pub(crate) scale: f64,
    /// The unit axes of the covariance, from the least spread to the most:
    /// `axes[1]` is the direction of the least-squares line.
    pub(crate) axes: [Vector2D; 2],
}

/// The spread of `points`, at least one.
pub(crate) fn spread(points: &[Point]) -> Spread {
    let count = points.len() as f64;
    let sum = points
        .iter()
        .fold(Vector::ZERO, |sum, &p| sum + (p - Point::ORIGIN));
    let centroid = Point::ORIGIN + sum / count;
    let mut covariance = [[0.0; 3]; 3];
    for &p in points {
        let d = p - centroid;
        let d = [d.x, d.y, d.z];
        for i in 0..3 {
            for j in i..3 {
                covariance[i][j] += d[i] * d[j];
            }
        }
    }
    let trace = covariance[0][0] + covariance[1][1] + covariance[2][2];
    let (_, vectors) = symmetric_eigen(covariance);
    Spread {
        centroid,
        scale: (trace / count).sqrt(),
        axes: vectors.map(|[x, y, z]| Vector::new(x, y, z)),
    }
}

/// The spread of `points` in the plane, at least one.
pub(crate) fn spread_2d(points: &[Point2D]) -> Spread2D {
    let count = points.len() as f64;
    let sum = points
        .iter()
        .fold(Vector2D::ZERO, |sum, &p| sum + (p - Point2D::ORIGIN));
    let centroid = Point2D::ORIGIN + sum / count;
    let mut covariance = [[0.0; 2]; 2];
    for &p in points {
        let d = p - centroid;
        covariance[0][0] += d.x * d.x;
        covariance[0][1] += d.x * d.y;
        covariance[1][1] += d.y * d.y;
    }
    let trace = covariance[0][0] + covariance[1][1];
    let (_, vectors) = symmetric_eigen(covariance);
    Spread2D {
        centroid,
        scale: (trace / count).sqrt(),
        axes: vectors.map(|[x, y]| Vector2D::new(x, y)),
    }
}

/// The algebraic circle through plane points `(x, y)`: the centre and the
/// radius of `x² + y² + D x + E y + F = 0` with the least sum of squared
/// left-hand sides, and that sum divided by four times the squared
/// radius, which is near the sum of squared distances to the circle.
/// `None` when the points leave the circle undetermined.
fn algebraic_circle(points: impl Iterator<Item = (f64, f64)>) -> Option<([f64; 2], f64, f64)> {
    let mut normal = [[0.0; 3]; 3];
    let mut right = [0.0; 3];
    let mut squares = 0.0;
    for (x, y) in points {
        let row = [x, y, 1.0];
        let z = x * x + y * y;
        for i in 0..3 {
            for j in 0..3 {
                normal[i][j] += row[i] * row[j];
            }
            right[i] -= z * row[i];
        }
        squares += z * z;
    }
    let [d, e, f] = solve_square(normal, right, 0.0)?;
    let center = [-0.5 * d, -0.5 * e];
    let squared_radius = center[0] * center[0] + center[1] * center[1] - f;
    if squared_radius.is_nan() || squared_radius <= 0.0 {
        return None;
    }
    let algebraic = squares - (right[0] * d + right[1] * e + right[2] * f);
    Some((
        center,
        squared_radius.sqrt(),
        algebraic.max(0.0) / (4.0 * squared_radius),
    ))
}

/// The geometric least-squares circle of plane `points`, at least three
/// and not collinear, spread as `spread`: its centre and radius.
pub(crate) fn circle_2d(points: &[Point2D], spread: &Spread2D) -> Option<(Point2D, f64)> {
    let local = |k: usize| (points[k] - spread.centroid) / spread.scale;
    let (center, radius, _) = algebraic_circle((0..points.len()).map(|k| {
        let q = local(k);
        (q.x, q.y)
    }))?;
    let [x, y, radius] = levenberg_marquardt(
        [center[0], center[1], radius],
        points.len(),
        |&[x, y, radius]: &[f64; 3], k| {
            let offset = local(k) - Vector2D::new(x, y);
            let distance = offset.magnitude();
            let outward = if distance > 0.0 {
                offset / distance
            } else {
                Vector2D::ZERO
            };
            [(distance - radius, [-outward.x, -outward.y, -1.0])]
        },
        |&[x, y, radius], [dx, dy, dr]| [x + dx, y + dy, radius + dr],
    );
    let center = spread.centroid + Vector2D::new(x, y) * spread.scale;
    Some((center, radius.abs() * spread.scale))
}

/// The geometric least-squares sphere of `points`, at least four and not
/// coplanar, spread as `spread`: its centre and radius.
pub(crate) fn sphere(points: &[Point], spread: &Spread) -> Option<(Point, f64)> {
    let local = |k: usize| (points[k] - spread.centroid) / spread.scale;
    // The algebraic sphere x² + y² + z² + D x + E y + F z + G = 0.
    let mut normal = [[0.0; 4]; 4];
    let mut right = [0.0; 4];
    for k in 0..points.len() {
        let q = local(k);
        let row = [q.x, q.y, q.z, 1.0];
        let squared = q.square_magnitude();
        for i in 0..4 {
            for j in 0..4 {
                normal[i][j] += row[i] * row[j];
            }
            right[i] -= squared * row[i];
        }
    }
    let [d, e, f, g] = solve_square(normal, right, 0.0)?;
    let center = Vector::new(-0.5 * d, -0.5 * e, -0.5 * f);
    let squared_radius = center.square_magnitude() - g;
    if squared_radius.is_nan() || squared_radius <= 0.0 {
        return None;
    }
    let [x, y, z, radius] = levenberg_marquardt(
        [center.x, center.y, center.z, squared_radius.sqrt()],
        points.len(),
        |&[x, y, z, radius]: &[f64; 4], k| {
            let offset = local(k) - Vector::new(x, y, z);
            let distance = offset.magnitude();
            let outward = if distance > 0.0 {
                offset / distance
            } else {
                Vector::ZERO
            };
            [(
                distance - radius,
                [-outward.x, -outward.y, -outward.z, -1.0],
            )]
        },
        |&[x, y, z, radius], [dx, dy, dz, dr]| [x + dx, y + dy, z + dz, radius + dr],
    );
    let center = spread.centroid + Vector::new(x, y, z) * spread.scale;
    Some((center, radius.abs() * spread.scale))
}

/// A unit direction with two unit vectors across it, all three
/// orthogonal: the local frame a fit turns the direction in.
#[derive(Debug, Clone, Copy)]
struct TurningDirection {
    direction: Direction,
    across: [Vector; 2],
}

impl TurningDirection {
    /// The direction turned by `first` along `across[0]` and `second`
    /// along `across[1]`, its vectors across carried along: each kept
    /// across the new direction, the one that stays further from it
    /// first.
    fn turned(&self, first: f64, second: f64) -> TurningDirection {
        let moved = *self.direction + self.across[0] * first + self.across[1] * second;
        let Some(direction) = moved.direction() else {
            return *self;
        };
        let [u, v] = self.across.map(|a| a.component_perpendicular_to(direction));
        let kept = if u.square_magnitude() >= v.square_magnitude() {
            u
        } else {
            v
        };
        let Some(kept) = kept.normalized() else {
            return *self;
        };
        let other = direction.cross(kept);
        let across = if u.square_magnitude() >= v.square_magnitude() {
            [kept, other]
        } else {
            [-other, kept]
        };
        TurningDirection { direction, across }
    }
}

/// A circle in space while it is fitted: centre, unit normal and radius.
#[derive(Debug, Clone, Copy)]
struct Ring {
    center: Vector,
    normal: TurningDirection,
    radius: f64,
}

/// The geometric least-squares circle in space of `points`, at least
/// three and not collinear, spread as `spread`: its centre, its unit
/// normal (either way) and its radius. Each point's distance to the
/// circle has two parts, across the circle's plane and across its
/// cylinder, each a residual.
pub(crate) fn circle(points: &[Point], spread: &Spread) -> Option<(Point, Direction, f64)> {
    let local = |k: usize| (points[k] - spread.centroid) / spread.scale;
    let [normal, first, second] = spread.axes;
    let (center, radius, _) = algebraic_circle((0..points.len()).map(|k| {
        let q = local(k);
        (q.dot(first), q.dot(second))
    }))?;
    let start = Ring {
        center: first * center[0] + second * center[1],
        normal: TurningDirection {
            direction: Direction::from_unit_vector(normal),
            across: [first, second],
        },
        radius,
    };
    let fitted = levenberg_marquardt(
        start,
        points.len(),
        |ring: &Ring, k| {
            let n = *ring.normal.direction;
            let [u, v] = ring.normal.across;
            let offset = local(k) - ring.center;
            let height = offset.dot(n);
            let radial = offset.component_perpendicular_to(ring.normal.direction);
            let distance = radial.magnitude();
            let outward = if distance > 0.0 {
                radial / distance
            } else {
                Vector::ZERO
            };
            let (along_u, along_v) = (offset.dot(u), offset.dot(v));
            let tilt = if distance > 0.0 {
                -height / distance
            } else {
                0.0
            };
            [
                (height, [-n.x, -n.y, -n.z, along_u, along_v, 0.0]),
                (
                    distance - ring.radius,
                    [
                        -outward.x,
                        -outward.y,
                        -outward.z,
                        tilt * along_u,
                        tilt * along_v,
                        -1.0,
                    ],
                ),
            ]
        },
        |ring, [dx, dy, dz, first, second, dr]| Ring {
            center: ring.center + Vector::new(dx, dy, dz),
            normal: ring.normal.turned(first, second),
            radius: ring.radius + dr,
        },
    );
    Some((
        spread.centroid + fitted.center * spread.scale,
        fitted.normal.direction,
        fitted.radius.abs() * spread.scale,
    ))
}

/// A cylinder while it is fitted: the point of its axis nearest the
/// centroid (the origin, once the points are moved to it), its unit axis
/// direction and its radius.
#[derive(Debug, Clone, Copy)]
struct Tube {
    origin: Vector,
    axis: TurningDirection,
    radius: f64,
}

/// The rings of the hemisphere of directions searched for a cylinder's
/// axis, from its pole to its equator: one every 15 degrees.
const DIRECTION_RINGS: usize = 6;

/// How many of the searched directions are refined, the best scored
/// ones that lie [`SEPARATION`] apart.
const REFINED_DIRECTIONS: usize = 3;

/// The least angle between two refined directions, in radians: 30
/// degrees, two steps of the search.
const SEPARATION: f64 = std::f64::consts::FRAC_PI_6;

/// The directions searched for a cylinder's axis: a hemisphere around the
/// axis of least spread, `axes[0]`, one direction every 15 degrees or so,
/// each with two unit vectors across it. The hemisphere holds each axis
/// line once, the three axes of `axes` among them, and maps onto itself
/// when an axis of `axes` turns the other way, so the search does not
/// depend on the signs the eigenvectors come out with.
fn searched_directions(axes: [Vector; 3]) -> impl Iterator<Item = TurningDirection> {
    let [pole, first, second] = axes;
    (0..=DIRECTION_RINGS).flat_map(move |ring| {
        let polar = ring as f64 * std::f64::consts::FRAC_PI_2 / DIRECTION_RINGS as f64;
        let (sine, cosine) = polar.sin_cos();
        // A multiple of four, so a turn of a quarter or a half maps the
        // ring onto itself; at the equator, half the ring, the other half
        // being the same lines.
        let around = 4 * ((DIRECTION_RINGS as f64 * sine).round() as usize).max(1);
        let count = if ring == 0 {
            1
        } else if ring == DIRECTION_RINGS {
            around / 2
        } else {
            around
        };
        (0..count).map(move |k| {
            let azimuth = k as f64 * std::f64::consts::TAU / around as f64;
            let (azimuth_sine, azimuth_cosine) = azimuth.sin_cos();
            let level = first * azimuth_cosine + second * azimuth_sine;
            TurningDirection {
                direction: Direction::from_unit_vector(pole * cosine + level * sine),
                across: [
                    level * cosine - pole * sine,
                    second * azimuth_cosine - first * azimuth_sine,
                ],
            }
        })
    })
}

/// The geometric least-squares cylinder of `points`, at least five and
/// not coplanar, spread as `spread`: the point of its axis nearest the
/// centroid, its unit axis direction (either way) and its radius.
pub(crate) fn cylinder(points: &[Point], spread: &Spread) -> Option<(Point, Direction, f64)> {
    let local = |k: usize| (points[k] - spread.centroid) / spread.scale;
    let mut candidates: Vec<(f64, Tube)> = searched_directions(spread.axes)
        .filter_map(|axis| {
            let [u, v] = axis.across;
            let (center, radius, score) = algebraic_circle((0..points.len()).map(|k| {
                let q = local(k);
                (q.dot(u), q.dot(v))
            }))?;
            let origin = u * center[0] + v * center[1];
            Some((
                score,
                Tube {
                    origin,
                    axis,
                    radius,
                },
            ))
        })
        .collect();
    candidates.sort_by(|a, b| a.0.total_cmp(&b.0));
    // The best scored directions, each at least `SEPARATION` from those
    // before it: few points on a long stretch of cylinder may score a
    // wrong direction best, whose refinement then settles in a local
    // minimum, while the right one scores next.
    let mut starts: Vec<Tube> = Vec::with_capacity(REFINED_DIRECTIONS);
    for (_, tube) in candidates {
        let apart =
            |other: &Tube| tube.axis.direction.dot(*other.axis.direction).abs() < SEPARATION.cos();
        if starts.iter().all(apart) {
            starts.push(tube);
            if starts.len() == REFINED_DIRECTIONS {
                break;
            }
        }
    }
    let residual = |tube: &Tube, k: usize| {
        let [u, v] = tube.axis.across;
        let offset = local(k) - tube.origin;
        let height = offset.dot(*tube.axis.direction);
        let radial = offset.component_perpendicular_to(tube.axis.direction);
        let distance = radial.magnitude();
        let (outward, tilt) = if distance > 0.0 {
            (radial / distance, -height / distance)
        } else {
            (Vector::ZERO, 0.0)
        };
        [(
            distance - tube.radius,
            [
                -outward.dot(u),
                -outward.dot(v),
                tilt * offset.dot(u),
                tilt * offset.dot(v),
                -1.0,
            ],
        )]
    };
    let advance = |tube: &Tube, [along_u, along_v, first, second, dr]: [f64; 5]| {
        let [u, v] = tube.axis.across;
        let axis = tube.axis.turned(first, second);
        let origin = tube.origin + u * along_u + v * along_v;
        Tube {
            // The axis point nearest the centroid, at the origin.
            origin: origin.component_perpendicular_to(axis.direction),
            axis,
            radius: tube.radius + dr,
        }
    };
    let squares = |tube: &Tube| -> f64 {
        (0..points.len())
            .map(|k| residual(tube, k)[0].0.powi(2))
            .sum()
    };
    let fitted = starts
        .into_iter()
        .map(|start| {
            let fitted = levenberg_marquardt(start, points.len(), residual, advance);
            (squares(&fitted), fitted)
        })
        .min_by(|a, b| a.0.total_cmp(&b.0))?
        .1;
    Some((
        spread.centroid + fitted.origin * spread.scale,
        fitted.axis.direction,
        fitted.radius.abs() * spread.scale,
    ))
}
