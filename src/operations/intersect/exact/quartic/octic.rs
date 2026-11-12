//! A torus against a cylinder or a cone in any position, along the
//! torus's tube circles.
//!
//! The torus's point `P(u, v) = C + (R + r cos v) e(u) + r sin v a` is of
//! degree one in the cosine and sine of each angle, and the cylinder's or
//! the cone's equation `G(P) = (P − Q)ᵀ M (P − Q) + c` is quadratic, so
//! `F(u, v) = G(P(u, v))` is a trigonometric polynomial of degree two in
//! each angle. On the tube circle at `u` it is
//! `A cos 2v + B sin 2v + C cos v + D sin v + E`, a quartic in
//! `tan(v / 2)`: up to four points per tube circle, each exact.
//!
//! The curve is of degree eight, and its points on the tube circles are
//! not one signed root as on the quadrics' generators: they appear and
//! vanish in pairs where a tube circle grazes the other surface, at the
//! points where `F = ∂F/∂v = 0`, and two branches cross where the
//! surfaces touch, at the points where `∂F/∂u = ∂F/∂v = 0` with `F = 0`.
//! Both are common zeros of two trigonometric polynomials, found by
//! subdivision (see [`trigonometric`](crate::math::numeric::trigonometric)).
//! Between the angles of those events the tube circles hold the same
//! number of points, in the same order around them from an angle no
//! branch crosses, and the events say how the branches join: a pair ends
//! where a tube circle grazes, two pairs swap where the surfaces touch and
//! the curve crosses itself, and the rest run on. Following the joins
//! gives the loops, each fitted through exact points as the quadrics'
//! quartics are.
//!
//! The surfaces touch along a curve only coaxially: a contact curve's
//! normals meet both the torus's centre circle and the other surface's
//! axis, which for a line and a circle happens along a curve only when
//! the line is the circle's axis. That case and the circles it gives are
//! read before.

use super::tube_circles::Tubes;
use super::*;
use crate::math::numeric::solve_trig;
use crate::math::numeric::trigonometric::{
    TrigonometricPolynomial, angular_distance, basis, common_zeros, wrap,
};
use crate::{Transform, Transformable};

/// `F(u, v)`, the other surface's equation on the torus's parameters.
type Equation = TrigonometricPolynomial<5>;

/// The products of the basis `(1, cos t, sin t)` in the basis
/// `(1, cos t, sin t, cos 2t, sin 2t)`.
const PRODUCTS: [[[f64; 5]; 3]; 3] = [
    [
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
    ],
    [
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.5, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.5],
    ],
    [
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.5],
        [0.5, 0.0, 0.0, -0.5, 0.0],
    ],
];

/// How far apart, in angle, the events of one tube circle may be: closer
/// ones are one event, as where the curve is symmetric.
const SAME_EVENT: f64 = 1e-9;

/// How far, in angle, a point where the curve turns back may be from a
/// point where the surfaces touch and still be read as part of it: the
/// rounding of a crossing splits it into two turns that close.
const NEAR_A_CROSSING: f64 = 1e-6;

/// How far, in angle, a point of a tube circle that grazes is kept from
/// the others: those nearer it on its own tube circle are its own.
const NEAR_A_GRAZE: f64 = 1e-5;

/// Where in the gap between two points on a tube circle the angle the
/// branches are ordered from lies: off its middle, where a symmetric
/// curve turns back.
const CUT_IN_A_GAP: f64 = 0.381_966_011_250_105_2;

/// How many times a stretch of tube circles is split where its branches
/// cross the angle they are ordered from.
const DEEPEST_CUT: u32 = 8;

/// The other surface: its equation's terms, and its own surface.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Quadric {
    Cylinder(Cylinder),
    Cone(Cone),
}

impl Quadric {
    /// `Q`, `M` and `c` of `(P − Q)ᵀ M (P − Q) + c`: a cylinder's squared
    /// distance to its axis less its radius squared, and a cone's
    /// `((P − V)·a)² − cos² φ |P − V|²`, both nappes.
    fn terms(&self) -> (Point, [[f64; 3]; 3], f64) {
        let (centre, axis, along, across, constant) = match self {
            Quadric::Cylinder(cylinder) => {
                let frame = cylinder.frame();
                let radius = cylinder.radius();
                (
                    frame.origin(),
                    frame.z_direction(),
                    -1.0,
                    1.0,
                    -radius * radius,
                )
            }
            Quadric::Cone(cone) => {
                let cos = cone.half_angle().cos();
                (
                    cone.apex(),
                    cone.frame().z_direction(),
                    1.0,
                    -cos * cos,
                    0.0,
                )
            }
        };
        let a = [axis.x, axis.y, axis.z];
        let matrix = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let identity = if i == j { across } else { 0.0 };
                identity + along * a[i] * a[j]
            })
        });
        (centre, matrix, constant)
    }

    /// The surface moved by `offset`.
    fn moved(&self, offset: Vector) -> Quadric {
        let transform = Transform::translation(offset);
        match self {
            Quadric::Cylinder(cylinder) => Quadric::Cylinder(
                cylinder
                    .transformed(transform)
                    .expect("a translation keeps a cylinder"),
            ),
            Quadric::Cone(cone) => Quadric::Cone(
                cone.transformed(transform)
                    .expect("a translation keeps a cone"),
            ),
        }
    }

    /// The signed distance from `point`, and the unit normal along which
    /// it grows.
    fn offset(&self, point: Point) -> (f64, Vector) {
        match self {
            Quadric::Cylinder(cylinder) => OtherCylinder::new(cylinder, Vector::ZERO).offset(point),
            Quadric::Cone(cone) => RuledCone::new(cone).offset(point),
        }
    }
}

/// `F(u, v)` for the torus's tube circles against `quadric`.
fn equation(tubes: &Tubes, quadric: &Quadric) -> Equation {
    let (centre, matrix, constant) = quadric.terms();
    let (major, minor) = (tubes.major, tubes.minor);
    // `P − Q` in the basis `(1, cos u, sin u) × (1, cos v, sin v)`.
    let mut terms = [[Vector::ZERO; 3]; 3];
    terms[0][0] = tubes.origin - centre;
    terms[1][0] = tubes.x * major;
    terms[2][0] = tubes.y * major;
    terms[1][1] = tubes.x * minor;
    terms[2][1] = tubes.y * minor;
    terms[0][2] = tubes.axis * minor;
    let apply = |w: Vector| {
        let w = [w.x, w.y, w.z];
        Vector::new(
            (0..3).map(|j| matrix[0][j] * w[j]).sum(),
            (0..3).map(|j| matrix[1][j] * w[j]).sum(),
            (0..3).map(|j| matrix[2][j] * w[j]).sum(),
        )
    };
    let mut c = [[0.0; 5]; 5];
    c[0][0] = constant;
    for (i1, row1) in terms.iter().enumerate() {
        for (j1, first) in row1.iter().enumerate() {
            if *first == Vector::ZERO {
                continue;
            }
            let mapped = apply(*first);
            for (i2, row2) in terms.iter().enumerate() {
                for (j2, second) in row2.iter().enumerate() {
                    let product = mapped.dot(*second);
                    if product == 0.0 {
                        continue;
                    }
                    let (along_u, along_v) = (PRODUCTS[i1][i2], PRODUCTS[j1][j2]);
                    for (a, &x) in along_u.iter().enumerate() {
                        if x == 0.0 {
                            continue;
                        }
                        for (b, &y) in along_v.iter().enumerate() {
                            c[a][b] += product * x * y;
                        }
                    }
                }
            }
        }
    }
    TrigonometricPolynomial { c }
}

/// The angles `v` where the tube circle at `u` meets the other surface,
/// in `(−π, π]`, a double one twice, but those within rounding of a point
/// where the surfaces touch and part, `alone`: there `F` is within
/// rounding of zero, and a tube circle through it may hold two points of
/// noise near it.
fn points_at(equation: &Equation, u: f64, alone: &[[f64; 2]]) -> InlineList<f64, 4> {
    let mut out = InlineList::new();
    for &v in &roots_at(equation, u) {
        if alone
            .iter()
            .all(|&point| angular_distance(point, [u, v]) > NEAR_A_CROSSING)
        {
            out.push(v);
        }
    }
    out
}

/// The angles `v` where the tube circle at `u` meets the other surface,
/// in `(−π, π]`, a double one twice.
fn roots_at(equation: &Equation, u: f64) -> InlineList<f64, 4> {
    let [e, c1, d1, a2, b2] = equation.along_v(u);
    let exact = Tolerance {
        confusion: 0.0,
        ..Tolerance::DEFAULT
    };
    let mut out = InlineList::new();
    for (v, multiplicity) in &solve_trig(a2, b2, c1, d1, e, exact, |_| false) {
        for _ in 0..*multiplicity {
            if out.len() < 4 {
                out.push(*v);
            }
        }
    }
    out
}

/// The angles `u` where the circle of the torus at `v` meets the other
/// surface.
fn roots_along(equation: &Equation, v: f64) -> InlineList<(f64, u32), 4> {
    let [e, c1, d1, a2, b2] = equation.along_u(v);
    let exact = Tolerance {
        confusion: 0.0,
        ..Tolerance::DEFAULT
    };
    solve_trig(a2, b2, c1, d1, e, exact, |_| false)
}

/// `values` sorted around the circle from `cut`.
fn ordered_from(values: &[f64], cut: f64) -> InlineList<f64, 4> {
    let mut out: InlineList<f64, 4> = InlineList::new();
    for &v in values {
        out.push(v);
    }
    out.sort_unstable_by(|a, b| {
        (a - cut)
            .rem_euclid(TAU)
            .total_cmp(&(b - cut).rem_euclid(TAU))
    });
    out
}

/// A stretch `[u₀, u₁]` of tube circles, each holding `count` points of
/// the curve, ordered around it from `cut`, which none crosses there.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Piece {
    u: [f64; 2],
    cut: f64,
    count: u8,
    /// Where its branches start in [`OcticCurve::branches`].
    first_branch: u32,
}

/// A branch of a [`Piece`]: its `v` at each end, and whether it turns back
/// there, joining the next branch on a tube circle that grazes.
#[derive(Debug, Clone, Copy, PartialEq)]
struct BranchEnds {
    v: [f64; 2],
    turning: [bool; 2],
}

/// A branch of a piece on a loop, run along `u` or back.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Visit {
    piece: u32,
    branch: u8,
    forward: bool,
}

/// A loop: its visits in [`OcticCurve::visits`], and whether it runs
/// through a point where the surfaces touch, where the normals turn
/// parallel.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Loop {
    first_visit: u32,
    visits: u32,
    crosses: bool,
}

/// The curve of degree eight where a torus meets a cylinder or a cone away
/// from their common axis: loops read along the torus's tube circles,
/// each a smooth piece of the curve, and the points where the surfaces
/// touch, at a point no loop reaches or where a loop crosses itself.
///
/// A loop is evaluated at a parameter in `[0, 2π]`, its two ends the same
/// point. The routed [`intersect`](crate::Intersect) fits each as a
/// [`BSplineCurve`](crate::BSplineCurve) through exact points.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Cylinder, Point, Tolerance, Torus, TorusCylinderIntersection, Vector};
/// // A pipe of radius 0.5 dropped through the tube.
/// let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
/// let pipe = Cylinder::new(Point::new(4.0, 0.0, 0.0), Vector::Z, 0.5).unwrap();
/// let TorusCylinderIntersection::Octic(curve) = torus.intersect_exact(&pipe, Tolerance::DEFAULT)
/// else {
///     panic!("expected the octic")
/// };
/// assert_eq!(curve.component_count(), 2);
/// let p = curve.point_at(0, 1.0);
/// assert!(((p.x - 4.0).hypot(p.y) - 0.5).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct OcticCurve {
    torus: Torus,
    /// The other surface as the curve is read on: moved by at most the
    /// confusion tolerance where it touches the torus within it.
    other: Quadric,
    /// Whether the torus is the first operand, its curve in parameters
    /// the first.
    torus_first: bool,
    equation: Equation,
    pieces: Vec<Piece>,
    branches: Vec<BranchEnds>,
    visits: Vec<Visit>,
    loops: Vec<Loop>,
    /// The points where the surfaces touch.
    touching: Vec<Point>,
    /// The angles of those where they touch and part.
    alone: Vec<[f64; 2]>,
    /// Those where the curve crosses itself.
    crossings: Vec<Crossing>,
}

/// A point where the surfaces touch and two branches cross: its angles,
/// and the slopes `dv/du` of the branches through it.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Crossing {
    angles: [f64; 2],
    slopes: [f64; 2],
}

/// How far in `u` from a crossing its branches are read from their slopes
/// there and refined: nearer, the quartic on the tube circle merges their
/// two points into a double one.
const NEAR_A_CROSSING_IN_U: f64 = 1e-4;

/// How far in `u` from a crossing a point of a branch is the crossing
/// itself: nearer, `F`'s rounding moves a point along the tube circle by
/// more than the branch's distance from it.
const AT_A_CROSSING: f64 = 1e-8;

impl OcticCurve {
    /// The number of loops.
    pub fn component_count(&self) -> usize {
        self.loops.len()
    }

    /// The points where the surfaces touch: where a loop crosses itself
    /// or two loops cross, and where they touch at a point alone.
    pub fn contact_points(&self) -> &[Point] {
        &self.touching
    }

    /// The point of loop `component` at `parameter` in `[0, 2π]`.
    ///
    /// # Panics
    ///
    /// When `component` is not below [`component_count`](Self::component_count).
    pub fn point_at(&self, component: usize, parameter: f64) -> Point {
        let tubes = Tubes::new(&self.torus);
        let [u, v] = self.angles_at(component, parameter);
        tubes.point(Point2D::new(u, v))
    }

    /// The same curve with the operands the other way round.
    pub(crate) fn operands_swapped(mut self: Box<Self>) -> Box<OcticCurve> {
        self.torus_first = !self.torus_first;
        self
    }

    /// Each loop fitted as cubic Hermite segments through exact points and
    /// tangents, within a quarter of `tol.confusion` of the curve, with its
    /// curves in the parameters of both; `None` when a fit fails.
    pub(crate) fn fitted(&self, tol: Tolerance) -> Option<Vec<IntersectionCurve>> {
        let components: Vec<usize> = (0..self.loops.len()).collect();
        let tubes = Tubes::new(&self.torus);
        match &self.other {
            Quadric::Cylinder(cylinder) => {
                let other = OtherCylinder::new(cylinder, Vector::ZERO);
                OcticFamily::new(self, tubes, other).fitted(
                    &components,
                    tol.confusion,
                    self.torus_first,
                )
            }
            Quadric::Cone(cone) => {
                let other = OtherCone::new(cone);
                OcticFamily::new(self, tubes, other).fitted(
                    &components,
                    tol.confusion,
                    self.torus_first,
                )
            }
        }
    }

    /// The angles `(u, v)` on the torus of the point of loop `component`
    /// at `parameter`.
    fn angles_at(&self, component: usize, parameter: f64) -> [f64; 2] {
        let lap = self.loops[component];
        let count = lap.visits as usize;
        let scaled = (parameter / TAU).clamp(0.0, 1.0) * count as f64;
        let k = (scaled.floor() as usize).min(count - 1);
        let visit = self.visits[lap.first_visit as usize + k];
        let s = scaled - k as f64;
        let s = if visit.forward { s } else { 1.0 - s };
        let piece = self.pieces[visit.piece as usize];
        let ends = self.branches[piece.first_branch as usize + visit.branch as usize];
        // Where a branch turns back, `v` goes as the square root of the
        // distance in `u`: `u` then goes as the square of the parameter.
        let h = match ends.turning {
            [true, true] => (0.5 * PI * s).sin().powi(2),
            [true, false] => 1.0 - (0.5 * PI * s).cos(),
            [false, true] => (0.5 * PI * s).sin(),
            [false, false] => s,
        };
        if h <= 0.0 {
            return [piece.u[0], ends.v[0]];
        }
        if h >= 1.0 {
            return [piece.u[1], ends.v[1]];
        }
        let u = piece.u[0] + (piece.u[1] - piece.u[0]) * h;
        let mut points = points_at(&self.equation, u, &self.alone);
        for crossing in &self.crossings {
            let [at_u, at_v] = crossing.angles;
            let from = (u - at_u + PI).rem_euclid(TAU) - PI;
            if from.abs() > NEAR_A_CROSSING_IN_U {
                continue;
            }
            // The two points nearest the crossing, from the branches.
            let steepest = crossing.slopes[0].abs().max(crossing.slopes[1].abs());
            let window = 2.0 * steepest * from.abs() + 1e-6;
            let mut kept = InlineList::new();
            for &v in &points {
                if angular_distance([0.0, v], [0.0, at_v]) > window {
                    kept.push(v);
                }
            }
            for slope in crossing.slopes {
                if kept.len() < 4 {
                    kept.push(if from.abs() <= AT_A_CROSSING {
                        at_v
                    } else {
                        self.on_tube_circle(u, at_v + slope * from)
                    });
                }
            }
            points = kept;
        }
        let roots = ordered_from(&points, piece.cut);
        if roots.len() == usize::from(piece.count) {
            [u, roots[visit.branch as usize]]
        } else {
            // Within rounding of an end, where two points merge.
            let end = usize::from(h >= 0.5);
            [piece.u[end], ends.v[end]]
        }
    }

    /// The point of the tube circle at `u` on the curve nearest `v`, by
    /// Newton's method along it.
    fn on_tube_circle(&self, u: f64, mut v: f64) -> f64 {
        let along = self.equation.along_v(u);
        let value = |coefficients: &[f64; 5], v: f64| {
            let at = basis::<5>(v);
            coefficients
                .iter()
                .zip(&at)
                .map(|(c, x)| c * x)
                .sum::<f64>()
        };
        // `(a cos kv + b sin kv)' = k b cos kv − k a sin kv`.
        let rate = [0.0, along[2], -along[1], 2.0 * along[4], -2.0 * along[3]];
        for _ in 0..8 {
            let slope = value(&rate, v);
            if slope == 0.0 {
                break;
            }
            let step = value(&along, v) / slope;
            v -= step;
            if step.abs() <= 1e-15 {
                break;
            }
        }
        v
    }
}

/// What a cylinder or a cone off the torus's axis cuts from it.
pub(in crate::operations::intersect::exact) enum OcticCut {
    Curve(Box<OcticCurve>),
    TangentPoint(Point),
    Empty,
}

/// A torus against a cylinder off its axis.
pub(in crate::operations::intersect::exact) fn torus_cylinder(
    torus: &Torus,
    cylinder: &Cylinder,
    tol: Tolerance,
) -> Option<OcticCut> {
    torus_quadric(torus, Quadric::Cylinder(*cylinder), tol)
}

/// A torus against a cone off its axis.
pub(in crate::operations::intersect::exact) fn torus_cone(
    torus: &Torus,
    cone: &Cone,
    tol: Tolerance,
) -> Option<OcticCut> {
    torus_quadric(torus, Quadric::Cone(*cone), tol)
}

/// A point where the surfaces touch: a crossing of two branches, or a
/// point alone.
struct Touching {
    angles: [f64; 2],
    crossing: bool,
    /// The slopes `dv/du` of the branches through a crossing.
    slopes: [f64; 2],
}

/// A point where a branch pair ends on a tube circle that grazes: on the
/// side of smaller `u` (`opens_back`) or of larger.
struct Turn {
    angles: [f64; 2],
    opens_back: bool,
}

/// The curve of a ring torus against `quadric`, `None` where it is not
/// read in closed form: a horn or spindle torus, a tube circle lying on
/// the other surface, surfaces touching to a higher order than a crossing
/// or a point, or touching at several points that no single move within
/// tolerance makes exact.
fn torus_quadric(torus: &Torus, mut quadric: Quadric, tol: Tolerance) -> Option<OcticCut> {
    let tubes = Tubes::of_ring(torus, tol)?;
    // A cone's apex within confusion of the torus is moved onto it.
    if let Quadric::Cone(cone) = quadric {
        let (gap, normal) = tubes.offset(cone.apex());
        if gap.abs() <= tol.confusion && gap != 0.0 {
            quadric = quadric.moved(normal * -gap);
        }
    }
    let reach = (tubes.origin - quadric.terms().0).magnitude() + tubes.major + tubes.minor;
    // `|F|` where the surfaces are within confusion: `|∇G| ≤ 2 |P − Q|`.
    let within = 4.0 * reach * tol.confusion;
    let mut touching = touching_points(&tubes, &quadric, within, tol)?;
    if let Some(snapped) = touching_exactly(&tubes, &quadric, &touching, tol) {
        quadric = snapped?;
        touching = touching_points(&tubes, &quadric, within, tol)?;
    }
    let equation = equation(&tubes, &quadric);
    let classified = classify(&equation, &touching)?;
    let turns = turns(&equation, &classified)?;
    let alone: Vec<[f64; 2]> = classified
        .iter()
        .filter(|t| !t.crossing)
        .map(|t| t.angles)
        .collect();
    let curve = sweep(torus, &tubes, quadric, equation, &classified, &turns, alone)?;
    Some(match (curve.loops.is_empty(), curve.touching.as_slice()) {
        (true, []) => OcticCut::Empty,
        (true, [point]) => OcticCut::TangentPoint(*point),
        _ => OcticCut::Curve(Box::new(curve)),
    })
}

/// The other surface moved to touch the torus exactly where it touches it
/// within confusion: by the one translation that closes every gap, or for
/// a cylinder by the one change of radius. `None` when nothing needs
/// moving, `Some(None)` when no single move closes every gap.
fn touching_exactly(
    tubes: &Tubes,
    quadric: &Quadric,
    touching: &[[f64; 2]],
    tol: Tolerance,
) -> Option<Option<Quadric>> {
    let gaps: Vec<(f64, Vector)> = touching
        .iter()
        .map(|&[u, v]| quadric.offset(tubes.point(Point2D::new(u, v))))
        .collect();
    if gaps.iter().all(|(gap, _)| *gap == 0.0) {
        return None;
    }
    let agree = SNAPS_AGREE * tol.confusion;
    let (first_gap, first_normal) = gaps[0];
    let first_move = first_normal * first_gap;
    if gaps
        .iter()
        .all(|&(gap, normal)| (normal * gap - first_move).magnitude() <= agree)
    {
        return Some(Some(quadric.moved(first_move)));
    }
    if let Quadric::Cylinder(cylinder) = quadric
        && gaps.iter().all(|(gap, _)| (gap - first_gap).abs() <= agree)
    {
        let widened = Cylinder::from_axis(cylinder.axis(), cylinder.radius() + first_gap).ok();
        return Some(widened.map(Quadric::Cylinder));
    }
    Some(None)
}

/// How near, in units of confusion, the moves that close each gap where
/// the surfaces touch must be to be one move: farther, the points would
/// stay that far from touching, read as touching.
const SNAPS_AGREE: f64 = 1e-3;

/// The angles of the points where the surfaces are within confusion of
/// touching: where `F` is stationary and the torus's point within
/// confusion of the other surface. `None` where they touch along more
/// than points.
fn touching_points(
    tubes: &Tubes,
    quadric: &Quadric,
    within: f64,
    tol: Tolerance,
) -> Option<Vec<[f64; 2]>> {
    let equation = equation(tubes, quadric);
    let (along_u, along_v) = (equation.derivative_u(), equation.derivative_v());
    let zeros = common_zeros(&along_u, &along_v, Some((&equation, within)))?;
    let near = |&[u, v]: &[f64; 2]| {
        quadric.offset(tubes.point(Point2D::new(u, v))).0.abs() <= tol.confusion
    };
    if zeros.clusters.iter().any(near) {
        return None;
    }
    Some(zeros.simple.into_iter().filter(near).collect())
}

/// The points where the surfaces touch, each a crossing or a point alone
/// by the sign of `F`'s Hessian there. `None` where it is singular, the
/// surfaces touching to a higher order, or where a branch through a
/// crossing runs along a tube circle.
fn classify(equation: &Equation, touching: &[[f64; 2]]) -> Option<Vec<Touching>> {
    let (along_u, along_v) = (equation.derivative_u(), equation.derivative_v());
    let (uu, uv, vv) = (
        along_u.derivative_u(),
        along_u.derivative_v(),
        along_v.derivative_v(),
    );
    touching
        .iter()
        .map(|&[u, v]| {
            let (a, b, c) = (uu.value(u, v), uv.value(u, v), vv.value(u, v));
            let size = a.abs().max(b.abs()).max(c.abs());
            let determinant = a * c - b * b;
            if determinant.abs() <= 1e-9 * size * size {
                return None;
            }
            let crossing = determinant < 0.0;
            if crossing && c.abs() <= 1e-6 * size {
                return None;
            }
            // `a + 2b m + c m² = 0` along each branch.
            let root = (-determinant).max(0.0).sqrt();
            Some(Touching {
                angles: [u, v],
                crossing,
                slopes: [(-b - root) / c, (-b + root) / c],
            })
        })
        .collect()
}

/// The points where a tube circle grazes the other surface and a pair of
/// branches ends, away from where the surfaces touch. `None` where one is
/// not simple: the curve turning to a higher order.
fn turns(equation: &Equation, touching: &[Touching]) -> Option<Vec<Turn>> {
    let along_v = equation.derivative_v();
    let zeros = common_zeros(equation, &along_v, None)?;
    let near_touching = |angles: [f64; 2]| {
        touching
            .iter()
            .any(|t| angular_distance(t.angles, angles) <= NEAR_A_CROSSING)
    };
    if zeros.clusters.iter().any(|&angles| !near_touching(angles)) {
        return None;
    }
    let (along_u, second_v) = (equation.derivative_u(), along_v.derivative_v());
    zeros
        .simple
        .into_iter()
        .filter(|&angles| !near_touching(angles))
        .map(|[u, v]| {
            let (slope, bend) = (along_u.value(u, v), second_v.value(u, v));
            if bend.abs() <= 1e-9 * slope.abs() {
                return None;
            }
            // `F ≈ F_u δu + F_vv δv² / 2` vanishes for `δu` of the sign
            // of `−F_u F_vv`.
            Some(Turn {
                angles: [u, v],
                opens_back: slope * bend > 0.0,
            })
        })
        .collect()
}

/// What happens to the branches at a point of an event's tube circle: how
/// many arrive from smaller `u` and leave to larger.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Slot {
    /// A branch runs through.
    Through,
    /// Two branches from smaller `u` join.
    ClosesBack,
    /// Two branches to larger `u` start.
    OpensAhead,
    /// Two branches cross where the surfaces touch.
    Crossing,
}

impl Slot {
    /// The branches on the side of smaller `u`, then of larger.
    fn counts(self) -> [usize; 2] {
        match self {
            Slot::Through => [1, 1],
            Slot::ClosesBack => [2, 0],
            Slot::OpensAhead => [0, 2],
            Slot::Crossing => [2, 2],
        }
    }
}

/// An end of a branch: `(piece, branch, end)`, `end` `0` at its smaller
/// `u`.
type End = (usize, usize, usize);

/// Splits the tube circles into pieces between the events, and follows
/// the branches through them into loops.
fn sweep(
    torus: &Torus,
    tubes: &Tubes,
    quadric: Quadric,
    equation: Equation,
    touching: &[Touching],
    turns: &[Turn],
    alone: Vec<[f64; 2]>,
) -> Option<OcticCurve> {
    // The events: each a tube circle and what happens on it.
    let mut critical: Vec<(f64, f64, Slot)> = turns
        .iter()
        .map(|turn| {
            let slot = if turn.opens_back {
                Slot::ClosesBack
            } else {
                Slot::OpensAhead
            };
            (wrap(turn.angles[0]), turn.angles[1], slot)
        })
        .chain(
            touching
                .iter()
                .filter(|t| t.crossing)
                .map(|t| (wrap(t.angles[0]), t.angles[1], Slot::Crossing)),
        )
        .collect();
    critical.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut events: Vec<(f64, Vec<(f64, Slot)>)> = Vec::new();
    for (u, v, slot) in critical {
        match events.last_mut() {
            Some((at, slots)) if u - *at <= SAME_EVENT => slots.push((v, slot)),
            _ => events.push((u, vec![(v, slot)])),
        }
    }
    if events.len() > 1 && events[0].0 + TAU - events[events.len() - 1].0 <= SAME_EVENT {
        let (_, slots) = events.pop().expect("more than one event");
        events[0].1.extend(slots);
    }
    if events.is_empty() {
        events.push((0.0, Vec::new()));
    }
    // The pieces between consecutive events, split where needed so that
    // no branch crosses the angle it is ordered from; each boundary is the
    // tube circle of an event, or one where nothing happens.
    let mut pieces: Vec<Piece> = Vec::new();
    let mut boundaries: Vec<(f64, Vec<(f64, Slot)>)> = Vec::new();
    for (k, (u, slots)) in events.iter().enumerate() {
        let next = if k + 1 < events.len() {
            events[k + 1].0
        } else {
            events[0].0 + TAU
        };
        boundaries.push((*u, slots.clone()));
        let before = pieces.len();
        split(&equation, *u, next, &alone, 0, &mut pieces)?;
        for piece in &pieces[before + 1..] {
            boundaries.push((piece.u[0], Vec::new()));
        }
    }
    let mut branches = Vec::new();
    for piece in &mut pieces {
        piece.first_branch = branches.len() as u32;
        branches.extend(std::iter::repeat_n(
            BranchEnds {
                v: [0.0; 2],
                turning: [false; 2],
            },
            usize::from(piece.count),
        ));
    }
    // The joins at each boundary, between the piece before and after.
    let mut joined: Vec<Option<End>> = vec![None; 2 * branches.len()];
    let slot_of =
        |(piece, branch, end): End| 2 * (pieces[piece].first_branch as usize + branch) + end;
    for (k, (u, critical)) in boundaries.iter().enumerate() {
        let (back, ahead) = ((k + pieces.len() - 1) % pieces.len(), k);
        let mut slots: Vec<(f64, Slot)> = critical.clone();
        for v in &points_at(&equation, *u, &alone) {
            if critical
                .iter()
                .all(|(w, _)| angular_distance([0.0, *v], [0.0, *w]) > NEAR_A_GRAZE)
            {
                slots.push((*v, Slot::Through));
            }
        }
        // The branches on each side, in each piece's order.
        let sides: [Vec<(f64, Slot, usize)>; 2] = std::array::from_fn(|side| {
            let piece = pieces[if side == 0 { back } else { ahead }];
            let mut ordered: Vec<(f64, Slot, usize)> = slots
                .iter()
                .enumerate()
                .map(|(index, &(v, slot))| (v, slot, index))
                .collect();
            ordered.sort_by(|a, b| {
                (a.0 - piece.cut)
                    .rem_euclid(TAU)
                    .total_cmp(&(b.0 - piece.cut).rem_euclid(TAU))
            });
            ordered
                .into_iter()
                .flat_map(|entry| std::iter::repeat_n(entry, entry.1.counts()[side]))
                .collect()
        });
        if sides[0].len() != usize::from(pieces[back].count)
            || sides[1].len() != usize::from(pieces[ahead].count)
        {
            return None;
        }
        // The first branch of each slot on each side.
        let first_of = |side: usize, index: usize| sides[side].iter().position(|e| e.2 == index);
        for (index, &(v, slot)) in slots.iter().enumerate() {
            let back_branch = first_of(0, index);
            let ahead_branch = first_of(1, index);
            let mut end = |end: End, turning: bool| {
                let ends = &mut branches[pieces[end.0].first_branch as usize + end.1];
                ends.v[end.2] = v;
                ends.turning[end.2] = turning;
            };
            let mut join = |a: End, b: End| {
                joined[slot_of(a)] = Some(b);
                joined[slot_of(b)] = Some(a);
            };
            match slot {
                Slot::Through => {
                    let (b, a) = (back_branch?, ahead_branch?);
                    end((back, b, 1), false);
                    end((ahead, a, 0), false);
                    join((back, b, 1), (ahead, a, 0));
                }
                Slot::ClosesBack => {
                    let b = back_branch?;
                    end((back, b, 1), true);
                    end((back, b + 1, 1), true);
                    join((back, b, 1), (back, b + 1, 1));
                }
                Slot::OpensAhead => {
                    let a = ahead_branch?;
                    end((ahead, a, 0), true);
                    end((ahead, a + 1, 0), true);
                    join((ahead, a, 0), (ahead, a + 1, 0));
                }
                Slot::Crossing => {
                    let (b, a) = (back_branch?, ahead_branch?);
                    for k in 0..2 {
                        end((back, b + k, 1), false);
                        end((ahead, a + k, 0), false);
                    }
                    // The lower branch from before runs on as the upper one.
                    join((back, b, 1), (ahead, a + 1, 0));
                    join((back, b + 1, 1), (ahead, a, 0));
                }
            }
        }
    }
    // The loops, following the joins.
    let crossings: Vec<[f64; 2]> = touching
        .iter()
        .filter(|t| t.crossing)
        .map(|t| t.angles)
        .collect();
    let mut visited = vec![false; branches.len()];
    let mut visits = Vec::new();
    let mut loops = Vec::new();
    let branch_of = |piece: usize, branch: usize| pieces[piece].first_branch as usize + branch;
    for (piece_index, piece) in pieces.iter().enumerate() {
        for branch in 0..usize::from(piece.count) {
            if visited[branch_of(piece_index, branch)] {
                continue;
            }
            let first_visit = visits.len() as u32;
            let mut crosses = false;
            let start: End = (piece_index, branch, 0);
            let mut at = start;
            loop {
                let (p, b, entered) = at;
                visited[branch_of(p, b)] = true;
                visits.push(Visit {
                    piece: p as u32,
                    branch: b as u8,
                    forward: entered == 0,
                });
                let leaving = (p, b, 1 - entered);
                let ends = branches[branch_of(p, b)];
                let leaving_angles = [pieces[p].u[1 - entered], ends.v[1 - entered]];
                crosses |= crossings
                    .iter()
                    .any(|&c| angular_distance(c, leaving_angles) <= NEAR_A_CROSSING);
                at = joined[slot_of(leaving)]?;
                if at == start {
                    break;
                }
                if visited[branch_of(at.0, at.1)] || visits.len() > 4 * branches.len() {
                    return None;
                }
            }
            loops.push(Loop {
                first_visit,
                visits: visits.len() as u32 - first_visit,
                crosses,
            });
        }
    }
    let touching_points = touching
        .iter()
        .map(|t| tubes.point(Point2D::new(t.angles[0], t.angles[1])))
        .collect();
    Some(OcticCurve {
        torus: *torus,
        other: quadric,
        torus_first: true,
        equation,
        pieces,
        branches,
        visits,
        loops,
        touching: touching_points,
        alone,
        crossings: touching
            .iter()
            .filter(|t| t.crossing)
            .map(|t| Crossing {
                angles: t.angles,
                slopes: t.slopes,
            })
            .collect(),
    })
}

/// Appends the pieces from `u₀` to `u₁`, split where the branches cross
/// the angle they are ordered from. The angle is the first of a few that
/// no branch crosses there: one in each gap between the points on the
/// middle tube circle, then angles spread around it, off the ones a
/// symmetric curve turns back on. Where every one is crossed, the stretch
/// is split where the one crossed least is. `None` where no angle is
/// found that the branches do not cross at an end.
fn split(
    equation: &Equation,
    start: f64,
    end: f64,
    alone: &[[f64; 2]],
    depth: u32,
    out: &mut Vec<Piece>,
) -> Option<()> {
    let middle = 0.5 * (start + end);
    let mut roots: Vec<f64> = points_at(equation, middle, alone).iter().copied().collect();
    roots.sort_by(f64::total_cmp);
    let count = roots.len() as u8;
    if roots.is_empty() {
        out.push(Piece {
            u: [start, end],
            cut: 0.0,
            count,
            first_branch: 0,
        });
        return Some(());
    }
    let gaps = (0..roots.len()).map(|k| {
        let (from, to) = (
            roots[k],
            if k + 1 < roots.len() {
                roots[k + 1]
            } else {
                roots[0] + TAU
            },
        );
        from + CUT_IN_A_GAP * (to - from)
    });
    let spread = (0..SPREAD_CUTS).map(|k| (k as f64 + CUT_IN_A_GAP) * TAU / SPREAD_CUTS as f64);
    // The fewest crossings inside, of an angle not crossed at an end.
    let mut fewest: Option<Vec<f64>> = None;
    for cut in gaps.chain(spread) {
        let Some(inside) = crossings(equation, cut, start, end) else {
            continue;
        };
        if inside.is_empty() {
            out.push(Piece {
                u: [start, end],
                cut: wrap(cut),
                count,
                first_branch: 0,
            });
            return Some(());
        }
        if fewest.as_ref().is_none_or(|best| inside.len() < best.len()) {
            fewest = Some(inside);
        }
    }
    if depth >= DEEPEST_CUT {
        return None;
    }
    let mut from = start;
    for u in fewest?.into_iter().chain(std::iter::once(end)) {
        split(equation, from, u, alone, depth + 1, out)?;
        from = u;
    }
    Some(())
}

/// How many angles spread around the tube circle are tried to order its
/// points from.
const SPREAD_CUTS: usize = 16;

/// The angles `u` strictly inside `(start, end)` where a branch crosses
/// `v = cut`, ascending, `None` when one does within [`SAME_EVENT`] of an
/// end.
fn crossings(equation: &Equation, cut: f64, start: f64, end: f64) -> Option<Vec<f64>> {
    let mut inside = Vec::new();
    for &(u, _) in &roots_along(equation, cut) {
        // In `[start, start + 2π)`, the stretch's end at most there.
        let u = start + (u - start).rem_euclid(TAU);
        if u <= start + SAME_EVENT || u >= start + TAU - SAME_EVENT || (u - end).abs() <= SAME_EVENT
        {
            return None;
        }
        if u < end {
            inside.push(u);
        }
    }
    inside.sort_by(f64::total_cmp);
    Some(inside)
}

/// The torus's tube circles read against the other surface along the
/// curve's loops.
#[derive(Clone, Copy)]
struct OcticFamily<'a, S> {
    curve: &'a OcticCurve,
    tubes: Tubes,
    other: S,
    other_parameters_fit: bool,
}

/// A point of a loop: its angles on the torus, `ρ(u)`, the point, and the
/// cosine and sine of its `v`.
struct OnTube {
    u: f64,
    v: f64,
    radial: Vector,
    point: Point,
    cos_v: f64,
    sin_v: f64,
}

impl<'a, S: SecondSurface> OcticFamily<'a, S> {
    fn new(curve: &'a OcticCurve, tubes: Tubes, other: S) -> Self {
        OcticFamily {
            curve,
            tubes,
            other,
            other_parameters_fit: true,
        }
    }

    #[inline]
    fn at(&self, component: usize, parameter: f64) -> OnTube {
        let [u, v] = self.curve.angles_at(component, parameter);
        let ((sin_u, cos_u), (sin_v, cos_v)) = (u.sin_cos(), v.sin_cos());
        let radial = self.tubes.radial(sin_u, cos_u);
        OnTube {
            u,
            v,
            radial,
            point: self.tubes.point_on(radial, cos_v, sin_v),
            cos_v,
            sin_v,
        }
    }

    /// The unit tangent of the curve at `on`, along the cross product of
    /// the two surfaces' normals, before orientation.
    #[inline]
    fn crossing(&self, on: &OnTube) -> Vector {
        let torus_normal = on.radial * on.cos_v + self.tubes.axis * on.sin_v;
        let along = torus_normal.cross(self.other.normal(on.point));
        along * (1.0 / along.magnitude())
    }

    fn crosses(&self, component: usize) -> bool {
        self.curve.loops[component].crosses
    }
}

impl<S: SecondSurface> CurveFamily for OcticFamily<'_, S> {
    type Component = usize;

    #[inline]
    fn slides_along(&self, component: usize) -> bool {
        self.crosses(component)
    }

    #[inline]
    fn point(&self, component: usize, parameter: f64) -> Point {
        self.at(component, parameter).point
    }

    #[inline]
    fn point_and_tangent(&self, component: usize, parameter: f64) -> (Point, Vector) {
        let on = self.at(component, parameter);
        let tangent = if self.crosses(component) {
            self.stepped_tangent(component, parameter)
        } else {
            self.crossing(&on)
        };
        (on.point, tangent)
    }

    #[inline]
    fn sample(&self, component: usize, parameter: f64, orientation: f64) -> Sample {
        let on = self.at(component, parameter);
        let tangent = if self.crosses(component) {
            self.stepped_tangent(component, parameter)
        } else {
            self.crossing(&on) * orientation
        };
        let (other_uv, other_derivative) = self.other.uv(on.point, tangent);
        Sample {
            parameter,
            point: on.point,
            tangent,
            uv: [Point2D::new(on.u, on.v), other_uv],
            uv_derivatives: [
                self.tubes
                    .uv_derivative(tangent, on.radial, on.cos_v, on.sin_v),
                other_derivative,
            ],
        }
    }

    /// The way the loop runs at its start, from its exact points on either
    /// side: the tangent through a crossing follows its parameter.
    fn orientation(&self, component: usize) -> f64 {
        if self.crosses(component) {
            return 1.0;
        }
        let along = self.crossing(&self.at(component, 0.0));
        if along.dot(self.stepped_tangent(component, 0.0)) >= 0.0 {
            1.0
        } else {
            -1.0
        }
    }

    #[inline]
    fn offsets(&self, point: Point) -> [(f64, Vector); 2] {
        [self.tubes.offset(point), self.other.offset(point)]
    }

    #[inline]
    fn surface_point(&self, side: usize, uv: Point2D) -> Point {
        if side == 0 {
            self.tubes.point(uv)
        } else {
            self.other.point(uv)
        }
    }

    /// Both of the torus's parameters are angles; a cylinder's or a
    /// cone's first is, its second runs along the axis.
    #[inline]
    fn angles(&self, side: usize) -> [bool; 2] {
        if side == 0 {
            [true, true]
        } else {
            [true, false]
        }
    }

    #[inline]
    fn second_parameters_fit(&self) -> bool {
        self.other_parameters_fit
    }

    fn checked(mut self, confusion: f64) -> Self {
        let tubes = self.tubes;
        self.other_parameters_fit = self
            .other
            .parameters_fit(|point| tubes.offset(point).0, confusion);
        self
    }
}
