//! The local minima of the distance between a curve and a surface
//! ([`CurveSurfaceExtrema`]).

use std::f64::consts::TAU;
use std::ops::Deref;

use crate::geometry::curves::{ElementaryCurve, ParametricCurve};
use crate::geometry::surfaces::{ElementarySurface, ParametricSurface};
use crate::math::numeric::newton_min_1d;
use crate::operations::extrema::curve_curve::{closest_from, curve_size, is_box_minimum};
use crate::operations::extrema::surface_search::{
    Carried, Foot, Samples, above_rounding, closed_gap, foot, holds_off_focus, inside,
    merge_basins, one_basin, outer_window, pair_scale, probe_length, reference_parameters,
    surface_domains, surface_size,
};
use crate::operations::extrema::{Extremum, parameters_apart};
use crate::operations::intersect::numeric::parameter_domain::ParameterDomain;
use crate::operations::intersect::{Completeness, CurveSurfaceRegion};
use crate::{Circle, Cylinder, Line, Plane, Point, Sphere, Tolerance};

/// One local minimum of the distance between a curve and a surface (see
/// [`curve_surface_extrema`](crate::curve_surface_extrema)).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfaceExtremum {
    /// Parameter on the curve.
    pub curve_parameter: f64,
    /// `(u, v)` parameters on the surface.
    pub surface_parameters: (f64, f64),
    /// The point of the curve at `curve_parameter`.
    pub curve_point: Point,
    /// The point of the surface at `surface_parameters`.
    pub surface_point: Point,
    /// Distance between the two points.
    pub distance: f64,
}

impl CurveSurfaceExtremum {
    /// The extremum between `curve_point`, at `curve_parameter` on the
    /// curve, and `surface_point`, at `surface_parameters` on the surface.
    pub(crate) fn new(
        curve_parameter: f64,
        curve_point: Point,
        [u, v]: [f64; 2],
        surface_point: Point,
    ) -> Self {
        CurveSurfaceExtremum {
            curve_parameter,
            surface_parameters: (u, v),
            curve_point,
            surface_point,
            distance: curve_point.distance(surface_point),
        }
    }

    fn surface_parameter_array(&self) -> [f64; 2] {
        [self.surface_parameters.0, self.surface_parameters.1]
    }
}

impl Extremum for CurveSurfaceExtremum {
    fn distance(&self) -> f64 {
        self.distance
    }

    fn parameters_apart(&self, other: &Self, tol: Tolerance) -> bool {
        parameters_apart(self.curve_parameter, other.curve_parameter, tol)
            || parameters_apart(self.surface_parameters.0, other.surface_parameters.0, tol)
            || parameters_apart(self.surface_parameters.1, other.surface_parameters.1, tol)
    }
}

/// The local minima of the distance between a curve and a surface,
/// closest first, and whether they are all of them (see
/// [`curve_surface_extrema`](crate::curve_surface_extrema)).
///
/// It reads as the slice of its minima. The first, [`closest`](Self::closest),
/// is the minimum distance between the two shapes and the two points
/// where it is reached.
///
/// # Examples
///
/// ```
/// use geomcore::{Line, Point, Sphere, Tolerance, Vector, curve_surface_extrema};
/// let line = Line::new(Point::new(0.0, 0.0, 3.0), Vector::X).unwrap();
/// let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
/// let extrema = curve_surface_extrema(&line, &sphere, Tolerance::DEFAULT);
/// assert!(extrema.is_complete());
/// let closest = extrema.closest().unwrap();
/// assert_eq!(closest.distance, 2.0);
/// assert_eq!(closest.curve_point, Point::new(0.0, 0.0, 3.0));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct CurveSurfaceExtrema {
    minima: Vec<CurveSurfaceExtremum>,
    completeness: Completeness<CurveSurfaceRegion>,
}

impl CurveSurfaceExtrema {
    /// `minima`, put in order of distance, with their completeness.
    pub(crate) fn new(
        mut minima: Vec<CurveSurfaceExtremum>,
        completeness: Completeness<CurveSurfaceRegion>,
    ) -> Self {
        minima.sort_by(|x, y| x.distance.total_cmp(&y.distance));
        CurveSurfaceExtrema {
            minima,
            completeness,
        }
    }

    /// The extrema made of the parts their getters return: for rebuilding
    /// stored extrema, as Python's `pickle` does. The parts are kept as
    /// they are, not checked or sorted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{CurveSurfaceExtrema, Line, Plane, Point, Tolerance, Vector, curve_surface_extrema};
    /// let line = Line::new(Point::ORIGIN, Vector::new(1.0, 0.0, 1.0)).unwrap();
    /// let plane = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
    /// let extrema = curve_surface_extrema(&line, &plane, Tolerance::DEFAULT);
    /// let rebuilt =
    ///     CurveSurfaceExtrema::from_parts(extrema.minima().to_vec(), extrema.completeness().clone());
    /// assert_eq!(rebuilt, extrema);
    /// ```
    pub fn from_parts(
        minima: Vec<CurveSurfaceExtremum>,
        completeness: Completeness<CurveSurfaceRegion>,
    ) -> Self {
        CurveSurfaceExtrema {
            minima,
            completeness,
        }
    }

    /// The minima, closest first.
    pub fn minima(&self) -> &[CurveSurfaceExtremum] {
        &self.minima
    }

    /// The closest minimum, the first of [`minima`](Self::minima): the
    /// minimum distance between the curve and the surface, with the points
    /// where it is reached. `None` when there is no minimum.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Plane, Point, Tolerance, Vector, curve_surface_extrema};
    /// // A line parallel to a plane, 2 above it.
    /// let line = Line::new(Point::new(0.0, 0.0, 2.0), Vector::X).unwrap();
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let extrema = curve_surface_extrema(&line, &plane, Tolerance::DEFAULT);
    /// assert_eq!(extrema.closest().unwrap().distance, 2.0);
    /// ```
    pub fn closest(&self) -> Option<&CurveSurfaceExtremum> {
        self.minima.first()
    }

    /// The minima, closest first, without the completeness.
    pub fn into_minima(self) -> Vec<CurveSurfaceExtremum> {
        self.minima
    }

    /// Whether the minima are all of them, and if not, why and where
    /// minima may be missing.
    pub fn completeness(&self) -> &Completeness<CurveSurfaceRegion> {
        &self.completeness
    }

    /// Whether the minima are all of them ([`Completeness::Complete`]).
    pub fn is_complete(&self) -> bool {
        self.completeness.is_complete()
    }
}

impl Deref for CurveSurfaceExtrema {
    type Target = [CurveSurfaceExtremum];

    fn deref(&self) -> &[CurveSurfaceExtremum] {
        &self.minima
    }
}

impl IntoIterator for CurveSurfaceExtrema {
    type Item = CurveSurfaceExtremum;
    type IntoIter = std::vec::IntoIter<CurveSurfaceExtremum>;

    fn into_iter(self) -> std::vec::IntoIter<CurveSurfaceExtremum> {
        self.minima.into_iter()
    }
}

impl<'a> IntoIterator for &'a CurveSurfaceExtrema {
    type Item = &'a CurveSurfaceExtremum;
    type IntoIter = std::slice::Iter<'a, CurveSurfaceExtremum>;

    fn into_iter(self) -> std::slice::Iter<'a, CurveSurfaceExtremum> {
        self.minima.iter()
    }
}

/// Every local minimum of the distance `|C(t) - S(u, v)|` between a curve
/// and a surface over their parameter domains, closest first, each with
/// the parameters and points on both and the distance. The first is the
/// minimum distance between them.
///
/// A whole line against a whole plane, sphere or cylinder, and a whole
/// circle against a whole plane, are solved in closed form, and their
/// minima are [`Completeness::Complete`](crate::Completeness::Complete):
/// the foot of the line on a sphere or a cylinder's axis, the point of a
/// circle nearest a plane, or the points where the curve crosses the
/// surface, each a minimum at distance `0`. A tangent contact, where the
/// curve stays within `tol.confusion` of the surface across its closest
/// point, is one minimum there.
///
/// Every other pair, a trimmed curve or surface included, is searched
/// numerically and its minima are
/// [`Completeness::Sampled`](crate::Completeness::Sampled): the curve and
/// a grid on the surface are sampled (an unbounded parameter
/// geometrically from its reference parameter, out to many times the
/// pair's own length scale), each sample is paired with the nearest
/// sample of the other shape, and the samples where that distance has a
/// discrete minimum, seen from the curve and from the surface, start a
/// bounded Newton minimisation of the squared distance. The surface
/// parameters follow the curve's as the foot of the curve point on the
/// surface. A result is kept when it is a minimum over the parameter box:
/// stationary, or held at a bound of a trimmed parameter by an outward
/// gradient, with no direction of descent.
///
/// A continuous family of equal minima (a line parallel to a plane, a
/// curve lying on the surface) is reported once, with a representative:
/// the curve at its reference parameter (`0` on a line or circle, the
/// start of a bounded domain) and its closest point on the surface, when
/// that pair belongs to the family, else a member found by the search.
///
/// # Examples
///
/// ```
/// use geomcore::{Circle, Plane, Point, Tolerance, Vector, curve_surface_extrema};
/// // A circle of radius 1 tilted 45° about the x axis, its centre 2 above a plane.
/// let circle = Circle::new(Point::new(0.0, 0.0, 2.0), Vector::new(0.0, 1.0, 1.0), 1.0).unwrap();
/// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
/// let extrema = curve_surface_extrema(&circle, &plane, Tolerance::DEFAULT);
/// assert_eq!(extrema.len(), 1);
/// assert!((extrema[0].distance - (2.0 - 0.5f64.sqrt())).abs() < 1e-12);
/// ```
pub fn curve_surface_extrema<C, S>(curve: &C, surface: &S, tol: Tolerance) -> CurveSurfaceExtrema
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    if let Some(parameters) = closed_form_parameters(curve, surface, tol) {
        let minima = parameters
            .into_iter()
            .map(|t| projected(curve, surface, t, tol))
            .collect();
        return CurveSurfaceExtrema::new(minima, Completeness::Complete);
    }
    sampled_curve_surface_extrema(curve, surface, tol)
}

/// [`curve_surface_extrema`] through the numeric search alone, never a
/// closed form: the reference the closed forms are tested against.
pub(crate) fn sampled_curve_surface_extrema<C, S>(
    curve: &C,
    surface: &S,
    tol: Tolerance,
) -> CurveSurfaceExtrema
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    CurveSurfaceExtrema::new(sampled_minima(curve, surface, tol), Completeness::Sampled)
}

/// The curve as a whole elementary curve: its elementary basis when the
/// curve spans all of it, not a trimmed piece.
pub(crate) fn whole_elementary_curve<C: ParametricCurve>(curve: &C) -> Option<ElementaryCurve> {
    let elementary = curve.as_elementary()?;
    (curve.bounds() == elementary.bounds() && curve.period() == elementary.period())
        .then_some(elementary)
}

/// The surface as a whole elementary surface: its elementary basis when
/// the surface spans all of it, not a trimmed patch.
pub(crate) fn whole_elementary_surface<S: ParametricSurface>(
    surface: &S,
) -> Option<ElementarySurface> {
    let elementary = surface.as_elementary()?;
    let whole = surface.u_bounds() == elementary.u_bounds()
        && surface.v_bounds() == elementary.v_bounds()
        && surface.u_period() == elementary.u_period()
        && surface.v_period() == elementary.v_period();
    whole.then_some(elementary)
}

/// The extremum of the curve point at `t` and its closest point on the
/// surface.
fn projected<C, S>(curve: &C, surface: &S, t: f64, tol: Tolerance) -> CurveSurfaceExtremum
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    let point = curve.point_at(t);
    let projection = surface.project_point(point, tol);
    CurveSurfaceExtremum::new(t, point, [projection.u, projection.v], projection.point)
}

/// The curve parameters of the minima of a pair with a closed form, each
/// projected onto the surface for its point there, or `None` when the
/// pair has none.
fn closed_form_parameters<C, S>(curve: &C, surface: &S, tol: Tolerance) -> Option<Vec<f64>>
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    let curve = whole_elementary_curve(curve)?;
    let surface = whole_elementary_surface(surface)?;
    match (curve, surface) {
        (ElementaryCurve::Line(line), ElementarySurface::Plane(plane)) => {
            Some(line_plane(&line, &plane, tol))
        }
        (ElementaryCurve::Line(line), ElementarySurface::Sphere(sphere)) => {
            Some(line_sphere(&line, &sphere, tol))
        }
        (ElementaryCurve::Line(line), ElementarySurface::Cylinder(cylinder)) => {
            Some(line_cylinder(&line, &cylinder, tol))
        }
        (ElementaryCurve::Circle(circle), ElementarySurface::Plane(plane)) => {
            Some(circle_plane(&circle, &plane, tol))
        }
        _ => None,
    }
}

/// The crossing of a line with a plane, or `0`, the line's reference
/// parameter, standing for a line parallel to the plane.
fn line_plane(line: &Line, plane: &Plane, tol: Tolerance) -> Vec<f64> {
    let normal = plane.axis().direction();
    let along = line.direction();
    if along.is_perpendicular(normal, tol.angular) {
        return vec![0.0];
    }
    vec![-plane.signed_distance(line.origin()) / along.dot(normal)]
}

/// The parameters about `foot` where a line at the distance `gap` from a
/// centre (a sphere's) or an axis (a cylinder's) of radius `radius` crosses
/// the surface, with `rate` the rate at which the line moves away from the
/// centre or axis: two crossings, or the foot alone where the line stays
/// out of the surface or touches it within `tol.confusion`.
fn crossings_about(foot: f64, gap: f64, radius: f64, rate: f64, tol: Tolerance) -> Vec<f64> {
    if radius - gap > tol.confusion {
        let half = ((radius - gap) * (radius + gap)).sqrt() / rate;
        vec![foot - half, foot + half]
    } else {
        vec![foot]
    }
}

/// The minima of a line and a sphere: the foot of the centre on the line,
/// or the two crossings.
fn line_sphere(line: &Line, sphere: &Sphere, tol: Tolerance) -> Vec<f64> {
    let foot = (sphere.center() - line.origin()).dot(line.direction());
    let gap = line.point_at(foot).distance(sphere.center());
    crossings_about(foot, gap, sphere.radius(), 1.0, tol)
}

/// The minima of a line and a cylinder: the line's point nearest the axis,
/// or the two crossings; `0`, the line's reference parameter, stands for a
/// line parallel to the axis.
fn line_cylinder(line: &Line, cylinder: &Cylinder, tol: Tolerance) -> Vec<f64> {
    let axis = cylinder.axis();
    let (along, across) = (line.direction(), axis.direction());
    if along.is_parallel(across, tol.angular) {
        return vec![0.0];
    }
    let normal = along.cross(across);
    let sine = normal.magnitude();
    let offset = axis.origin() - line.origin();
    let foot = offset.cross(across).dot(normal) / (sine * sine);
    let gap = offset.dot(normal).abs() / sine;
    crossings_about(foot, gap, cylinder.radius(), sine, tol)
}

/// The minima of a circle and a plane: the circle's point nearest the
/// plane, or the two crossings; `0`, the circle's reference parameter,
/// stands for a circle parallel to the plane. The circle's height above
/// the plane is `h + A cos(θ - φ)`.
fn circle_plane(circle: &Circle, plane: &Plane, tol: Tolerance) -> Vec<f64> {
    let normal = plane.axis().direction();
    if circle.normal().is_parallel(normal, tol.angular) {
        return vec![0.0];
    }
    let frame = circle.frame();
    let (a, b) = (
        frame.x_direction().dot(normal),
        frame.y_direction().dot(normal),
    );
    let amplitude = circle.radius() * a.hypot(b);
    let phase = b.atan2(a);
    let height = plane.signed_distance(circle.center());
    let angles = if amplitude - height.abs() > tol.confusion {
        let half = (-height / amplitude).acos();
        vec![phase - half, phase + half]
    } else if height > 0.0 {
        vec![phase + TAU / 2.0]
    } else {
        vec![phase]
    };
    angles
        .into_iter()
        .map(|angle| angle.rem_euclid(TAU))
        .collect()
}

/// The minima of the numeric search (see [`curve_surface_extrema`]).
fn sampled_minima<C, S>(curve: &C, surface: &S, tol: Tolerance) -> Vec<CurveSurfaceExtremum>
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    let curve_domain = ParameterDomain::of(curve);
    let surface_domains = surface_domains(surface);
    let [u0, v0] = reference_parameters(surface_domains);
    let scale = pair_scale(
        curve.point_at(curve_domain.reference()),
        curve_size(curve, curve_domain),
        surface.point_at(u0, v0),
        surface_size(surface, surface_domains),
    );
    let curve_samples = Samples::of_curve(curve, curve_domain, scale);
    let surface_samples = Samples::of_surface(surface, surface_domains, scale);
    let search = Search {
        curve,
        curve_domain,
        surface,
        surface_domains,
        scale,
        tol,
    };
    // Seen from the curve, each sample measured to its foot on the surface
    // from the nearest surface sample: a distance exact at every sample, so
    // that minima a few samples apart are told apart.
    let from_curve = curve_samples.discrete_minima(|point| {
        let (nearest, _) = surface_samples.nearest(point);
        let foot = foot(
            surface,
            surface_domains,
            point,
            surface_samples.parameters[nearest],
        );
        (foot.point.distance(point), foot.parameters)
    });
    let from_surface = surface_samples.discrete_minima(|point| {
        let (nearest, distance) = curve_samples.nearest(point);
        (distance, nearest)
    });
    let seeds = from_curve
        .into_iter()
        .map(|(on_curve, start)| (curve_samples.parameters[on_curve][0], start))
        .chain(from_surface.into_iter().map(|(on_surface, on_curve)| {
            (
                curve_samples.parameters[on_curve][0],
                surface_samples.parameters[on_surface],
            )
        }));
    let found: Vec<CurveSurfaceExtremum> = seeds
        .filter_map(|(t, start)| search.refine(t, start))
        .collect();
    merge_basins(
        found,
        |p, q| search.same_basin(p, q),
        |e| search.representative(e),
        tol,
    )
}

/// A curve and a surface with their domains, for the numeric search.
struct Search<'a, C, S> {
    curve: &'a C,
    curve_domain: ParameterDomain,
    surface: &'a S,
    surface_domains: [ParameterDomain; 2],
    /// The configuration's length scale, for the outer window.
    scale: f64,
    tol: Tolerance,
}

impl<C: ParametricCurve, S: ParametricSurface> Search<'_, C, S> {
    /// The extremum of the curve point at `t` and its foot on the surface
    /// from `start`, with half the squared distance and its first and
    /// reduced second derivatives along `t`, the curve's speed and the
    /// foot itself.
    fn local(&self, t: f64, start: [f64; 2]) -> (CurveSurfaceExtremum, (f64, f64, f64), f64, Foot) {
        let (point, first, second) = self.curve.second_derivatives(t);
        let foot = foot(self.surface, self.surface_domains, point, start);
        let r = point - foot.point;
        let [[reduced]] = foot.reduced(
            [[first.dot(first) + r.dot(second)]],
            [[-first.dot(foot.first[0]), -first.dot(foot.first[1])]],
        );
        let extremum = CurveSurfaceExtremum::new(t, point, foot.parameters, foot.point);
        (
            extremum,
            (
                0.5 * r.dot(r),
                above_rounding(r.dot(first), point, foot.point, first),
                reduced,
            ),
            first.magnitude(),
            foot,
        )
    }

    /// The minimum reached from the seed `(t, start)`, when it is one over
    /// the parameter box.
    fn refine(&self, t: f64, start: [f64; 2]) -> Option<CurveSurfaceExtremum> {
        let mut warm = start;
        let t = newton_min_1d(
            |t| {
                let (extremum, objective, ..) = self.local(t, warm);
                warm = extremum.surface_parameter_array();
                objective
            },
            t,
            {
                let speed = self.curve.first_derivatives(t).1.magnitude();
                outer_window(self.curve_domain, t, speed, self.scale)
            },
            1e-14,
            100,
        );
        let t = self.curve_domain.wrap(t);
        if !t.is_finite() {
            return None;
        }
        let (extremum, (_, gradient, reduced), speed, found) = self.local(t, warm);
        let along_curve = is_box_minimum(
            [t, 0.0],
            [gradient, 0.0],
            [[reduced, 0.0], [0.0, 0.0]],
            [speed, 0.0],
            [self.curve_domain.bounds(), None],
            self.tol,
        );
        let finite = extremum.surface_point.x.is_finite()
            && extremum.surface_point.y.is_finite()
            && extremum.surface_point.z.is_finite();
        let on_surface = found.is_minimum(self.surface, self.surface_domains, self.tol);
        let holds = !found.is_degenerate() || {
            let step = probe_length(self.tol, self.scale) / speed;
            let probes = [t - step, t + step]
                .into_iter()
                .filter(|&probe| inside(self.curve_domain, probe).is_some())
                .map(|probe| {
                    let point = self.curve.point_at(probe);
                    let start = extremum.surface_parameter_array();
                    foot(self.surface, self.surface_domains, point, start)
                        .point
                        .distance(point)
                });
            holds_off_focus(extremum.distance, probes, self.tol)
        };
        (finite && along_curve && on_surface && holds)
            .then(|| self.closed(&extremum).unwrap_or(extremum))
    }

    /// The crossing a minimum lies next to, when the gap between the curve
    /// and the surface closes from it (see `closed_gap`).
    fn closed(&self, e: &CurveSurfaceExtremum) -> Option<CurveSurfaceExtremum> {
        let (u, v) = e.surface_parameters;
        let [t, u, v] = closed_gap(
            [e.curve_parameter, u, v],
            |[t, u, v]| {
                let (c, along) = self.curve.first_derivatives(t);
                let (s, [su, sv]) = self.surface.first_derivatives(u, v);
                (c - s, [along, -su, -sv])
            },
            self.tol,
        )?;
        let t = inside(self.curve_domain, t)?;
        let surface = [
            inside(self.surface_domains[0], u)?,
            inside(self.surface_domains[1], v)?,
        ];
        Some(CurveSurfaceExtremum::new(
            t,
            self.curve.point_at(t),
            surface,
            self.surface.point_at(surface[0], surface[1]),
        ))
    }

    /// Whether two minima lie in one basin: walking the curve from one to
    /// the other, each foot warm-started from the last, finds no ridge
    /// (see `one_basin`) and arrives on the other's surface point; either
    /// end may start the walk, as the feet follow one sheet of the surface.
    /// Or walking the surface along the straight way in its parameters,
    /// each closest point on the curve found from the last, for a valley
    /// along a piece of the surface that the shorter way round a closed
    /// curve leaves (a whole circle against a piece of a coaxial
    /// cylinder). Or, for a valley of the surface over one curve point, no
    /// ridge along the straight way in the parameters of both shapes.
    fn same_basin(&self, p: &CurveSurfaceExtremum, q: &CurveSurfaceExtremum) -> bool {
        self.walks_to(p, q)
            || self.walks_to(q, p)
            || self.walks_on_surface_to(p, q)
            || self.straight_between(p, q)
    }

    /// Whether the walk on the surface from `p`'s surface point to `q`'s,
    /// along the straight way in its parameters, keeps the closest point
    /// on the curve in one basin and arrives on `q`'s curve point.
    fn walks_on_surface_to(&self, p: &CurveSurfaceExtremum, q: &CurveSurfaceExtremum) -> bool {
        let (from, to) = (p.surface_parameter_array(), q.surface_parameter_array());
        let way = [0, 1].map(|i| self.surface_domains[i].step(from[i], to[i]));
        let mut warm = p.curve_parameter;
        let mut last = p.curve_point;
        let no_ridge = one_basin(
            p.distance,
            q.distance,
            |f| {
                let on_surface = self
                    .surface
                    .point_at(from[0] + f * way[0], from[1] + f * way[1]);
                let (t, on_curve) = closest_from(self.curve, self.curve_domain, on_surface, warm);
                (warm, last) = (t, on_curve);
                on_curve.distance(on_surface)
            },
            self.tol,
        );
        no_ridge && {
            let end = self.surface.point_at(from[0] + way[0], from[1] + way[1]);
            let (_, on_curve) = closest_from(self.curve, self.curve_domain, end, warm);
            last = on_curve;
            last.distance(q.curve_point) <= self.tol.confusion
        }
    }

    /// Whether no ridge parts `p` from `q` along the straight way in the
    /// parameters of both shapes, each periodic one the shorter way round.
    fn straight_between(&self, p: &CurveSurfaceExtremum, q: &CurveSurfaceExtremum) -> bool {
        let dt = self.curve_domain.step(p.curve_parameter, q.curve_parameter);
        let (from, to) = (p.surface_parameter_array(), q.surface_parameter_array());
        let way = [0, 1].map(|i| self.surface_domains[i].step(from[i], to[i]));
        one_basin(
            p.distance,
            q.distance,
            |f| {
                let on_surface = self
                    .surface
                    .point_at(from[0] + f * way[0], from[1] + f * way[1]);
                self.curve
                    .point_at(p.curve_parameter + f * dt)
                    .distance(on_surface)
            },
            self.tol,
        )
    }

    /// Whether the walk from `p` to `q` stays in one basin and arrives.
    fn walks_to(&self, p: &CurveSurfaceExtremum, q: &CurveSurfaceExtremum) -> bool {
        let dt = self.curve_domain.step(p.curve_parameter, q.curve_parameter);
        let mut feet = Carried::new(
            self.surface,
            self.surface_domains,
            |f| self.curve.point_at(p.curve_parameter + f * dt),
            p.surface_parameter_array(),
            p.curve_point.distance(q.curve_point),
            self.scale,
        );
        let no_ridge = one_basin(
            p.distance,
            q.distance,
            |f| {
                let foot = feet.foot_at(f);
                foot.point
                    .distance(self.curve.point_at(p.curve_parameter + f * dt))
            },
            self.tol,
        );
        no_ridge && feet.foot_at(1.0).point.distance(q.surface_point) <= self.tol.confusion
    }

    /// A valley's representative: the curve at its reference parameter and
    /// its foot on the surface, when that pair lies in the valley; else the
    /// curve's closest point to the surface at its reference parameters,
    /// and its foot, for a valley that starts where a piece of the surface
    /// starts (a whole circle against a piece of a coaxial cylinder); else
    /// the minimum the search reaches from the curve's reference: along a
    /// valley the gradient is rounding alone, so the curve parameter stays
    /// there. Each is read from the shapes alone, not from where the search
    /// entered the valley (ADR 0051).
    fn representative(&self, e: &CurveSurfaceExtremum) -> CurveSurfaceExtremum {
        let in_valley = |pair: &CurveSurfaceExtremum| {
            pair.distance <= e.distance + self.tol.confusion && self.same_basin(e, pair)
        };
        let t = self.curve_domain.reference();
        let projection = self.surface.project_point(self.curve.point_at(t), self.tol);
        let start = [projection.u, projection.v];
        let pair = self.local(t, start).0;
        if in_valley(&pair) {
            return pair;
        }
        let [u, v] = reference_parameters(self.surface_domains);
        let surface_reference = self.surface.point_at(u, v);
        let (near, on_curve) = closest_from(
            self.curve,
            self.curve_domain,
            surface_reference,
            e.curve_parameter,
        );
        let projection = self.surface.project_point(on_curve, self.tol);
        let from_surface = self.local(near, [projection.u, projection.v]).0;
        if in_valley(&from_surface) {
            return from_surface;
        }
        self.refine(t, start).unwrap_or(*e)
    }
}

#[cfg(test)]
mod tests;
