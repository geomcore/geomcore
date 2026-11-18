//! The local minima of the distance between two surfaces
//! ([`SurfaceSurfaceExtrema`]).

use std::ops::Deref;

use crate::geometry::surfaces::{ElementarySurface, ParametricSurface};
use crate::math::numeric::newton_min_2d;
use crate::operations::extrema::curve_curve::is_box_minimum;
use crate::operations::extrema::curve_surface::whole_elementary_surface;
use crate::operations::extrema::surface_search::{
    Carried, Foot, Objective, Samples, above_rounding, closed_gap, collapsed_minimum, foot,
    holds_off_focus, inside, merge_basins, one_basin, outer_window, pair_scale, probe_length,
    reference_parameters, seeds, surface_domains, surface_size,
};
use crate::operations::extrema::{Extremum, parameters_apart};
use crate::operations::intersect::numeric::parameter_domain::ParameterDomain;
use crate::operations::intersect::{Completeness, SurfaceSurfaceRegion};
use crate::{Cylinder, Plane, Point, Sphere, Tolerance, Vector};

/// One local minimum of the distance between two surfaces (see
/// [`surface_surface_extrema`](crate::surface_surface_extrema)).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceSurfaceExtremum {
    /// `(u, v)` parameters on the first surface.
    pub first_parameters: (f64, f64),
    /// `(u, v)` parameters on the second surface.
    pub second_parameters: (f64, f64),
    /// The point of the first surface at `first_parameters`.
    pub first_point: Point,
    /// The point of the second surface at `second_parameters`.
    pub second_point: Point,
    /// Distance between the two points.
    pub distance: f64,
}

impl SurfaceSurfaceExtremum {
    /// The extremum between `first_point`, at `first_parameters` on the
    /// first surface, and `second_point`, at `second_parameters` on the
    /// second.
    pub(crate) fn new(
        [u1, v1]: [f64; 2],
        first_point: Point,
        [u2, v2]: [f64; 2],
        second_point: Point,
    ) -> Self {
        SurfaceSurfaceExtremum {
            first_parameters: (u1, v1),
            second_parameters: (u2, v2),
            first_point,
            second_point,
            distance: first_point.distance(second_point),
        }
    }

    fn first_parameter_array(&self) -> [f64; 2] {
        [self.first_parameters.0, self.first_parameters.1]
    }

    fn second_parameter_array(&self) -> [f64; 2] {
        [self.second_parameters.0, self.second_parameters.1]
    }
}

impl Extremum for SurfaceSurfaceExtremum {
    fn distance(&self) -> f64 {
        self.distance
    }

    fn parameters_apart(&self, other: &Self, tol: Tolerance) -> bool {
        let (p, q) = (self, other);
        parameters_apart(p.first_parameters.0, q.first_parameters.0, tol)
            || parameters_apart(p.first_parameters.1, q.first_parameters.1, tol)
            || parameters_apart(p.second_parameters.0, q.second_parameters.0, tol)
            || parameters_apart(p.second_parameters.1, q.second_parameters.1, tol)
    }
}

/// The local minima of the distance between two surfaces, closest first,
/// and whether they are all of them (see
/// [`surface_surface_extrema`](crate::surface_surface_extrema)).
///
/// It reads as the slice of its minima. The first, [`closest`](Self::closest),
/// is the minimum distance between the two surfaces and the two points
/// where it is reached.
///
/// # Examples
///
/// ```
/// use geomcore::{Point, Sphere, Tolerance, surface_surface_extrema};
/// let first = Sphere::new(Point::ORIGIN, 1.0).unwrap();
/// let second = Sphere::new(Point::new(4.0, 0.0, 0.0), 2.0).unwrap();
/// let extrema = surface_surface_extrema(&first, &second, Tolerance::DEFAULT);
/// assert!(extrema.is_complete());
/// let closest = extrema.closest().unwrap();
/// assert_eq!(closest.distance, 1.0);
/// assert_eq!(closest.first_point, Point::new(1.0, 0.0, 0.0));
/// assert_eq!(closest.second_point, Point::new(2.0, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceSurfaceExtrema {
    minima: Vec<SurfaceSurfaceExtremum>,
    completeness: Completeness<SurfaceSurfaceRegion>,
}

impl SurfaceSurfaceExtrema {
    /// `minima`, put in order of distance, with their completeness.
    pub(crate) fn new(
        mut minima: Vec<SurfaceSurfaceExtremum>,
        completeness: Completeness<SurfaceSurfaceRegion>,
    ) -> Self {
        minima.sort_by(|x, y| x.distance.total_cmp(&y.distance));
        SurfaceSurfaceExtrema {
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
    /// use geomcore::{Plane, Point, SurfaceSurfaceExtrema, Tolerance, Vector, surface_surface_extrema};
    /// let first = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let second = Plane::new(Point::new(0.0, 0.0, 3.0), Vector::Z).unwrap();
    /// let extrema = surface_surface_extrema(&first, &second, Tolerance::DEFAULT);
    /// let rebuilt =
    ///     SurfaceSurfaceExtrema::from_parts(extrema.minima().to_vec(), extrema.completeness().clone());
    /// assert_eq!(rebuilt, extrema);
    /// ```
    pub fn from_parts(
        minima: Vec<SurfaceSurfaceExtremum>,
        completeness: Completeness<SurfaceSurfaceRegion>,
    ) -> Self {
        SurfaceSurfaceExtrema {
            minima,
            completeness,
        }
    }

    /// The minima, closest first.
    pub fn minima(&self) -> &[SurfaceSurfaceExtremum] {
        &self.minima
    }

    /// The closest minimum, the first of [`minima`](Self::minima): the
    /// minimum distance between the surfaces, with the points where it is
    /// reached. `None` when there is no minimum.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Tolerance, Vector, surface_surface_extrema};
    /// let first = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let second = Plane::new(Point::new(0.0, 0.0, 3.0), Vector::Z).unwrap();
    /// let extrema = surface_surface_extrema(&first, &second, Tolerance::DEFAULT);
    /// assert_eq!(extrema.closest().unwrap().distance, 3.0);
    /// ```
    pub fn closest(&self) -> Option<&SurfaceSurfaceExtremum> {
        self.minima.first()
    }

    /// The minima, closest first, without the completeness.
    pub fn into_minima(self) -> Vec<SurfaceSurfaceExtremum> {
        self.minima
    }

    /// Whether the minima are all of them, and if not, why and where
    /// minima may be missing.
    pub fn completeness(&self) -> &Completeness<SurfaceSurfaceRegion> {
        &self.completeness
    }

    /// Whether the minima are all of them ([`Completeness::Complete`]).
    pub fn is_complete(&self) -> bool {
        self.completeness.is_complete()
    }
}

impl Deref for SurfaceSurfaceExtrema {
    type Target = [SurfaceSurfaceExtremum];

    fn deref(&self) -> &[SurfaceSurfaceExtremum] {
        &self.minima
    }
}

impl IntoIterator for SurfaceSurfaceExtrema {
    type Item = SurfaceSurfaceExtremum;
    type IntoIter = std::vec::IntoIter<SurfaceSurfaceExtremum>;

    fn into_iter(self) -> std::vec::IntoIter<SurfaceSurfaceExtremum> {
        self.minima.into_iter()
    }
}

impl<'a> IntoIterator for &'a SurfaceSurfaceExtrema {
    type Item = &'a SurfaceSurfaceExtremum;
    type IntoIter = std::slice::Iter<'a, SurfaceSurfaceExtremum>;

    fn into_iter(self) -> std::slice::Iter<'a, SurfaceSurfaceExtremum> {
        self.minima.iter()
    }
}

/// Every local minimum of the distance `|S1(u1, v1) - S2(u2, v2)|`
/// between two surfaces over their parameter domains, closest first, each
/// with the parameters and points on both and the distance. The first is
/// the minimum distance between them.
///
/// **Surfaces that meet** have a minimum at distance `0` all along each
/// curve where they meet. Those are reported as one minimum, a point
/// where they meet: the distance is `0`, and
/// [`Intersect`](crate::Intersect) gives the curves.
///
/// Two whole planes, two whole spheres, and two whole cylinders with
/// parallel axes are solved in closed form, and their minima are
/// [`Completeness::Complete`](crate::Completeness::Complete). Every other
/// pair, a trimmed surface included, is searched numerically and its
/// minima are [`Completeness::Sampled`](crate::Completeness::Sampled): a
/// grid on each surface is sampled (an unbounded parameter geometrically
/// from its reference parameter, out to many times the pair's own length
/// scale), each sample is paired with the nearest sample of the other
/// surface, and the samples where that distance has a discrete minimum,
/// seen from either surface, start a bounded Newton minimisation of the
/// squared distance. The second surface's parameters follow the first's
/// as the foot of the first surface's point on the second. A result is
/// kept when it is a minimum over the parameter boxes: stationary, or held
/// at a bound of a trimmed parameter by an outward gradient, with no
/// direction of descent.
///
/// A continuous family of equal minima (two parallel planes, concentric
/// spheres, parallel cylinders) is reported once, with a representative:
/// the first surface at its reference parameters and its closest point on
/// the second, when that pair belongs to the family, else a member found
/// by the closed form or the search.
///
/// # Examples
///
/// ```
/// use geomcore::{Cylinder, Point, Tolerance, Vector, surface_surface_extrema};
/// // Two cylinders of radius 1 along z, their axes 5 apart.
/// let first = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
/// let second = Cylinder::new(Point::new(5.0, 0.0, 0.0), Vector::Z, 1.0).unwrap();
/// let extrema = surface_surface_extrema(&first, &second, Tolerance::DEFAULT);
/// assert_eq!(extrema.len(), 1);
/// assert_eq!(extrema[0].distance, 3.0);
/// ```
pub fn surface_surface_extrema<A, B>(first: &A, second: &B, tol: Tolerance) -> SurfaceSurfaceExtrema
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    if let Some(minima) = closed_form(first, second, tol) {
        return SurfaceSurfaceExtrema::new(minima, Completeness::Complete);
    }
    sampled_surface_surface_extrema(first, second, tol)
}

/// [`surface_surface_extrema`] through the numeric search alone, never a
/// closed form: the reference the closed forms are tested against.
pub(crate) fn sampled_surface_surface_extrema<A, B>(
    first: &A,
    second: &B,
    tol: Tolerance,
) -> SurfaceSurfaceExtrema
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    SurfaceSurfaceExtrema::new(sampled_minima(first, second, tol), Completeness::Sampled)
}

/// The two surfaces, for the closed forms: each point is read back into
/// its surface's parameters.
struct Pair<'a, A, B> {
    first: &'a A,
    second: &'a B,
    tol: Tolerance,
}

impl<A: ParametricSurface, B: ParametricSurface> Pair<'_, A, B> {
    /// The extremum of `first_point` on the first surface and
    /// `second_point` on the second.
    fn at(&self, first_point: Point, second_point: Point) -> SurfaceSurfaceExtremum {
        let (u1, v1) = self.first.parameters_of(first_point, self.tol);
        let (u2, v2) = self.second.parameters_of(second_point, self.tol);
        SurfaceSurfaceExtremum::new([u1, v1], first_point, [u2, v2], second_point)
    }

    /// The first surface's point at its reference parameters.
    fn reference_point(&self) -> Point {
        let [u, v] = reference_parameters(surface_domains(self.first));
        self.first.point_at(u, v)
    }

    /// The representative of a family of equal minima: the first surface
    /// at its reference parameters and its closest point on the second.
    fn reference_pair(&self) -> SurfaceSurfaceExtremum {
        let [u, v] = reference_parameters(surface_domains(self.first));
        let point = self.first.point_at(u, v);
        let projection = self.second.project_point(point, self.tol);
        SurfaceSurfaceExtremum::new(
            [u, v],
            point,
            [projection.u, projection.v],
            projection.point,
        )
    }
}

/// The minima of a pair with a closed form, or `None` when it has none.
fn closed_form<A, B>(first: &A, second: &B, tol: Tolerance) -> Option<Vec<SurfaceSurfaceExtremum>>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let pair = Pair { first, second, tol };
    match (
        whole_elementary_surface(first)?,
        whole_elementary_surface(second)?,
    ) {
        (ElementarySurface::Plane(a), ElementarySurface::Plane(b)) => {
            Some(vec![plane_plane(&pair, &a, &b)])
        }
        (ElementarySurface::Sphere(a), ElementarySurface::Sphere(b)) => {
            Some(vec![sphere_sphere(&pair, &a, &b)])
        }
        (ElementarySurface::Cylinder(a), ElementarySurface::Cylinder(b)) => {
            parallel_cylinders(&pair, &a, &b).map(|minimum| vec![minimum])
        }
        _ => None,
    }
}

/// Two planes: parallel ones are a family at their distance apart; others
/// meet, at the point of their line nearest the first plane's origin.
fn plane_plane<A, B>(pair: &Pair<'_, A, B>, first: &Plane, second: &Plane) -> SurfaceSurfaceExtremum
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (n1, n2) = (first.axis().direction(), second.axis().direction());
    if n1.is_parallel(n2, pair.tol.angular) {
        return pair.reference_pair();
    }
    // The point `origin + a n1 + b n2` on both planes, nearest the origin:
    // `a + b c = 0` and `a c + b = e`, `c` the cosine between the normals.
    let origin = first.frame().origin();
    let cosine = n1.dot(n2);
    let offset = n2.dot(second.frame().origin() - origin);
    let b = offset / n1.cross(n2).square_magnitude();
    let point = origin + n1 * (-cosine * b) + n2 * b;
    pair.at(point, point)
}

/// Where two balls (spheres, or a cylinder's cross-section) of radii `r1`
/// and `r2`, their centres `gap` apart along `towards`, come closest: the
/// two points as offsets from the first centre and from the second, or
/// `None` where they meet, or are concentric within `tol.confusion`.
fn closest_offsets(
    gap: f64,
    towards: Vector,
    r1: f64,
    r2: f64,
    tol: Tolerance,
) -> Option<(Vector, Vector)> {
    if gap <= tol.confusion {
        return None;
    }
    if gap >= r1 + r2 - tol.confusion {
        // Apart, or touching from outside.
        Some((towards * r1, towards * -r2))
    } else if gap <= (r1 - r2).abs() + tol.confusion {
        // One inside the other: both points on the far side of the outer.
        let side = if r1 >= r2 { 1.0 } else { -1.0 };
        Some((towards * (side * r1), towards * (side * r2)))
    } else {
        None
    }
}

/// Two spheres: concentric ones are a family; apart, nested or touching
/// ones come closest along the line of centres; others meet on a circle,
/// represented by its point nearest the first sphere's reference point.
fn sphere_sphere<A, B>(
    pair: &Pair<'_, A, B>,
    first: &Sphere,
    second: &Sphere,
) -> SurfaceSurfaceExtremum
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (c1, c2) = (first.center(), second.center());
    let (r1, r2) = (first.radius(), second.radius());
    let gap = c1.distance(c2);
    let Some(towards) = (c2 - c1).direction().filter(|_| gap > pair.tol.confusion) else {
        return pair.reference_pair();
    };
    if let Some((on_first, on_second)) = closest_offsets(gap, towards.vector(), r1, r2, pair.tol) {
        return pair.at(c1 + on_first, c2 + on_second);
    }
    let along = (gap * gap + r1 * r1 - r2 * r2) / (2.0 * gap);
    let radius = (r1 * r1 - along * along).max(0.0).sqrt();
    let frame = first.frame();
    let across = [
        pair.reference_point() - c1,
        frame.x_direction().vector(),
        frame.y_direction().vector(),
    ]
    .map(|v| v.component_perpendicular_to(towards))
    .into_iter()
    .max_by(|a, b| a.magnitude().total_cmp(&b.magnitude()))
    .and_then(Vector::direction);
    let Some(across) = across else {
        return pair.reference_pair();
    };
    let point = c1 + towards * along + across * radius;
    pair.at(point, point)
}

/// Two cylinders with parallel axes, read in the cross-section through
/// the first cylinder's reference point: coaxial ones are a family at
/// the difference of their radii; apart, nested or touching ones come
/// closest along the line between the axes, all along them; others meet
/// along two lines, represented by the point of one in that section.
/// `None` for axes that are not parallel.
fn parallel_cylinders<A, B>(
    pair: &Pair<'_, A, B>,
    first: &Cylinder,
    second: &Cylinder,
) -> Option<SurfaceSurfaceExtremum>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (axis1, axis2) = (first.axis(), second.axis());
    let along_axes = axis1.direction();
    if !along_axes.is_parallel(axis2.direction(), pair.tol.angular) {
        return None;
    }
    let reference = pair.reference_point();
    let centre1 = axis1.closest_point(reference);
    let between = (axis2.origin() - centre1).component_perpendicular_to(along_axes);
    let centre2 = centre1 + between;
    let gap = between.magnitude();
    let (r1, r2) = (first.radius(), second.radius());
    let Some(towards) = between.direction().filter(|_| gap > pair.tol.confusion) else {
        return Some(pair.reference_pair());
    };
    if let Some((on_first, on_second)) = closest_offsets(gap, towards.vector(), r1, r2, pair.tol) {
        return Some(pair.at(centre1 + on_first, centre2 + on_second));
    }
    let along = (gap * gap + r1 * r1 - r2 * r2) / (2.0 * gap);
    let radius = (r1 * r1 - along * along).max(0.0).sqrt();
    let across = along_axes.cross(towards);
    let side = if (reference - centre1).dot(across) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    let point = centre1 + towards * along + across * (side * radius);
    Some(pair.at(point, point))
}

/// The minima of the numeric search (see [`surface_surface_extrema`]).
fn sampled_minima<A, B>(first: &A, second: &B, tol: Tolerance) -> Vec<SurfaceSurfaceExtremum>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let first_domains = surface_domains(first);
    let second_domains = surface_domains(second);
    let [u1, v1] = reference_parameters(first_domains);
    let [u2, v2] = reference_parameters(second_domains);
    let scale = pair_scale(
        first.point_at(u1, v1),
        surface_size(first, first_domains),
        second.point_at(u2, v2),
        surface_size(second, second_domains),
    );
    let first_samples = Samples::of_surface(first, first_domains, scale);
    let second_samples = Samples::of_surface(second, second_domains, scale);
    let search = Search {
        first,
        first_domains,
        second,
        second_domains,
        scale,
        tol,
    };
    let found: Vec<SurfaceSurfaceExtremum> = seeds(&first_samples, &second_samples)
        .into_iter()
        .filter_map(|(on_first, on_second)| {
            search.refine(
                first_samples.parameters[on_first],
                second_samples.parameters[on_second],
            )
        })
        .collect();
    merge_basins(
        found,
        |p, q| search.same_basin(p, q),
        |e| search.representative(e),
        tol,
    )
}

/// Two surfaces with their domains, for the numeric search.
struct Search<'a, A, B> {
    first: &'a A,
    first_domains: [ParameterDomain; 2],
    second: &'a B,
    second_domains: [ParameterDomain; 2],
    /// The configuration's length scale, for the outer windows.
    scale: f64,
    tol: Tolerance,
}

impl<A: ParametricSurface, B: ParametricSurface> Search<'_, A, B> {
    /// The extremum of the first surface's point at `parameters` and its
    /// foot on the second from `start`, with half the squared distance and
    /// its gradient and reduced Hessian over the first surface's
    /// parameters, the first surface's speeds and the foot itself.
    fn local(
        &self,
        [u, v]: [f64; 2],
        start: [f64; 2],
    ) -> (SurfaceSurfaceExtremum, Objective, [f64; 2], Foot) {
        let (point, [pu, pv], [puu, puv, pvv]) = self.first.second_derivatives(u, v);
        let foot = foot(self.second, self.second_domains, point, start);
        let r = point - foot.point;
        let mixed = pu.dot(pv) + r.dot(puv);
        let [su, sv] = foot.first;
        let hessian = foot.reduced(
            [
                [pu.dot(pu) + r.dot(puu), mixed],
                [mixed, pv.dot(pv) + r.dot(pvv)],
            ],
            [[-pu.dot(su), -pu.dot(sv)], [-pv.dot(su), -pv.dot(sv)]],
        );
        let extremum = SurfaceSurfaceExtremum::new([u, v], point, foot.parameters, foot.point);
        (
            extremum,
            (
                0.5 * r.dot(r),
                [pu, pv].map(|along| above_rounding(r.dot(along), point, foot.point, along)),
                hessian,
            ),
            [pu.magnitude(), pv.magnitude()],
            foot,
        )
    }

    /// The outer window from `parameters` (see `outer_window`).
    fn window(&self, parameters: [f64; 2]) -> [Option<(f64, f64)>; 2] {
        let (_, speeds) = self.first.first_derivatives(parameters[0], parameters[1]);
        [0, 1].map(|i| {
            let speed = speeds[i].magnitude();
            outer_window(self.first_domains[i], parameters[i], speed, self.scale)
        })
    }

    /// The minimum reached from the seed `(parameters, start)`, when it
    /// is one over the parameter boxes.
    fn refine(&self, parameters: [f64; 2], start: [f64; 2]) -> Option<SurfaceSurfaceExtremum> {
        let mut warm = start;
        let [u, v] = newton_min_2d(
            |u, v| {
                let (extremum, objective, ..) = self.local([u, v], warm);
                warm = extremum.second_parameter_array();
                objective
            },
            parameters,
            self.window(parameters),
            1e-14,
            100,
        );
        let parameters = [self.first_domains[0].wrap(u), self.first_domains[1].wrap(v)];
        if !parameters.iter().all(|x| x.is_finite()) {
            return None;
        }
        let (extremum, (_, gradient, hessian), speeds, found) = self.local(parameters, warm);
        let residual = extremum.first_point - extremum.second_point;
        let on_first = collapsed_minimum(
            self.first,
            self.first_domains,
            parameters,
            residual,
            self.tol,
        )
        .unwrap_or_else(|| {
            let bounds = self.first_domains.map(ParameterDomain::bounds);
            is_box_minimum(parameters, gradient, hessian, speeds, bounds, self.tol)
        });
        let finite = extremum.distance.is_finite();
        let on_second = found.is_minimum(self.second, self.second_domains, self.tol);
        let holds = !found.is_degenerate() || {
            let length = probe_length(self.tol, self.scale);
            let probes = [(0, -1.0), (0, 1.0), (1, -1.0), (1, 1.0)]
                .into_iter()
                .filter_map(|(i, side)| {
                    let mut probe = parameters;
                    probe[i] += side * length / speeds[i];
                    inside(self.first_domains[i], probe[i]).map(|_| probe)
                })
                .map(|[u, v]| {
                    let point = self.first.point_at(u, v);
                    let start = extremum.second_parameter_array();
                    foot(self.second, self.second_domains, point, start)
                        .point
                        .distance(point)
                });
            holds_off_focus(extremum.distance, probes, self.tol)
        };
        (finite && on_first && on_second && holds)
            .then(|| self.closed(&extremum).unwrap_or(extremum))
    }

    /// The meeting point a minimum lies next to, when the gap between the
    /// surfaces closes from it (see `closed_gap`).
    fn closed(&self, e: &SurfaceSurfaceExtremum) -> Option<SurfaceSurfaceExtremum> {
        let ([u1, v1], [u2, v2]) = (e.first_parameter_array(), e.second_parameter_array());
        let [u1, v1, u2, v2] = closed_gap(
            [u1, v1, u2, v2],
            |[u1, v1, u2, v2]| {
                let (p, [pu, pv]) = self.first.first_derivatives(u1, v1);
                let (s, [su, sv]) = self.second.first_derivatives(u2, v2);
                (p - s, [pu, pv, -su, -sv])
            },
            self.tol,
        )?;
        let first = [
            inside(self.first_domains[0], u1)?,
            inside(self.first_domains[1], v1)?,
        ];
        let second = [
            inside(self.second_domains[0], u2)?,
            inside(self.second_domains[1], v2)?,
        ];
        Some(SurfaceSurfaceExtremum::new(
            first,
            self.first.point_at(first[0], first[1]),
            second,
            self.second.point_at(second[0], second[1]),
        ))
    }

    /// Whether two minima lie in one basin: both at distance `0`, where
    /// the surfaces meet, or walking the first surface from one to the
    /// other, each foot on the second warm-started from the last, finds no
    /// ridge (see `one_basin`) and arrives on the other's second point.
    /// Either end may start the walk, as the feet follow one sheet. Or,
    /// for a valley of the second surface over one point of the first, no
    /// ridge along the straight way in the parameters of both surfaces.
    fn same_basin(&self, p: &SurfaceSurfaceExtremum, q: &SurfaceSurfaceExtremum) -> bool {
        let confusion = self.tol.confusion;
        (p.distance <= confusion && q.distance <= confusion)
            || self.walks_to(p, q)
            || self.walks_to(q, p)
            || self.straight_between(p, q)
    }

    /// Whether no ridge parts `p` from `q` along the straight way in the
    /// parameters of both surfaces, each periodic one the shorter way round.
    fn straight_between(&self, p: &SurfaceSurfaceExtremum, q: &SurfaceSurfaceExtremum) -> bool {
        let way = |domains: [ParameterDomain; 2], from: [f64; 2], to: [f64; 2]| {
            [0, 1].map(|i| domains[i].step(from[i], to[i]))
        };
        let (first, second) = (p.first_parameter_array(), p.second_parameter_array());
        let on_first = way(self.first_domains, first, q.first_parameter_array());
        let on_second = way(self.second_domains, second, q.second_parameter_array());
        one_basin(
            p.distance,
            q.distance,
            |f| {
                let a = self
                    .first
                    .point_at(first[0] + f * on_first[0], first[1] + f * on_first[1]);
                let b = (self.second)
                    .point_at(second[0] + f * on_second[0], second[1] + f * on_second[1]);
                a.distance(b)
            },
            self.tol,
        )
    }

    /// Whether the walk from `p` to `q` stays in one basin and arrives.
    fn walks_to(&self, p: &SurfaceSurfaceExtremum, q: &SurfaceSurfaceExtremum) -> bool {
        let to = q.first_parameter_array();
        let mut feet = self.feet_from(p, to);
        let from = p.first_parameter_array();
        let way = [0, 1].map(|i| self.first_domains[i].step(from[i], to[i]));
        let no_ridge = one_basin(
            p.distance,
            q.distance,
            |f| {
                let point = self
                    .first
                    .point_at(from[0] + f * way[0], from[1] + f * way[1]);
                feet.foot_at(f).point.distance(point)
            },
            self.tol,
        );
        no_ridge && feet.foot_at(1.0).point.distance(q.second_point) <= self.tol.confusion
    }

    /// The feet on the second surface carried from `e`'s along the straight
    /// way on the first from `e`'s parameters to `to`.
    fn feet_from(
        &self,
        e: &SurfaceSurfaceExtremum,
        to: [f64; 2],
    ) -> Carried<'_, B, impl Fn(f64) -> Point + '_> {
        let from = e.first_parameter_array();
        let way = [0, 1].map(|i| self.first_domains[i].step(from[i], to[i]));
        let chord = e.first_point.distance(self.first.point_at(to[0], to[1]));
        Carried::new(
            self.second,
            self.second_domains,
            move |f| {
                self.first
                    .point_at(from[0] + f * way[0], from[1] + f * way[1])
            },
            e.second_parameter_array(),
            chord,
            self.scale,
        )
    }

    /// A valley's representative: the first surface at its reference
    /// parameters and its foot on the second, when that pair lies in the
    /// valley; else `e` with one parameter of the first surface moved to
    /// its reference, for a valley along that parameter, its foot carried
    /// along from `e`'s; else the minimum the search reaches from the
    /// reference. Along a valley the gradient is rounding alone, so where
    /// the search stops along it would hang on the rounding.
    fn representative(&self, e: &SurfaceSurfaceExtremum) -> SurfaceSurfaceExtremum {
        let reference = reference_parameters(self.first_domains);
        let point = self.first.point_at(reference[0], reference[1]);
        let projection = self.second.project_point(point, self.tol);
        let start = [projection.u, projection.v];
        let [u, v] = e.first_parameter_array();
        let in_valley = |pair: &SurfaceSurfaceExtremum| {
            pair.distance <= e.distance + self.tol.confusion && self.same_basin(e, pair)
        };
        [
            self.local(reference, start).0,
            self.carried(e, [reference[0], v]),
            self.carried(e, [u, reference[1]]),
        ]
        .into_iter()
        .find(in_valley)
        .unwrap_or_else(|| self.refine(reference, start).unwrap_or(*e))
    }

    /// The pair at the first surface's parameters `to`, its foot on the
    /// second carried from `e`'s along the straight way there, so that it
    /// stays on `e`'s sheet.
    fn carried(&self, e: &SurfaceSurfaceExtremum, to: [f64; 2]) -> SurfaceSurfaceExtremum {
        let warm = self.feet_from(e, to).foot_at(1.0).parameters;
        self.local(to, warm).0
    }
}

#[cfg(test)]
mod tests;
