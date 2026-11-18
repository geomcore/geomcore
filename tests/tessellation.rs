//! Tessellation through its public interface: on every sample surface the
//! reported deviation is checked against the distances it bounds, both
//! ways, and the triangles against the grid's rules (no triangle without
//! area, unit normals on the side the triangles turn to, every edge shared
//! by two triangles but on the boundary of the ranges).

mod common;

use common::samples::{quartic_branch, surfaces, tilted_frame};
use geomcore::prelude::*;
use geomcore::{
    BSplineSurface, Cone, Cylinder, Direction, ExtrusionSurface, Plane, Point, Point2D,
    RevolutionSurface, Sphere, Surface, Tessellation, TessellationBuffers, TessellationCriteria,
    TessellationError, Tolerance, Torus, TrimmedSurfaceConstructionError, Vector,
};
use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

const TOL: Tolerance = Tolerance::DEFAULT;

// ---------------------------------------------------------------- checker

/// The point of the triangle `[a, b, c]` closest to `p`.
fn closest_on_triangle(p: Point, [a, b, c]: [Point; 3]) -> Point {
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return a + ab * (d1 / (d1 - d3));
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return a + ac * (d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denominator = 1.0 / (va + vb + vc);
    a + ab * (vb * denominator) + ac * (vc * denominator)
}

fn corners(t: &Tessellation, triangle: [u32; 3]) -> [Point; 3] {
    triangle.map(|k| t.vertices[k as usize])
}

/// The rounding allowed on top of a deviation: a few units in the last
/// place of coordinates as large as `size`.
fn rounding(size: f64) -> f64 {
    1e-12 * size.max(1.0)
}

/// The largest coordinate of the tessellation.
fn size(t: &Tessellation) -> f64 {
    t.vertices
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
        .fold(0.0, f64::max)
}

/// Checks every rule a tessellation of `surface` keeps; `context` names
/// the case in a failure.
fn check<S: ParametricSurface>(surface: &S, t: &Tessellation, max_deviation: f64, context: &str) {
    let slack = rounding(size(t));
    assert!(
        t.deviation <= max_deviation,
        "{context}: deviation {} over {max_deviation}",
        t.deviation
    );
    let n = t.vertices.len();
    assert_eq!(t.normals.len(), n, "{context}");
    assert_eq!(t.uvs.len(), n, "{context}");
    // Unit normals; vertices on the surface at their parameters.
    for k in 0..n {
        let normal = t.normals[k].vector();
        assert!(
            (normal.magnitude() - 1.0).abs() < 1e-12,
            "{context}: normal {k}"
        );
        let (u, v) = (t.uvs[k].x, t.uvs[k].y);
        let on = surface.point_at(u, v).distance(t.vertices[k]);
        assert!(
            on <= t.deviation + slack,
            "{context}: vertex {k} off by {on}"
        );
    }
    // Triangles: valid indices, an area, turning with the normals.
    let mut around: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (index, &triangle) in t.triangles.iter().enumerate() {
        assert!(triangle.iter().all(|&k| (k as usize) < n), "{context}");
        let [a, b, c] = corners(t, triangle);
        assert!(
            a != b && b != c && c != a,
            "{context}: triangle {index} has no area"
        );
        let facing = (b - a).cross(c - a);
        for k in triangle {
            assert!(
                facing.dot(t.normals[k as usize].vector()) > 0.0,
                "{context}: triangle {index} turns against the normal of vertex {k}"
            );
            around[k as usize].push(index);
        }
    }
    // Every point of the triangles lies within the deviation of the
    // surface, and every point of the surface within it of the triangles.
    let barycentric = [
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.6, 0.2, 0.2],
        [0.2, 0.6, 0.2],
        [0.2, 0.2, 0.6],
    ];
    for &triangle in &t.triangles {
        let points = corners(t, triangle);
        let uvs = triangle.map(|k| t.uvs[k as usize]);
        let neighbours: Vec<usize> = triangle
            .iter()
            .flat_map(|&k| around[k as usize].iter().copied())
            .collect();
        for weights in barycentric {
            let mesh_point = Point::ORIGIN
                + (0..3)
                    .map(|i| (points[i] - Point::ORIGIN) * weights[i])
                    .fold(Vector::ZERO, |sum, term| sum + term);
            let (u, v) = (0..3).fold((0.0, 0.0), |(u, v), i| {
                (u + uvs[i].x * weights[i], v + uvs[i].y * weights[i])
            });
            let surface_point = surface.point_at(u, v);
            // The point at the same parameters bounds the distance; the
            // projection is searched only where that bound is not enough.
            let mut to_surface = surface_point.distance(mesh_point);
            if to_surface > t.deviation + slack {
                to_surface = to_surface.min(surface.project_point(mesh_point, TOL).distance);
            }
            assert!(
                to_surface <= t.deviation + slack,
                "{context}: a triangle point lies {to_surface} from the surface, over {}",
                t.deviation
            );
            let to_mesh = neighbours
                .iter()
                .map(|&other| {
                    let corners = corners(t, t.triangles[other]);
                    surface_point.distance(closest_on_triangle(surface_point, corners))
                })
                .fold(f64::INFINITY, f64::min);
            assert!(
                to_mesh <= t.deviation + slack,
                "{context}: a surface point lies {to_mesh} from the triangles, over {}",
                t.deviation
            );
        }
    }
}

/// Checks that every edge joins two triangles, but the edges along the
/// boundary of the ranges, which join one; edges are keyed by their end
/// points, so a seam's vertices and a collapsed row's count as one.
fn check_watertight(t: &Tessellation, (u0, u1): (f64, f64), (v0, v1): (f64, f64), context: &str) {
    let key = |p: Point| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
    let mut edges: HashMap<_, Vec<(u32, u32)>> = HashMap::new();
    for &[a, b, c] in &t.triangles {
        for (p, q) in [(a, b), (b, c), (c, a)] {
            let (kp, kq) = (key(t.vertices[p as usize]), key(t.vertices[q as usize]));
            let edge = if kp < kq { (kp, kq) } else { (kq, kp) };
            edges.entry(edge).or_default().push((p, q));
        }
    }
    let on_boundary = |k: u32| {
        let Point2D { x: u, y: v } = t.uvs[k as usize];
        u == u0 || u == u1 || v == v0 || v == v1
    };
    for uses in edges.values() {
        match uses.len() {
            2 => {}
            1 => {
                let (p, q) = uses[0];
                assert!(
                    on_boundary(p) && on_boundary(q),
                    "{context}: an inner edge has one triangle"
                );
            }
            count => panic!("{context}: an edge joins {count} triangles"),
        }
    }
}

/// The ranges of a sample: its bounds, or an interval around 0 where they
/// are infinite.
fn ranges<S: ParametricSurface>(surface: &S) -> ((f64, f64), (f64, f64)) {
    let range = |(first, last): (f64, f64)| {
        if first.is_finite() && last.is_finite() {
            (first, last)
        } else {
            (-1.5, 0.7)
        }
    };
    (range(surface.u_bounds()), range(surface.v_bounds()))
}

fn tessellated<S: ParametricSurface + Tessellate>(
    surface: &S,
    max_deviation: f64,
) -> (Tessellation, (f64, f64), (f64, f64)) {
    let (u_range, v_range) = ranges(surface);
    let criteria = TessellationCriteria::new(max_deviation)
        .unwrap()
        .with_u_range(u_range)
        .with_v_range(v_range);
    let t = surface.tessellate(criteria, TOL).unwrap();
    let size = surface.tessellation_size(criteria, TOL).unwrap();
    assert_eq!(
        (size.vertices, size.triangles),
        (t.vertices.len(), t.triangles.len())
    );
    (t, u_range, v_range)
}

// ---------------------------------------------------------------- samples

#[test]
fn every_sample_holds_its_deviation_both_ways() {
    for (index, surface) in surfaces().iter().enumerate() {
        for max_deviation in [1e-1, 1e-2] {
            let context = format!("sample {index} at {max_deviation}");
            let (t, u_range, v_range) = tessellated(surface, max_deviation);
            assert!(t.proven, "{context}");
            check(surface, &t, max_deviation, &context);
            check_watertight(&t, u_range, v_range, &context);
        }
    }
}

#[test]
fn every_trimmed_sample_holds_its_deviation_both_ways() {
    for (index, surface) in surfaces().iter().enumerate() {
        let (u_range, v_range) = ranges(surface);
        let inner = |(first, last): (f64, f64)| {
            (first + 0.15 * (last - first), first + 0.8 * (last - first))
        };
        let trimmed = surface
            .clone()
            .trimmed(inner(u_range), inner(v_range))
            .unwrap();
        let context = format!("trimmed sample {index}");
        let (t, u_range, v_range) = tessellated(&trimmed, 1e-2);
        check(&trimmed, &t, 1e-2, &context);
        check_watertight(&t, u_range, v_range, &context);
    }
}

#[test]
fn a_finer_deviation_gives_more_triangles() {
    for surface in surfaces() {
        let coarse = tessellated(&surface, 1e-1).0.triangles.len();
        let fine = tessellated(&surface, 1e-3).0.triangles.len();
        assert!(fine >= coarse);
    }
}

// ---------------------------------------------------------------- elementary

#[test]
fn a_plane_is_one_cell() {
    let plane = Plane::from_frame(tilted_frame());
    let criteria = TessellationCriteria::new(1e-6)
        .unwrap()
        .with_u_range((-1.0, 2.0))
        .with_v_range((0.0, 3.0));
    let t = plane.tessellate(criteria, TOL).unwrap();
    assert_eq!((t.vertices.len(), t.triangles.len()), (4, 2));
    assert_eq!(t.deviation, 0.0);
    check(&plane, &t, 1e-6, "plane");
}

#[test]
fn a_cylinder_steps_by_the_sagitta() {
    // A unit circle's arc of angle h is 1 - cos(h/2) from its chord: 1e-3
    // allows h = 0.0894, so 71 steps around.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let criteria = TessellationCriteria::new(1e-3)
        .unwrap()
        .with_v_range((0.0, 5.0));
    let t = cylinder.tessellate(criteria, TOL).unwrap();
    assert_eq!((t.vertices.len(), t.triangles.len()), (72 * 2, 71 * 2));
    let step = TAU / 71.0;
    assert!((t.deviation - (1.0 - (0.5 * step).cos())).abs() < 1e-15);
    check(&cylinder, &t, 1e-3, "cylinder");
    check_watertight(&t, (0.0, TAU), (0.0, 5.0), "cylinder");
}

#[test]
fn a_whole_sphere_is_closed_and_has_no_triangle_at_its_poles_without_area() {
    let sphere = Sphere::from_frame(tilted_frame(), 3.0).unwrap();
    let criteria = TessellationCriteria::new(1e-3).unwrap();
    let t = sphere.tessellate(criteria, TOL).unwrap();
    check(&sphere, &t, 1e-3, "sphere");
    // Closed: no edge on a boundary at all.
    check_watertight(&t, (f64::NAN, f64::NAN), (f64::NAN, f64::NAN), "sphere");
    // The poles' normals are the axis.
    let axis = tilted_frame().z_direction().vector();
    for (k, uv) in t.uvs.iter().enumerate() {
        if uv.y == FRAC_PI_2 {
            assert!((t.normals[k].vector() - axis).magnitude() < 1e-12);
        }
    }
}

#[test]
fn a_whole_torus_is_closed() {
    let torus = Torus::from_frame(tilted_frame(), 4.0, 1.0).unwrap();
    let t = torus
        .tessellate(TessellationCriteria::new(1e-2).unwrap(), TOL)
        .unwrap();
    check(&torus, &t, 1e-2, "torus");
    check_watertight(&t, (f64::NAN, f64::NAN), (f64::NAN, f64::NAN), "torus");
}

#[test]
fn a_spindle_torus_splits_where_its_tube_crosses_the_axis() {
    let torus = Torus::from_frame(tilted_frame(), 1.0, 2.0).unwrap();
    let t = torus
        .tessellate(TessellationCriteria::new(1e-2).unwrap(), TOL)
        .unwrap();
    check(&torus, &t, 1e-2, "spindle torus");
    check_watertight(
        &t,
        (f64::NAN, f64::NAN),
        (f64::NAN, f64::NAN),
        "spindle torus",
    );
}

#[test]
fn a_cone_across_its_apex_keeps_each_nappes_normals() {
    let cone = Cone::from_frame(tilted_frame(), 0.5, 1.0).unwrap();
    let apex = -1.0 / 0.5f64.sin();
    let criteria = TessellationCriteria::new(1e-3)
        .unwrap()
        .with_v_range((apex - 2.0, apex + 3.0));
    let t = cone.tessellate(criteria, TOL).unwrap();
    check(&cone, &t, 1e-3, "cone across its apex");
    check_watertight(
        &t,
        (0.0, TAU),
        (apex - 2.0, apex + 3.0),
        "cone across its apex",
    );
    // The apex rows: one for each nappe, their normals opposite at each u.
    let at_apex: Vec<usize> = (0..t.uvs.len()).filter(|&k| t.uvs[k].y == apex).collect();
    let columns = at_apex.len() / 2;
    assert_eq!(at_apex.len(), 2 * columns);
    for i in 0..columns {
        let (below, above) = (t.normals[at_apex[i]], t.normals[at_apex[columns + i]]);
        assert!(below.vector().dot(above.vector()) < 0.0);
    }
}

#[test]
fn a_cone_ending_at_its_apex() {
    let cone = Cone::from_frame(tilted_frame(), 0.3, 2.0).unwrap();
    let apex = -2.0 / 0.3f64.sin();
    let criteria = TessellationCriteria::new(1e-3)
        .unwrap()
        .with_v_range((apex, 0.0));
    let t = cone.tessellate(criteria, TOL).unwrap();
    check(&cone, &t, 1e-3, "cone to its apex");
    check_watertight(&t, (0.0, TAU), (apex, 0.0), "cone to its apex");
}

#[test]
fn a_range_across_the_seam_keeps_its_parameters() {
    let cylinder = Cylinder::from_frame(tilted_frame(), 2.0).unwrap();
    let criteria = TessellationCriteria::new(1e-3)
        .unwrap()
        .with_u_range((1.5 * PI, 2.5 * PI))
        .with_v_range((0.0, 1.0));
    let t = cylinder.tessellate(criteria, TOL).unwrap();
    assert!(t.uvs.iter().all(|uv| (1.5 * PI..=2.5 * PI).contains(&uv.x)));
    check(&cylinder, &t, 1e-3, "cylinder across the seam");
}

// ---------------------------------------------------------------- swept and spline

#[test]
fn a_sampled_profile_reports_an_estimated_deviation() {
    let extrusion =
        ExtrusionSurface::new(quartic_branch().into(), Vector::new(0.2, 0.1, 1.0)).unwrap();
    let criteria = TessellationCriteria::new(1e-2)
        .unwrap()
        .with_u_range((-1.5, 0.7))
        .with_v_range((0.0, 1.0));
    let t = extrusion.tessellate(criteria, TOL).unwrap();
    assert!(!t.proven);
    assert!(t.deviation <= 1e-2);
}

#[test]
fn a_revolution_of_a_skew_line_holds_its_deviation() {
    // A line skew to the axis turns into a hyperboloid: its cells twist,
    // by as much as both steps, so the turn takes steps as fine as the
    // deviation (not its square root, as a sagitta does).
    let line = geomcore::Line::new(Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 1.0)).unwrap();
    let axis = geomcore::Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let hyperboloid = RevolutionSurface::new(line.into(), axis).unwrap();
    let criteria = TessellationCriteria::new(1e-3)
        .unwrap()
        .with_v_range((-0.5, 0.5));
    let t = hyperboloid.tessellate(criteria, TOL).unwrap();
    check(&hyperboloid, &t, 1e-3, "hyperboloid");
}

#[test]
fn a_bspline_with_a_collapsed_boundary_leaves_its_triangles_without_area_out() {
    // A quadratic patch whose first row of poles is one point: a cap.
    let apex = Point::new(0.0, 0.0, 1.0);
    let poles = vec![
        vec![apex, apex, apex],
        vec![
            Point::new(-1.0, 0.0, 0.5),
            Point::new(0.0, 1.0, 0.6),
            Point::new(1.0, 0.0, 0.5),
        ],
        vec![
            Point::new(-2.0, 0.0, 0.0),
            Point::new(0.0, 2.0, 0.2),
            Point::new(2.0, 0.0, 0.0),
        ],
    ];
    let knots = (vec![0.0, 1.0], vec![3, 3]);
    let surface = BSplineSurface::new(
        2,
        2,
        poles,
        knots.0.clone(),
        knots.1.clone(),
        knots.0,
        knots.1,
        false,
        false,
    )
    .unwrap();
    let t = surface
        .tessellate(TessellationCriteria::new(1e-3).unwrap(), TOL)
        .unwrap();
    check(&surface, &t, 1e-3, "B-spline cap");
    check_watertight(&t, (0.0, 1.0), (0.0, 1.0), "B-spline cap");
}

// ---------------------------------------------------------------- buffers and errors

#[test]
fn the_buffers_take_what_the_size_says() {
    let sphere: Surface = Sphere::new(Point::ORIGIN, 1.0).unwrap().into();
    let criteria = TessellationCriteria::new(1e-2).unwrap();
    let size = sphere.tessellation_size(criteria, TOL).unwrap();
    let mut vertices = vec![Point::ORIGIN; size.vertices];
    let mut normals = vec![Direction::Z; size.vertices];
    let mut uvs = vec![Point2D::ORIGIN; size.vertices];
    let mut triangles = vec![[0; 3]; size.triangles];
    let reached = sphere
        .tessellate_into(
            criteria,
            TOL,
            TessellationBuffers {
                vertices: &mut vertices,
                normals: &mut normals,
                uvs: &mut uvs,
                triangles: &mut triangles,
            },
        )
        .unwrap();
    let t = sphere.tessellate(criteria, TOL).unwrap();
    assert_eq!(t.vertices, vertices);
    assert_eq!(t.normals, normals);
    assert_eq!(t.uvs, uvs);
    assert_eq!(t.triangles, triangles);
    assert_eq!((t.deviation, t.proven), (reached.deviation, reached.proven));
}

#[test]
#[should_panic(expected = "triangles: wrong length")]
fn a_buffer_of_the_wrong_length_panics() {
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let criteria = TessellationCriteria::new(1e-2).unwrap();
    let size = sphere.tessellation_size(criteria, TOL).unwrap();
    let mut vertices = vec![Point::ORIGIN; size.vertices];
    let mut normals = vec![Direction::Z; size.vertices];
    let mut uvs = vec![Point2D::ORIGIN; size.vertices];
    let mut triangles = vec![[0; 3]; size.triangles + 1];
    let _ = sphere.tessellate_into(
        criteria,
        TOL,
        TessellationBuffers {
            vertices: &mut vertices,
            normals: &mut normals,
            uvs: &mut uvs,
            triangles: &mut triangles,
        },
    );
}

#[test]
fn the_errors() {
    assert_eq!(
        TessellationCriteria::new(0.0),
        Err(TessellationError::InvalidDeviation)
    );
    assert_eq!(
        TessellationCriteria::new(f64::NAN),
        Err(TessellationError::InvalidDeviation)
    );
    let criteria = TessellationCriteria::new(1e-3).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(
        plane.tessellation_size(criteria, TOL),
        Err(TessellationError::UnboundedU)
    );
    assert_eq!(
        plane.tessellation_size(criteria.with_u_range((0.0, 1.0)), TOL),
        Err(TessellationError::UnboundedV)
    );
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    assert_eq!(
        sphere.tessellation_size(criteria.with_v_range((0.0, 2.0)), TOL),
        Err(TessellationError::InvalidVRange(
            TrimmedSurfaceConstructionError::OutsideBasis
        ))
    );
    assert_eq!(
        sphere.tessellation_size(criteria.with_u_range((0.0, 7.0)), TOL),
        Err(TessellationError::InvalidURange(
            TrimmedSurfaceConstructionError::ExceedsPeriod
        ))
    );
    assert_eq!(
        sphere.tessellation_size(criteria.with_u_range((1.0, 1.0)), TOL),
        Err(TessellationError::InvalidURange(
            TrimmedSurfaceConstructionError::EmptyRange
        ))
    );
    // A trimmed surface's ranges bound it.
    let band = sphere.trimmed((0.0, 1.0), (0.0, 1.0)).unwrap();
    assert_eq!(
        band.tessellation_size(criteria.with_u_range((0.5, 1.5)), TOL),
        Err(TessellationError::InvalidURange(
            TrimmedSurfaceConstructionError::OutsideBasis
        ))
    );
    assert!(
        band.tessellation_size(criteria.with_u_range((0.5, 0.75)), TOL)
            .is_ok()
    );
    // A deviation far below the size asks for more vertices than a u32 reaches.
    assert_eq!(
        sphere.tessellation_size(TessellationCriteria::new(1e-300).unwrap(), TOL),
        Err(TessellationError::TooManyVertices)
    );
}
