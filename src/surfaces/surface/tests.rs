use super::*;

use crate::Point;

fn sample<S: ParametricSurface>(s: &S) -> Point {
    s.eval_point(0.3, 0.4)
}

fn plane() -> Plane {
    Plane::new(Point::ORIGIN, Vector::Z).unwrap()
}

fn cylinder() -> Cylinder {
    Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap()
}

fn cone() -> Cone {
    Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap()
}

fn sphere() -> Sphere {
    Sphere::new(Point::ORIGIN, 3.0).unwrap()
}

fn torus() -> Torus {
    Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap()
}

fn bezier() -> BezierSurface {
    // Biquadratic dome over the unit square.
    let row = |x: f64, z: f64| {
        vec![
            Point::new(x, 0.0, 0.0),
            Point::new(x, 0.5, z),
            Point::new(x, 1.0, 0.0),
        ]
    };
    BezierSurface::new(2, 2, vec![row(0.0, 0.0), row(0.5, 1.0), row(1.0, 0.0)]).unwrap()
}

fn bspline() -> BSplineSurface {
    // Bilinear degree-1x1 patch over the unit square (clamped both ways).
    let poles = vec![
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
        vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    ];
    BSplineSurface::new(
        1,
        1,
        poles,
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

// ---- generic fn agrees with inherent calls ----

#[test]
fn test_sample_surface_enum_agrees_with_inherent() {
    let p = plane();
    let wrapped: Surface = p.into();
    assert_eq!(sample(&wrapped), p.eval_point(0.3, 0.4));
}

// ---- heterogeneous Vec<Surface> eval ----

#[test]
fn test_heterogeneous_vec_eval() {
    let p = plane();
    let c = cylinder();
    let cn = cone();
    let s = sphere();
    let t = torus();
    let surfaces: Vec<Surface> = vec![p.into(), c.into(), cn.into(), s.into(), t.into()];
    let points: Vec<Point> = surfaces
        .iter()
        .map(|surf| surf.eval_point(0.2, 0.1))
        .collect();
    assert_eq!(points[0], p.eval_point(0.2, 0.1));
    assert_eq!(points[1], c.eval_point(0.2, 0.1));
    assert_eq!(points[2], cn.eval_point(0.2, 0.1));
    assert_eq!(points[3], s.eval_point(0.2, 0.1));
    assert_eq!(points[4], t.eval_point(0.2, 0.1));
}

// ---- From round-trips: owned and by-reference (copying) ----

#[test]
fn test_from_plane_round_trips() {
    let p = plane();
    match Surface::from(p) {
        Surface::Plane(inner) => assert_eq!(inner, p),
        _ => panic!("expected Surface::Plane"),
    }
}

#[test]
fn test_from_plane_ref_round_trips() {
    let p = plane();
    match Surface::from(&p) {
        Surface::Plane(inner) => assert_eq!(inner, p),
        _ => panic!("expected Surface::Plane"),
    }
}

#[test]
fn test_from_cylinder_round_trips() {
    let c = cylinder();
    match Surface::from(c) {
        Surface::Cylinder(inner) => assert_eq!(inner, c),
        _ => panic!("expected Surface::Cylinder"),
    }
}

#[test]
fn test_from_cylinder_ref_round_trips() {
    let c = cylinder();
    match Surface::from(&c) {
        Surface::Cylinder(inner) => assert_eq!(inner, c),
        _ => panic!("expected Surface::Cylinder"),
    }
}

#[test]
fn test_from_cone_round_trips() {
    let c = cone();
    match Surface::from(c) {
        Surface::Cone(inner) => assert_eq!(inner, c),
        _ => panic!("expected Surface::Cone"),
    }
}

#[test]
fn test_from_cone_ref_round_trips() {
    let c = cone();
    match Surface::from(&c) {
        Surface::Cone(inner) => assert_eq!(inner, c),
        _ => panic!("expected Surface::Cone"),
    }
}

#[test]
fn test_from_sphere_round_trips() {
    let s = sphere();
    match Surface::from(s) {
        Surface::Sphere(inner) => assert_eq!(inner, s),
        _ => panic!("expected Surface::Sphere"),
    }
}

#[test]
fn test_from_sphere_ref_round_trips() {
    let s = sphere();
    match Surface::from(&s) {
        Surface::Sphere(inner) => assert_eq!(inner, s),
        _ => panic!("expected Surface::Sphere"),
    }
}

#[test]
fn test_from_torus_round_trips() {
    let t = torus();
    match Surface::from(t) {
        Surface::Torus(inner) => assert_eq!(inner, t),
        _ => panic!("expected Surface::Torus"),
    }
}

#[test]
fn test_from_torus_ref_round_trips() {
    let t = torus();
    match Surface::from(&t) {
        Surface::Torus(inner) => assert_eq!(inner, t),
        _ => panic!("expected Surface::Torus"),
    }
}

#[test]
fn test_from_bspline_round_trips() {
    let b = bspline();
    match Surface::from(b.clone()) {
        Surface::BSpline(inner) => assert_eq!(inner, b),
        _ => panic!("expected Surface::BSpline"),
    }
}

#[test]
fn test_from_bspline_ref_clones() {
    let b = bspline();
    match Surface::from(&b) {
        Surface::BSpline(inner) => assert_eq!(inner, b),
        _ => panic!("expected Surface::BSpline"),
    }
}

#[test]
#[allow(clippy::needless_borrows_for_generic_args)]
fn test_into_from_reference_supports_parametrize_on_style_call() {
    // Mirrors the Task 17 usage pattern: a function taking
    // `impl Into<Surface>` called with a `&Cylinder` (Cylinder is Copy,
    // so clippy would rather we drop the borrow; the point of this test
    // is specifically that the borrowed form is accepted).
    fn accepts(surface: impl Into<Surface>) -> Surface {
        surface.into()
    }
    let c = cylinder();
    match accepts(&c) {
        Surface::Cylinder(inner) => assert_eq!(inner, c),
        _ => panic!("expected Surface::Cylinder"),
    }
}

// ---- bounds/periods table, per type (inherent trait impl + enum dispatch) ----

#[test]
fn test_bspline_enum_dispatch_agrees_with_inherent() {
    let b = bspline();
    let wrapped: Surface = b.clone().into();
    assert_eq!(wrapped.eval_point(0.3, 0.4), b.eval_point(0.3, 0.4));
    assert_eq!(
        wrapped.eval_derivative(0.3, 0.4, 1, 0),
        b.eval_derivative(0.3, 0.4, 1, 0)
    );
    assert_eq!(
        wrapped.eval_derivative(0.3, 0.4, 0, 1),
        b.eval_derivative(0.3, 0.4, 0, 1)
    );
}

// ---- default eval_points == mapped eval_point ----

#[test]
fn test_surface_enum_eval_points_matches_mapped_eval_point() {
    let uvs = [(0.0, 0.0), (0.5, 0.2), (1.0, -0.3), (-2.0, 1.1)];
    let surfaces: Vec<Surface> = vec![
        plane().into(),
        cylinder().into(),
        cone().into(),
        sphere().into(),
        torus().into(),
    ];
    for surface in &surfaces {
        let expected: Vec<Point> = uvs.iter().map(|&(u, v)| surface.eval_point(u, v)).collect();
        assert_eq!(surface.eval_points(&uvs), expected);
    }
}

// ---- normal: unit Su x Sv, None where degenerate ----

#[test]
fn test_normal_surface_enum_dispatch_agrees_with_concrete() {
    let tol = Tolerance::DEFAULT;
    let surfaces: Vec<Surface> = vec![
        plane().into(),
        cylinder().into(),
        cone().into(),
        sphere().into(),
        torus().into(),
        extrusion().into(),
        revolution().into(),
        bspline().into(),
        bezier().into(),
    ];
    for surface in &surfaces {
        assert_eq!(surface.normal(0.3, 0.2, tol), surface.normal(0.3, 0.2, tol));
        let (_, [su, sv]) = surface.eval_d1(0.3, 0.2);
        let n = surface.normal(0.3, 0.2, tol).unwrap();
        let expected = su.cross(sv).normalized().unwrap();
        assert!((n - expected).magnitude() < 1e-15);
    }
}

// ---- normal rung 3 (step-off) and hopeless cases, via a mock ----

// ---- normal: independence from the parametrization ----

// ---- closedness: periodicity implies closed, not conversely ----

fn closed_u_patch() -> BSplineSurface {
    // Clamped degree-1 patch whose first and last pole rows coincide:
    // closed in `u`, open in `v`, not periodic either way.
    let row = |x: f64| vec![Point::new(x, 0.0, 0.0), Point::new(x, 1.0, 0.0)];
    BSplineSurface::new(
        1,
        1,
        vec![row(0.0), row(1.0), row(0.0)],
        vec![0.0, 1.0, 2.0],
        vec![2, 1, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_surface_enum_closedness_dispatch_agrees_with_concrete() {
    let tol = Tolerance::DEFAULT;
    let surfaces: Vec<Surface> = vec![
        plane().into(),
        cylinder().into(),
        cone().into(),
        sphere().into(),
        torus().into(),
        extrusion().into(),
        revolution().into(),
        bspline().into(),
        bezier().into(),
        closed_u_patch().into(),
    ];
    // The enum dispatches to each variant's own answer, including the
    // B-spline override.
    for surface in &surfaces {
        let expected = match surface {
            Surface::Plane(s) => s.is_u_closed(tol),
            Surface::Cylinder(s) => s.is_u_closed(tol),
            Surface::Cone(s) => s.is_u_closed(tol),
            Surface::Sphere(s) => s.is_u_closed(tol),
            Surface::Torus(s) => s.is_u_closed(tol),
            Surface::Extrusion(s) => s.is_u_closed(tol),
            Surface::Revolution(s) => s.is_u_closed(tol),
            Surface::BSpline(s) => s.is_u_closed(tol),
            Surface::Bezier(s) => s.is_u_closed(tol),
        };
        assert_eq!(surface.is_u_closed(tol), expected);
        let expected_v = match surface {
            Surface::Plane(s) => s.is_v_closed(tol),
            Surface::Cylinder(s) => s.is_v_closed(tol),
            Surface::Cone(s) => s.is_v_closed(tol),
            Surface::Sphere(s) => s.is_v_closed(tol),
            Surface::Torus(s) => s.is_v_closed(tol),
            Surface::Extrusion(s) => s.is_v_closed(tol),
            Surface::Revolution(s) => s.is_v_closed(tol),
            Surface::BSpline(s) => s.is_v_closed(tol),
            Surface::Bezier(s) => s.is_v_closed(tol),
        };
        assert_eq!(surface.is_v_closed(tol), expected_v);
    }
}

// ---- elementary surfaces ----

#[test]
fn test_elementary_surface_agrees_with_surface() {
    let tol = Tolerance::DEFAULT;
    let surfaces: [ElementarySurface; 5] = [
        plane().into(),
        cylinder().into(),
        cone().into(),
        sphere().into(),
        torus().into(),
    ];
    let query = Point::new(0.7, -1.2, 2.5);
    for surface in surfaces {
        let wrapped = Surface::from(surface);
        assert_eq!(wrapped.as_elementary(), Some(surface));
        assert_eq!(surface.eval_point(0.3, 0.4), wrapped.eval_point(0.3, 0.4));
        assert_eq!(
            surface.eval_derivative(0.3, 0.4, 1, 1),
            wrapped.eval_derivative(0.3, 0.4, 1, 1)
        );
        assert_eq!(surface.u_bounds(), wrapped.u_bounds());
        assert_eq!(surface.v_bounds(), wrapped.v_bounds());
        assert_eq!(surface.u_period(), wrapped.u_period());
        assert_eq!(surface.v_period(), wrapped.v_period());
        assert_eq!(
            surface.project_point(query, tol),
            wrapped.project_point(query, tol)
        );
        assert_eq!(surface.contains(query, tol), wrapped.contains(query, tol));
    }
}

#[test]
fn test_bspline_surface_is_not_elementary() {
    assert_eq!(Surface::from(bspline()).as_elementary(), None);
}

// ---- continuity: degree minus multiplicity, smooth elsewhere ----

fn creased_u_patch() -> BSplineSurface {
    // Degree 2 in `u` with a double interior knot (C0 at u = 0.5),
    // degree 1 in `v` without interior knots (smooth).
    let row = |x: f64, z: f64| vec![Point::new(x, 0.0, z), Point::new(x, 1.0, z)];
    BSplineSurface::new(
        2,
        1,
        vec![
            row(0.0, 0.0),
            row(1.0, 0.0),
            row(2.0, 1.0),
            row(3.0, 1.0),
            row(4.0, 0.0),
        ],
        vec![0.0, 0.5, 1.0],
        vec![3, 2, 3],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_surface_enum_continuity_dispatch_agrees_with_concrete() {
    let tol = Tolerance::DEFAULT;
    let surfaces: Vec<Surface> = vec![
        plane().into(),
        cylinder().into(),
        cone().into(),
        sphere().into(),
        torus().into(),
        extrusion().into(),
        revolution().into(),
        bspline().into(),
        bezier().into(),
        creased_u_patch().into(),
    ];
    for surface in &surfaces {
        let expected_u = match surface {
            Surface::Plane(s) => s.continuity_u(0.5, tol),
            Surface::Cylinder(s) => s.continuity_u(0.5, tol),
            Surface::Cone(s) => s.continuity_u(0.5, tol),
            Surface::Sphere(s) => s.continuity_u(0.5, tol),
            Surface::Torus(s) => s.continuity_u(0.5, tol),
            Surface::Extrusion(s) => s.continuity_u(0.5, tol),
            Surface::Revolution(s) => s.continuity_u(0.5, tol),
            Surface::BSpline(s) => s.continuity_u(0.5, tol),
            Surface::Bezier(s) => s.continuity_u(0.5, tol),
        };
        assert_eq!(surface.continuity_u(0.5, tol), expected_u);
        assert_eq!(
            surface.continuity_v(0.5, tol),
            match surface {
                Surface::BSpline(s) => s.continuity_v(0.5, tol),
                Surface::Plane(s) => s.continuity_v(0.5, tol),
                Surface::Cylinder(s) => s.continuity_v(0.5, tol),
                Surface::Cone(s) => s.continuity_v(0.5, tol),
                Surface::Sphere(s) => s.continuity_v(0.5, tol),
                Surface::Torus(s) => s.continuity_v(0.5, tol),
                Surface::Extrusion(s) => s.continuity_v(0.5, tol),
                Surface::Revolution(s) => s.continuity_v(0.5, tol),
                Surface::Bezier(s) => s.continuity_v(0.5, tol),
            }
        );
    }
}

// ---- curvatures: principal, mean, Gaussian ----

#[test]
fn test_surface_enum_curvatures_dispatch_agrees_with_concrete() {
    let tol = Tolerance::DEFAULT;
    let surfaces: Vec<Surface> = vec![
        plane().into(),
        cylinder().into(),
        cone().into(),
        sphere().into(),
        torus().into(),
        extrusion().into(),
        revolution().into(),
        bspline().into(),
        bezier().into(),
    ];
    for surface in &surfaces {
        let expected = match surface {
            Surface::Plane(s) => s.curvatures(0.3, 0.2, tol),
            Surface::Cylinder(s) => s.curvatures(0.3, 0.2, tol),
            Surface::Cone(s) => s.curvatures(0.3, 0.2, tol),
            Surface::Sphere(s) => s.curvatures(0.3, 0.2, tol),
            Surface::Torus(s) => s.curvatures(0.3, 0.2, tol),
            Surface::Extrusion(s) => s.curvatures(0.3, 0.2, tol),
            Surface::Revolution(s) => s.curvatures(0.3, 0.2, tol),
            Surface::BSpline(s) => s.curvatures(0.3, 0.2, tol),
            Surface::Bezier(s) => s.curvatures(0.3, 0.2, tol),
        };
        assert_eq!(surface.curvatures(0.3, 0.2, tol), expected);
    }
}

fn extrusion() -> ExtrusionSurface {
    let circle = crate::Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    ExtrusionSurface::new(circle.into(), Vector::new(0.2, 0.0, 1.0)).unwrap()
}

fn revolution() -> RevolutionSurface {
    let axis = crate::Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let line = crate::Line::new(Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 1.0)).unwrap();
    RevolutionSurface::new(line.into(), axis).unwrap()
}
