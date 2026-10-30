use super::*;
use crate::curves::Line;
use crate::surfaces::{Cone, Cylinder, Plane, Sphere, Torus};
use crate::{Axis, Ellipse};
use std::f64::consts::FRAC_PI_4;

const TOL: Tolerance = Tolerance::DEFAULT;

/// The box from `-half` to `half` along each axis.
fn cube(half: f64) -> Aabb {
    [-half, -half, -half, half, half, half]
}

/// Asserts that every point of a 9 by 9 grid over each patch of `patches`
/// is the point of `surface` at the parameters it maps back to, and that
/// those lie within the surface's bounds (a periodic parameter within one
/// period from its start).
fn assert_maps_back<S: ParametricSurface>(surface: &S, patches: &SurfacePatches) {
    let within = |x: f64, (lo, hi): (f64, f64), period: Option<f64>| match period {
        Some(_) => (lo..hi).contains(&x),
        None => (lo..=hi).contains(&x),
    };
    for patch in patches.patches() {
        let ((u0, u1), (v0, v1)) = (patch.u_bounds(), patch.v_bounds());
        for i in 0..=8 {
            for j in 0..=8 {
                let u = u0 + (u1 - u0) * f64::from(i) / 8.0;
                let v = v0 + (v1 - v0) * f64::from(j) / 8.0;
                let (su, sv) = patches.parameters((u, v));
                let (on_patch, on_surface) = (patch.eval_point(u, v), surface.eval_point(su, sv));
                assert!(
                    on_patch.distance(on_surface) < 1e-9,
                    "({u}, {v}) -> ({su}, {sv}): {on_patch:?} vs {on_surface:?}"
                );
                assert!(within(su, surface.u_bounds(), surface.u_period()), "u {su}");
                assert!(within(sv, surface.v_bounds(), surface.v_period()), "v {sv}");
            }
        }
    }
}

/// A cubic B-spline curve through the plane `y = 0`, away from the z axis.
fn meridian() -> Curve {
    let poles = [(1.0, -1.0), (2.0, -0.5), (1.5, 0.5), (2.5, 1.0)]
        .map(|(x, z)| Point::new(x, 0.0, z))
        .to_vec();
    BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false)
        .unwrap()
        .into()
}

#[test]
fn test_a_bspline_is_its_own_patch() {
    let bspline = Plane::new(Point::ORIGIN, Vector::Z)
        .unwrap()
        .to_bspline_patch((0.0, 2.0), (-1.0, 1.0))
        .unwrap();
    let patches = SurfacePatches::new(&bspline, || -> Option<Aabb> {
        unreachable!("a B-spline needs no cover")
    })
    .unwrap();
    let only: Vec<_> = patches.patches().collect();
    assert_eq!(only.len(), 1);
    assert!(std::ptr::eq(only[0], &bspline));
    assert_eq!(patches.parameters((0.3, 0.7)), (0.3, 0.7));
    assert_eq!(patches.carried_range(1, (-0.5, 0.25)), (-0.5, 0.25));
}

#[test]
fn test_every_elementary_surface_maps_back() {
    let surfaces: [ElementarySurface; 5] = [
        Plane::new(Point::new(0.1, 0.2, 0.3), Vector::new(0.2, 0.1, 1.0))
            .unwrap()
            .into(),
        Cylinder::new(Point::ORIGIN, Vector::new(0.0, 0.3, 1.0), 1.5)
            .unwrap()
            .into(),
        // The apex at `z = -4 / tan(0.3)`, outside the cover.
        Cone::new(Point::ORIGIN, Vector::Z, 0.3, 4.0)
            .unwrap()
            .into(),
        Sphere::new(Point::new(0.5, 0.0, 0.0), 2.0).unwrap().into(),
        Torus::new(Point::ORIGIN, Vector::X, 2.0, 0.5)
            .unwrap()
            .into(),
    ];
    for surface in surfaces {
        let patches = SurfacePatches::new(&surface, || Some(cube(3.0))).unwrap();
        assert_eq!(patches.patches().count(), 1);
        assert_maps_back(&surface, &patches);
    }
}

#[test]
fn test_only_unbounded_surfaces_need_a_cover() {
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 2.0, 1.0).unwrap();
    assert!(SurfacePatches::new(&sphere, || None).is_some());
    assert!(SurfacePatches::new(&torus, || None).is_some());
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    assert!(SurfacePatches::new(&plane, || None).is_none());
    assert!(SurfacePatches::new(&cylinder, || None).is_none());
}

#[test]
fn test_a_cone_holding_its_apex_splits_at_it() {
    // The apex at the origin: `v = -1 / sin(PI / 4)`.
    let cone = Cone::new(Point::new(0.0, 0.0, 1.0), Vector::Z, FRAC_PI_4, 1.0).unwrap();
    let apex = -std::f64::consts::SQRT_2;
    let patches = SurfacePatches::new(&cone, || Some(cube(2.0))).unwrap();
    let nappes: Vec<_> = patches.patches().collect();
    assert_eq!(nappes.len(), 2);
    assert!((nappes[0].v_bounds().1 - apex).abs() < 1e-12);
    assert!((nappes[1].v_bounds().0 - apex).abs() < 1e-12);
    assert_maps_back(&cone, &patches);
    // Every point of the shared boundary is the apex, and maps to it.
    for i in 0..=8 {
        let u = std::f64::consts::TAU * f64::from(i) / 8.0;
        let (su, sv) = patches.parameters((u, apex));
        assert!((sv - apex).abs() < 1e-12);
        assert!(cone.eval_point(su, sv).distance(Point::ORIGIN) < 1e-12);
    }
    // A cover above the apex leaves one nappe.
    let above = [-2.0, -2.0, 0.5, 2.0, 2.0, 2.0];
    let patches = SurfacePatches::new(&cone, || Some(above)).unwrap();
    assert_eq!(patches.patches().count(), 1);
    assert_maps_back(&cone, &patches);
}

#[test]
fn test_angles_come_back_on_their_own_branch() {
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let patches = SurfacePatches::new(&sphere, || None).unwrap();
    // The patch's latitude follows the sphere's only at its knots: a
    // southern latitude comes back negative, not a turn up.
    for latitude in [-1.5, -0.7, -0.1] {
        let (_, v) = patches.parameters((1.0, latitude));
        assert!(v < 0.0 && (v - latitude).abs() < 0.1, "{latitude} -> {v}");
    }
    // The seam: the end of the turn reads as its start, just before it
    // stays just before it.
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 2.0, 1.0).unwrap();
    let patches = SurfacePatches::new(&torus, || None).unwrap();
    let (u, v) = patches.parameters((TAU, TAU));
    assert!(u.abs() < 1e-12 && v.abs() < 1e-12, "({u}, {v})");
    let (u, _) = patches.parameters((TAU - 1e-6, 1.0));
    assert!(TAU - 1e-5 < u && u < TAU, "{u}");
    // A range ending on the seam keeps its end there rather than flipping
    // to the start.
    let (lo, hi) = patches.carried_range(0, (6.0, TAU));
    assert!(lo < hi && (hi - TAU).abs() < 1e-12, "({lo}, {hi})");
}

#[test]
fn test_swept_surfaces_map_back() {
    let circle = Circle::new(Point::new(0.5, 0.0, 0.0), Vector::Z, 1.0).unwrap();
    let ellipse = Ellipse::new(Point::new(3.0, 0.0, 0.0), Vector::Y, Vector::X, 1.0, 0.5).unwrap();
    let oblique = Vector::new(0.3, 0.0, 1.0);
    let extrusions = [
        ExtrusionSurface::new(circle.into(), oblique).unwrap(),
        ExtrusionSurface::new(meridian(), oblique).unwrap(),
    ];
    for extrusion in &extrusions {
        assert!(SurfacePatches::new(extrusion, || None).is_none());
        let patches = SurfacePatches::new(extrusion, || Some(cube(3.0))).unwrap();
        assert_maps_back(extrusion, &patches);
        // Every point of the surface inside the cover lies on the patch.
        // The sweep runs at `rise` to the z axis.
        let rise = 1.0 / oblique.magnitude();
        let (v0, v1) = patches.patches().next().unwrap().v_bounds();
        assert!(v0 <= -3.0 / rise && 3.0 / rise <= v1, "({v0}, {v1})");
    }
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    for profile in [ellipse.into(), meridian()] {
        let revolution = RevolutionSurface::new(profile, axis).unwrap();
        let patches = SurfacePatches::new(&revolution, || None).unwrap();
        assert_maps_back(&revolution, &patches);
    }
}

#[test]
fn test_unbounded_profiles_do_not_convert() {
    let line = Line::new(Point::new(1.0, 0.0, 0.0), Vector::new(0.2, 0.0, 1.0)).unwrap();
    let extrusion = ExtrusionSurface::new(line.into(), Vector::Y).unwrap();
    assert!(SurfacePatches::new(&extrusion, || Some(cube(1.0))).is_none());
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let revolution = RevolutionSurface::new(line.into(), axis).unwrap();
    assert!(SurfacePatches::new(&revolution, || Some(cube(1.0))).is_none());
}

#[test]
fn test_an_unbounded_surface_is_covered_by_the_other() {
    let sphere = Sphere::new(Point::new(5.0, 0.0, 0.0), 1.0).unwrap();
    let plane = Plane::new(Point::new(5.0, 0.0, 0.0), Vector::new(0.0, 0.3, 1.0)).unwrap();
    let (on_plane, on_sphere) = pair_patches(&plane, &sphere, TOL).unwrap();
    let (cover, held) = (on_plane.bounds(), on_sphere.bounds());
    // The plane's patch reaches past the sphere's box on every side the
    // plane runs along.
    for c in 0..2 {
        assert!(cover[c] < held[c] && held[c + 3] < cover[c + 3]);
    }
    assert_maps_back(&plane, &on_plane);
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    assert!(pair_patches(&plane, &cylinder, TOL).is_none());
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 2.0, 1.0).unwrap();
    let (whole, _) = pair_patches(&torus, &sphere, TOL).unwrap();
    assert_maps_back(&torus, &whole);
}
