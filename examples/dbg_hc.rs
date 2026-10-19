use geomcore::IntersectExact;
use geomcore::{Cone, Frame, Hyperbola, Point, Tolerance, Vector};
fn main() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    match h.intersect_exact(&cone, tol) {
        geomcore::ConicSurfaceIntersection::Hits(hits) => {
            for h in &hits {
                println!(
                    "t={:.6} mult={} pt={:?} cone-resid={:.3e}",
                    h.parameter,
                    h.multiplicity,
                    h.point,
                    cone.parameters_of(h.point).0
                );
                let (u, v) = cone.parameters_of(h.point);
                let back = cone.eval_point(u, v);
                println!("   back={back:?} resid={:.3e}", back.distance(h.point));
            }
        }
        other => println!("{other:?}"),
    }
}
