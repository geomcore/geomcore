//! Bounding boxes allocate nothing, so culling many pairs costs no heap
//! traffic. A test binary of its own, as it replaces the global allocator.

use geomcore::{
    BSplineCurve, BSplineSurface, Circle, Curve, Ellipse, HasBoundingBox, Hyperbola, Parabola,
    Point, Sphere, Surface, Torus, TrimmedCurve, Vector,
};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

struct CountingAllocator;

thread_local! {
    static ALLOCATIONS: Cell<usize> = const { Cell::new(0) };
}

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCATIONS.with(|count| count.set(count.get() + 1));
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
}

#[global_allocator]
static GLOBAL: CountingAllocator = CountingAllocator;

fn assert_boxes_without_allocating(shape: &impl HasBoundingBox) {
    let before = ALLOCATIONS.with(Cell::get);
    let b = shape.bounding_box();
    assert_eq!(ALLOCATIONS.with(Cell::get) - before, 0);
    assert!(b.is_some());
}

#[test]
fn test_boxes_do_not_allocate() {
    let n = Vector::new(0.3, -0.5, 0.8);
    let x = n.cross(Vector::X);
    let circle = Circle::new(Point::ORIGIN, n, 2.0).unwrap();
    let ellipse = Ellipse::new(Point::ORIGIN, n, x, 3.0, 1.0).unwrap();
    let parabola = Parabola::new(Point::ORIGIN, n, x, 0.5).unwrap();
    let hyperbola = Hyperbola::new(Point::ORIGIN, n, x, 1.0, 2.0).unwrap();
    let poles = (0..6)
        .map(|i| {
            let a = f64::from(i);
            Point::new(a, a.sin(), 0.5 * a.cos())
        })
        .collect();
    let weights = vec![1.0, 0.5, 2.0, 1.0, 3.0, 1.0];
    let spline = BSplineCurve::new_rational(
        3,
        poles,
        weights,
        vec![0.0, 1.0, 2.0, 3.0],
        vec![4, 1, 1, 4],
        false,
    )
    .unwrap();
    let net = vec![
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 1.0)],
        vec![Point::new(1.0, 0.0, 2.0), Point::new(1.0, 1.0, 0.0)],
    ];
    let surface = BSplineSurface::new(
        1,
        1,
        net,
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap();

    assert_boxes_without_allocating(&circle);
    assert_boxes_without_allocating(&ellipse);
    assert_boxes_without_allocating(&TrimmedCurve::new(circle, 5.0, 7.0).unwrap());
    assert_boxes_without_allocating(&TrimmedCurve::new(ellipse, 1.0, 3.0).unwrap());
    assert_boxes_without_allocating(&TrimmedCurve::new(parabola, -2.0, 3.0).unwrap());
    assert_boxes_without_allocating(&TrimmedCurve::new(hyperbola, -1.0, 2.0).unwrap());
    assert_boxes_without_allocating(
        &TrimmedCurve::segment(Point::ORIGIN, Point::ORIGIN + n).unwrap(),
    );
    assert_boxes_without_allocating(&spline);
    assert_boxes_without_allocating(&TrimmedCurve::new(&spline, 0.5, 2.5).unwrap());
    let run_time = TrimmedCurve::new(Curve::from(spline.clone()), 0.5, 2.5).unwrap();
    assert_boxes_without_allocating(&run_time);
    assert_boxes_without_allocating(&Sphere::new(Point::ORIGIN, 1.0).unwrap());
    assert_boxes_without_allocating(&Torus::new(Point::ORIGIN, n, 3.0, 1.0).unwrap());
    assert_boxes_without_allocating(&surface);
    assert_boxes_without_allocating(&Surface::from(surface.clone()));
}
