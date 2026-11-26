//! Arcs a few tolerances across near a cone (#516), graded against the
//! distance to the cone worked out here from the shapes' own data, not by
//! the crate's evaluation or distance routines.

use crate::geometry::curves::ParametricCurve;
use crate::{
    Circle, Cone, CurveSurfaceIntersection, Ellipse, Frame, Intersect, Point, Tolerance,
    TrimmedCurve, Vector,
};
use std::f64::consts::TAU;

const TOL: Tolerance = Tolerance::DEFAULT;

/// Samples of the reference distance along a range.
const SAMPLES: usize = 20_000;

/// How far beyond the tolerance a reported value may sit for rounding:
/// the shapes are placed about 100 units from the origin, where a
/// coordinate rounds at `1.4e-14`.
const ROUNDING: f64 = 1e-5;

/// A distance clearly within the tolerance, which the answer must hold a
/// point for; between it and the tolerance rounding may decide either way.
const CLEARLY_WITHIN: f64 = 0.995;

/// `a + b` as an unevaluated sum `(sum, error)`, exactly.
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let sum = a + b;
    let b_virtual = sum - a;
    let a_virtual = sum - b_virtual;
    (sum, (a - a_virtual) + (b - b_virtual))
}

/// A circle or an ellipse `center + a cos t x + b sin t y` against the
/// double cone of apex `apex`, unit `axis` and half-angle `angle`.
///
/// The arc's points are taken relative to the apex: the center less the
/// apex is kept exactly as a sum of two doubles, so an arc `1e-7` across
/// placed 100 units out keeps every digit of its offset from the cone.
struct Reference {
    from_apex: [(f64, f64); 3],
    x: [f64; 3],
    y: [f64; 3],
    a: f64,
    b: f64,
    axis: [f64; 3],
    cosine: f64,
    sine: f64,
    apex: [f64; 3],
}

impl Reference {
    fn new(frame: Frame, a: f64, b: f64, cone: &Cone) -> Self {
        let v = |d: Vector| [d.x, d.y, d.z];
        let (center, apex) = (frame.origin(), cone.apex());
        let apex = [apex.x, apex.y, apex.z];
        let center = [center.x, center.y, center.z];
        let (sine, cosine) = cone.half_angle().sin_cos();
        Self {
            from_apex: [0, 1, 2].map(|i| two_sum(center[i], -apex[i])),
            x: v(frame.x_direction().vector()),
            y: v(frame.y_direction().vector()),
            a,
            b,
            axis: v(cone.frame().z_direction().vector()),
            cosine,
            sine,
            apex,
        }
    }

    fn circle(circle: &Circle, cone: &Cone) -> Self {
        Self::new(circle.frame(), circle.radius(), circle.radius(), cone)
    }

    fn ellipse(ellipse: &Ellipse, cone: &Cone) -> Self {
        Self::new(
            ellipse.frame(),
            ellipse.major_radius(),
            ellipse.minor_radius(),
            cone,
        )
    }

    /// The arc's point at `t` less the apex.
    fn offset(&self, t: f64) -> [f64; 3] {
        let (sine, cosine) = t.sin_cos();
        [0, 1, 2].map(|i| {
            let (high, low) = self.from_apex[i];
            high + (low + (self.a * cosine * self.x[i] + self.b * sine * self.y[i]))
        })
    }

    /// `point` less the apex.
    fn point_offset(&self, point: Point) -> [f64; 3] {
        let point = [point.x, point.y, point.z];
        [0, 1, 2].map(|i| {
            let (high, low) = two_sum(point[i], -self.apex[i]);
            high + low
        })
    }

    /// The signed distance from the point `apex + w` to the double cone,
    /// `ρ cos φ − |h| sin φ`: exact for a double cone, whose nearest
    /// generator is never more than a right angle away.
    fn signed(&self, w: [f64; 3]) -> f64 {
        let height = w[0] * self.axis[0] + w[1] * self.axis[1] + w[2] * self.axis[2];
        let radial = [0, 1, 2].map(|i| w[i] - height * self.axis[i]);
        let rho = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
        rho * self.cosine - height.abs() * self.sine
    }

    fn distance_at(&self, t: f64) -> f64 {
        self.signed(self.offset(t))
    }

    /// How far the arc's point at `t` lies from `point`.
    fn gap(&self, t: f64, point: Point) -> f64 {
        let (on, off) = (self.offset(t), self.point_offset(point));
        let d = [0, 1, 2].map(|i| on[i] - off[i]);
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }
}

/// Checks the answer for the arc over `(low, high)` against the reference
/// distance sampled along it:
///
/// - each point's parameter lies in the range, the arc there within the
///   tolerance of the cone and of the point;
/// - no parameter twice, and no two points joined by a stretch clearly
///   within the tolerance: those are one contact (ADR 0034);
/// - each stretch clearly within the tolerance, a sign change among them,
///   holds a point;
/// - an overlap only where the arc is within the tolerance all along.
fn check(
    reference: &Reference,
    result: &CurveSurfaceIntersection,
    (low, high): (f64, f64),
    context: &str,
) {
    let band = TOL.confusion;
    let limit = band * (1.0 + ROUNDING);
    let parameter = |k: usize| low + (high - low) * k as f64 / SAMPLES as f64;
    let values: Vec<f64> = (0..=SAMPLES)
        .map(|k| reference.distance_at(parameter(k)))
        .collect();
    let hits = result.hits();
    let found: Vec<(f64, u32)> = hits
        .iter()
        .map(|h| (h.curve_parameter, h.multiplicity))
        .collect();
    if !result.overlaps().is_empty() {
        let worst = values.iter().fold(0f64, |m, v| m.max(v.abs()));
        assert!(
            worst <= limit,
            "{context}: overlap, the arc reaches {worst}"
        );
        return;
    }
    for h in hits {
        let t = h.curve_parameter;
        assert!(
            (low..=high).contains(&t),
            "{context}: point at {t} outside [{low}, {high}], {found:?}"
        );
        let distance = reference.distance_at(t).abs();
        assert!(
            distance <= limit,
            "{context}: the arc at {t} is {} tol from the cone, {found:?}",
            distance / band
        );
        let gap = reference.gap(t, h.point);
        assert!(
            gap <= limit,
            "{context}: the point at {t} is {} tol off the arc, {found:?}",
            gap / band
        );
    }
    let mut parameters: Vec<f64> = hits.iter().map(|h| h.curve_parameter).collect();
    parameters.sort_by(f64::total_cmp);
    for pair in parameters.windows(2) {
        assert!(pair[0] < pair[1], "{context}: {} twice, {found:?}", pair[0]);
        let between = values
            .iter()
            .enumerate()
            .filter(|&(k, _)| (pair[0]..=pair[1]).contains(&parameter(k)))
            .fold(0f64, |m, (_, v)| m.max(v.abs()));
        assert!(
            between > CLEARLY_WITHIN * band,
            "{context}: points at {} and {} joined within the tolerance, {found:?}",
            pair[0],
            pair[1]
        );
    }
    // The stretches within the tolerance, each a run of samples
    // `start..end`, with the nearest it comes and whether it changes sign.
    let mut runs = Vec::new();
    let mut k = 0;
    while k <= SAMPLES {
        if values[k].abs() > band {
            k += 1;
            continue;
        }
        let start = k;
        let (mut nearest, mut sign_change) = (f64::INFINITY, false);
        while k <= SAMPLES && values[k].abs() <= band {
            nearest = nearest.min(values[k].abs());
            sign_change |= k > start && values[k] * values[k - 1] < 0.0;
            k += 1;
        }
        runs.push((start, k, nearest, sign_change));
    }
    let holds = |&(start, end, _, _): &(usize, usize, f64, bool), t: f64| {
        (parameter(start.saturating_sub(1))..=parameter(end.min(SAMPLES))).contains(&t)
    };
    // Over a whole turn the range's ends are one point: a run through it
    // is one stretch.
    let whole_turn = high - low >= TAU * (1.0 - 1e-12);
    let joined = whole_turn
        && runs.len() > 1
        && runs[0].0 == 0
        && runs.last().is_some_and(|run| run.1 > SAMPLES);
    for (index, run) in runs.iter().enumerate() {
        let partner = match (joined, index) {
            (true, 0) => runs.last(),
            (true, i) if i == runs.len() - 1 => runs.first(),
            _ => None,
        };
        let nearest = partner.map_or(run.2, |other| run.2.min(other.2));
        let sign_change = run.3 || partner.is_some_and(|other| other.3);
        let inside = |t: f64| holds(run, t) || partner.is_some_and(|other| holds(other, t));
        if nearest <= CLEARLY_WITHIN * band || sign_change {
            assert!(
                parameters.iter().any(|&t| inside(t)),
                "{context}: nothing in run {:?}, nearest {} tol, {found:?}",
                (parameter(run.0), parameter(run.1.min(SAMPLES))),
                nearest / band
            );
        }
        // A run inside the range crosses the cone when its ends lie on
        // either side of it, and touches it otherwise.
        let (start, end) = (run.0, run.1);
        let in_run: Vec<u32> = hits
            .iter()
            .filter(|h| inside(h.curve_parameter))
            .map(|h| h.multiplicity)
            .collect();
        if start > 0 && end <= SAMPLES && in_run.len() == 1 {
            let (first, last) = (values[start], values[end - 1]);
            if first.abs().min(last.abs()) > 0.5 * band {
                let crosses = first.signum() != last.signum();
                assert_eq!(
                    in_run[0],
                    if crosses { 1 } else { 2 },
                    "{context}: multiplicity in run {:?}, {found:?}",
                    (parameter(start), parameter(end))
                );
            }
        }
    }
}

fn frame(o: [f64; 3], x: [f64; 3], y: [f64; 3]) -> Frame {
    let v = |c: [f64; 3]| Vector::new(c[0], c[1], c[2]);
    let z = [
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    ];
    Frame::from_axes(Point::new(o[0], o[1], o[2]), v(x), v(y), v(z)).unwrap()
}

/// A curve of center `o`, axes `x`, `y` and radii, a cone of frame
/// `(co, cx, cy, cz)`, half-angle `angle` and apex at `co`, and a range.
type Case = (
    ([f64; 3], [f64; 3], [f64; 3], &'static [f64]),
    ([f64; 3], [f64; 3], [f64; 3], [f64; 3], f64),
    (f64, f64),
);

/// The answer for a case of #516, checked by [`check`].
fn solve(((o, x, y, radii), (co, cx, cy, cz, angle), range): Case, name: &str) -> Vec<(f64, u32)> {
    let curve_frame = frame(o, x, y);
    let v = |c: [f64; 3]| Vector::new(c[0], c[1], c[2]);
    let cone_frame =
        Frame::from_axes(Point::new(co[0], co[1], co[2]), v(cx), v(cy), v(cz)).unwrap();
    let cone = Cone::from_frame(cone_frame, angle, 0.0).unwrap();
    let (reference, result) = if let [radius] = *radii {
        let circle = Circle::from_frame(curve_frame, radius).unwrap();
        let arc = TrimmedCurve::new(circle, range.0, range.1).unwrap();
        (Reference::circle(&circle, &cone), arc.intersect(&cone, TOL))
    } else {
        let ellipse = Ellipse::from_frame(curve_frame, radii[0], radii[1]).unwrap();
        let arc = TrimmedCurve::new(ellipse, range.0, range.1).unwrap();
        (
            Reference::ellipse(&ellipse, &cone),
            arc.intersect(&cone, TOL),
        )
    };
    assert!(result.is_complete(), "{name}");
    check(&reference, &result, range, name);
    result
        .hits()
        .iter()
        .map(|h| (h.curve_parameter, h.multiplicity))
        .collect()
}

/// Whether one of `found` lies within `gap` of `parameter`.
fn has_near(found: &[(f64, u32)], parameter: f64, gap: f64) -> bool {
    found.iter().any(|&(t, _)| (t - parameter).abs() <= gap)
}

/// A full circle of radius `4.6e-8` comes within 0.763 tolerances of the
/// cone at 4.4611 and is 1.065 away at 0: the contact is at 4.4611, not
/// snapped to the seam.
#[test]
fn test_a_dip_of_a_tiny_circle_keeps_its_parameter() {
    let found = solve(
        (
            (
                [-89.46658225129956, 89.12630495042073, -75.55804831524729],
                [-0.41341057110554047, 0.7608706381800793, 0.5001675435827776],
                [-0.748785866294052, -0.5966300931440257, 0.2887079118992365],
                &[4.61006204769539e-08],
            ),
            (
                [-89.46658234507585, 89.12630524693964, -75.55804866350746],
                [
                    -0.8231947691980463,
                    0.35278093748861056,
                    0.44485501245870257,
                ],
                [
                    -0.5536825716785023,
                    -0.6722068409033569,
                    -0.4915013457379434,
                ],
                [0.12564227705854258, -0.6509098041874821, 0.7486858119586336],
                0.34114123384775524,
            ),
            (0.0, TAU),
        ),
        "dip missed",
    );
    assert_eq!(found.len(), 1, "{found:?}");
    assert!(has_near(&found, 4.4611, 0.05), "{found:?}");
}

/// An ellipse coming within 0.982 tolerances at 5.7633 and 1.012 away at
/// 2π touches at 5.7633.
#[test]
fn test_a_dip_of_a_tiny_ellipse_keeps_its_parameter() {
    let found = solve(
        (
            (
                [-92.48995297057903, 99.59448370620004, -10.53640886248977],
                [
                    0.35012282425458247,
                    -0.4552675298200658,
                    -0.8186241409997206,
                ],
                [
                    -0.058332576237364535,
                    0.8616477854063403,
                    -0.5041432380324676,
                ],
                &[3.228660451335302e-08, 1.5622376868941468e-08],
            ),
            (
                [-92.48995284260262, 99.5944836612388, -10.536408852457484],
                [0.4020695488111843, 0.1988899465408449, -0.893746534026147],
                [
                    -0.8472242445931784,
                    0.45096680701477376,
                    -0.2807846476294924,
                ],
                [0.3472047771737727, 0.8700986887335712, 0.34982440534049547],
                0.4214648832933235,
            ),
            (3.246440154447605, 6.9927630566223575),
        ),
        "ellipse dip missed",
    );
    assert_eq!(found.len(), 1, "{found:?}");
    assert!(has_near(&found, 5.7633, 0.05), "{found:?}");
}

/// The arc is within the tolerance only at its two ends, 0.778 and 0.895
/// tolerances from the cone, and 1.191 away at its nearest inside: a
/// point at each end, joined only outside the range.
#[test]
fn test_a_tiny_arc_within_the_tolerance_only_at_its_ends_has_a_point_at_each() {
    let range = (4.449137490072062, 8.900353306761264);
    let found = solve(
        (
            (
                [-86.82018361701176, -23.31432763903799, 28.445888462052682],
                [0.18591460040104832, 0.0915726627682782, -0.9782894299701116],
                [
                    -0.8211347096854305,
                    -0.532309154986639,
                    -0.20587557423656264,
                ],
                &[3.1941602501763135e-08],
            ),
            (
                [-86.82018333950091, -23.314327838575693, 28.445888325620245],
                [
                    -0.24005290813756042,
                    0.19544561027809124,
                    -0.950881493519421,
                ],
                [
                    -0.5793737536921014,
                    -0.8147860378147721,
                    -0.02120768999275998,
                ],
                [-0.7789099144493055, 0.5458248125590567, 0.30882781475677756],
                0.3557493799286992,
            ),
            range,
        ),
        "point beyond tol inside the range",
    );
    assert_eq!(found.len(), 2, "{found:?}");
    assert!(has_near(&found, range.0, 1e-3), "{found:?}");
    assert!(has_near(&found, range.1, 1e-3), "{found:?}");
}

/// Crossings at 5.6477 and 6.2742 are not both reported at 2π.
#[test]
fn test_two_crossings_of_a_tiny_arc_do_not_share_a_parameter() {
    solve(
        (
            (
                [95.98954985657625, -9.490615073881084, -41.796108374620786],
                [0.7656457043397121, -0.5443283378232753, -0.3427729803655715],
                [
                    -0.4682398800845533,
                    -0.10623791852365298,
                    -0.8771914952655235,
                ],
                &[1.4787741032793394e-07],
            ),
            (
                [95.98954934789944, -9.490615728742398, -41.796107264569486],
                [
                    0.7248975367237498,
                    -0.6776851445924829,
                    -0.12355729865331921,
                ],
                [
                    -0.5981389528437625,
                    -0.5302552369499964,
                    -0.6008819990464604,
                ],
                [0.3416918997325034, 0.5094823142031222, -0.7897305978442422],
                0.1440482556363024,
            ),
            (4.58081127651212, 9.497667993364324),
        ),
        "same parameter twice",
    );
}

/// Within the tolerance only at the range's start, 0.514 tolerances from
/// the cone: the point is the arc's own there, not one from outside the
/// range.
#[test]
fn test_a_tiny_arc_touching_at_its_start_reports_its_own_point() {
    let range = (0.5452076360168637, 1.8495507328218532);
    let found = solve(
        (
            (
                [-28.632644171743802, 10.695273488523924, 91.15290315177226],
                [0.39249788586055107, 0.4111306035330702, 0.8227496802998657],
                [0.4556790960443249, 0.6901083973210353, -0.5622339027266311],
                &[1.3110703191722558e-07],
            ),
            (
                [-28.632643775408397, 10.695273134721447, 91.15290302310751],
                [0.6603948067239213, 0.7477562546393164, -0.0688424498391548],
                [
                    -0.1320776495411709,
                    0.025418023417337184,
                    -0.9909134263785285,
                ],
                [-0.739211873378502, 0.6634866296567984, 0.11554781920446668],
                0.3822050050352508,
            ),
            range,
        ),
        "point off",
    );
    assert_eq!(found.len(), 1, "{found:?}");
    assert!(has_near(&found, range.0, 1e-3), "{found:?}");
}

/// Crossings at 2.2377 and 2.6952, 0.711 tolerances apart with the arc
/// within the tolerance between them, are one contact (ADR 0034).
#[test]
fn test_crossings_joined_within_the_tolerance_are_one_contact() {
    let found = solve(
        (
            (
                [-35.5050204793883, 5.142934443148261, 74.8863518147076],
                [0.5807861887105095, -0.668775064267898, 0.46414148318868204],
                [
                    -0.35309185904187235,
                    0.30677923594790435,
                    0.8838623419229799,
                ],
                &[1.5678848204451716e-07],
            ),
            (
                [-35.50502047305211, 5.142934751698775, 74.88635182341255],
                [-0.680616239243626, 0.3169286822118132, -0.6605435226155404],
                [0.5464603823714951, -0.3809309044757732, -0.7458369101309407],
                [
                    -0.48799855058828073,
                    -0.8685895788039095,
                    0.08607879075001221,
                ],
                0.15178361808809515,
            ),
            (0.9476204055231849, 4.594057080021418),
        ),
        "two crossings closer than tol",
    );
    assert_eq!(found, [(found[0].0, 2)], "{found:?}");
    assert!(found[0].0 > 2.2 && found[0].0 < 2.75, "{found:?}");
}

/// A uniform sample in `[0, 1)` from a fixed xorshift sequence.
fn sequence() -> impl FnMut() -> f64 {
    let mut seed = 0x2545_F491_4F6C_DD1Du64;
    move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A random unit vector.
fn unit(random: &mut impl FnMut() -> f64) -> Vector {
    loop {
        let v = Vector::new(
            2.0 * random() - 1.0,
            2.0 * random() - 1.0,
            2.0 * random() - 1.0,
        );
        let length = v.magnitude();
        if length > 0.1 && length <= 1.0 {
            return v / length;
        }
    }
}

/// Circles and ellipses a few tolerances across, cut from cones placed
/// at random up to 100 units out, each of shift, size and tilt moved by
/// up to two tolerances, trimmed to random ranges: see [`check`].
#[test]
fn test_tiny_arcs_near_a_cone_agree_with_the_reference_distance() {
    let mut random = sequence();
    let band = TOL.confusion;
    for case in 0..2000 {
        let angle = 0.1 + 1.3 * random();
        let height = 10f64.powf(-7.0 + 1.5 * random()) * if random() < 0.5 { -1.0 } else { 1.0 };
        let axis = unit(&mut random);
        let side = unit(&mut random).cross(axis);
        let Some(first) = side.normalized() else {
            continue;
        };
        let second = axis.cross(first);
        let apex = Point::new(
            200.0 * random() - 100.0,
            200.0 * random() - 100.0,
            200.0 * random() - 100.0,
        );
        let Ok(cone_frame) = Frame::from_axes(apex, first, second, axis) else {
            continue;
        };
        let cone = Cone::from_frame(cone_frame, angle, 0.0).unwrap();
        let jitter = |r: f64| (2.0 * r - 1.0) * 2.0 * band;
        // In the cone's frame, the plane `z = h + x tan β` cuts the ellipse
        // of the near-elliptic sweep of the parent module; `β = 0` a circle.
        let ellipse = case % 2 == 1;
        let tilt = if ellipse {
            0.9 * (std::f64::consts::FRAC_PI_2 - angle) * random()
        } else {
            0.0
        };
        let (sine, cosine) = tilt.sin_cos();
        let k = angle.tan().powi(2);
        let a = cosine * cosine - k * sine * sine;
        let offset = height * k * sine / a;
        let minor = (k * height * height + a * offset * offset).sqrt() + jitter(random());
        let major = if ellipse { minor / a.sqrt() } else { minor };
        if minor <= 0.05 * band {
            continue;
        }
        let local = |x: f64, y: f64, z: f64| first * x + second * y + axis * z;
        let center = apex
            + local(
                offset * cosine + jitter(random()),
                jitter(random()),
                height + offset * sine + jitter(random()),
            );
        let lean = (2.0 * band / minor).min(0.3);
        let normal = local(
            -sine + lean * (2.0 * random() - 1.0),
            lean * (2.0 * random() - 1.0),
            cosine,
        );
        let x_direction = local(cosine, 0.0, sine);
        let low = TAU * random();
        let high = if random() < 0.2 {
            low + TAU
        } else {
            low + (0.05 + 0.95 * random()) * TAU
        };
        let context = format!("arc {case}: angle {angle}, height {height}, range ({low}, {high})");
        if ellipse {
            let Ok(curve) = Ellipse::new(center, normal, x_direction, major, minor) else {
                continue;
            };
            let arc = TrimmedCurve::new(curve, low, high).unwrap();
            let result = arc.intersect(&cone, TOL);
            assert!(result.is_complete(), "{context}");
            check(
                &Reference::ellipse(&curve, &cone),
                &result,
                (low, high),
                &context,
            );
        } else {
            let Ok(curve) = Circle::new(center, normal, minor) else {
                continue;
            };
            let arc = TrimmedCurve::new(curve, low, high).unwrap();
            let result = arc.intersect(&cone, TOL);
            assert!(result.is_complete(), "{context}");
            check(
                &Reference::circle(&curve, &cone),
                &result,
                (low, high),
                &context,
            );
        }
    }
}

/// The reference agrees with the arc's own points far from any cone
/// rounding question: a guard on [`Reference`] itself.
#[test]
fn test_the_reference_follows_the_arc() {
    let circle = Circle::from_frame(
        frame([3.0, -2.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        0.5,
    )
    .unwrap();
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 0.0).unwrap();
    let reference = Reference::circle(&circle, &cone);
    for t in [0.0, 1.0, 2.5, 4.0] {
        assert!(reference.gap(t, circle.point_at(t)) < 1e-14);
    }
}
