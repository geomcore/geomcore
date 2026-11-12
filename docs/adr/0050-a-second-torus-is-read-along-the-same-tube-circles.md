# A second torus is read along the same tube circles

Two tori met in closed form only with collinear axes, and every other position was marched (#272): two pipe bends meeting at a corner, a fillet against a bend, two rings linked like a chain. The documentation called the general curve "up to degree 16". But a torus is bicircular, and that bound is not reached. On the first torus's tube circles, `w = P − C₂` gives `|w|² = |D|² + R² + r² + 2Rr cos v + 2(R + r cos v) D·e(u) + 2r sin v D·a₁`, with `D = C₁ − C₂`. The parts of `(R + r cos v)² + r² sin² v` of degree two cancel, so `|w|²` is of degree one in each angle, and so is `w·a₂`. The second torus's equation `(|w|² + R₂² − r₂²)² − 4R₂²(|w|² − (w·a₂)²)` is then of degree two in each, like a quadric's. So each tube circle meets the second torus in up to four points, and the curve is of degree eight.

So the sweep of ADR 0049 reads it unchanged. The other surface became an enum of a cylinder, a cone or a torus. Each gives its equation `F(u, v)`, its signed distance, a bound on `|F|` within confusion from the largest gradient of its equation, and its moves for snapping. Like a cylinder's radius, a torus's minor radius is offset exactly (an offset torus is a torus), so a torus touching within confusion is moved by one translation or one change of minor radius. The fit reads the second torus's angles as angles on both sides.

**A shared tube circle.** Pipe bends joined end to end share a tube circle: the minor radii are equal, and the centre circles touch at its centre. The tori touch all along it, and `F` vanishes to second order on the line `u = u₀`. The event search then meets a curve of common zeros and gives up. But each coefficient of `v`, of degree two in `u` with a double zero at `u₀`, factors as `(1 − cos(u − u₀))(p + q cos(u − u₀) + s sin(u − u₀))`. The quotient is again a trigonometric polynomial, of degree one in `u`, and its zeros are the rest of the curve. So the circle is found and divided out, and the rest is swept as before:

- **Finding it.** The circle's centre lies on the first centre circle where `e(u₀)` is across `(C₁ − C₂) × a₁` and `a₂ × a₁`, giving at most four candidates. The circle is accepted where sixteen of its points lie within confusion of the second torus, with parallel normals.
- **Touching points.** Points where the rest touches are read from the quotient alone, within the bound `F` has across from the circle, where `1 − cos(u − u₀)` is two. They are not moved, since the tori touch along the circle anyway, and points on the circle are dropped from `contact_points`.
- **The fit.** Where the rest of the curve crosses the circle the normals are parallel, so every loop's tangent there is stepped, as through a crossing.

The circle comes back from `OcticCurve::shared_circle` and is routed as a tangent curve. With no loops left, it is a `TangentCircle`.

A circle of one torus can lie on the other in two more ways: as a latitude circle of the other, or as a Villarceau circle of either. On a latitude circle the tori cross, and `F` vanishes only to first order on `u = u₀`, which no trigonometric polynomial factor takes out. A Villarceau circle is not a line of either's parameters. Both stay `NotAnalytic` and are marched.

**Distance near crossings.** The loops through a crossing slid to every distance (ADR 0049), at three evaluations of the curve per step. Two tori crossing where they touch made this the dearest call in the fit: side by side, 157M instructions. Away from the crossing the surfaces cross at an angle and the first order is right, so `CurveFamily::slides_along` now takes the point. The octic family slides only where the normals are within a sine of 0.25 of parallel, which brought the call down to 74M and the pipe bends from 46M to 15M. The quadrics' families ignore the point.

## Considered options

- **A resultant of degree 16**, the bound the documentation gave: rejected. The equation along a tube circle is a quartic, as for a quadric, and the sweep already reads it.
- **Divide out every shared circle**: deferred. A latitude circle leaves a factor that is not a trigonometric polynomial, and a Villarceau circle is a curve in `(u, v)`. Each needs its own quotient, and pipe bends are the case parts meet.
- **Slide everywhere on a loop through a crossing**: kept as the fallback near crossings only, for the cost above.
