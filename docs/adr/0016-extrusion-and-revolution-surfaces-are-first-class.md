# Extrusion and revolution surfaces are first-class; offset is an operation

`ExtrusionSurface` and `RevolutionSurface` are surface types of their own. They are not constructors that return a `BSplineSurface`. Each holds its **profile**, a `Curve`, and `Surface` gains `Extrusion` and `Revolution` variants.

- An extrusion sweeps the profile along a unit direction `d`: `S(u, v) = C(u) + v d`. `u` is the profile's parameter and `v` the signed distance travelled, so a circle extruded along its normal is parametrized as a `Cylinder` is.
- A revolution turns the profile about an `Axis`: `S(u, v)` is `C(v)` rotated by the angle `u`, and `u` is periodic with period `2π`. `u` goes around the axis and `v` along the meridian, as on a `Torus`, a `Sphere`, a `Cone` and a `Cylinder`. STEP's `surface_of_revolution` puts the angle in `v`, so an importer swaps the two parameters.

A surface offset is an operation, `offset(distance)`, not a stored type. The offset of an elementary surface is elementary, so the operation returns it exactly:

- a plane moves along its normal;
- a cylinder or sphere changes its radius by `distance`;
- a torus changes its minor radius by `distance`;
- a cone keeps its half-angle and moves its apex along the axis.

The sign follows the surface normal, and an offset that would leave a radius at or below zero is an error. An offset of a B-spline, extrusion or revolution surface reports that it is not representable until an offset surface type exists. That type waits for a caller that needs it.

## The structure is what makes them fast

Both surfaces convert exactly to rational B-splines, so a converting constructor would lose no points. It would lose the structure, and the structure turns their queries into smaller problems:

- **Projection onto an extrusion is a 2D projection.** For a fixed `u` the closest `v` is `(p - C(u)) · d`, which leaves a point-curve projection in the plane perpendicular to `d`. A parallel projection is affine, so the profile's image in that plane is a curve of the same kind: a line, a conic, or a B-spline with the same knots and weights. The 2D tier (ADR 0014) already projects onto each of them, in closed form for lines and conics.
- **Projection onto a revolution is a one-parameter search.** For a fixed `v`, the closest `u` is the angle of the point about the axis. That leaves the distance from the point's `(radius, height)` about the axis to the profile's meridian curve `v ↦ (radius(v), height(v))`. When the profile lies in a plane through the axis, this curve is the profile itself in 2D, and the 2D tier projects onto it. Otherwise it is a one-parameter Newton search, which replaces the two-parameter search on a surface.
- **Common intersections have closed forms.**
  - A plane perpendicular to the revolution axis cuts the surface in the circles through the profile's crossings with that plane.
  - A plane through the axis cuts a profile that lies in a plane through the axis in two rotated copies of it.
  - A plane parallel to the extrusion direction cuts the surface in lines through the profile's crossings.
  - Any other plane cuts an extrusion in the profile carried along `d` onto that plane, which is again a line, a conic or a B-spline.

  ADR 0011 routes these pairs exact; converted to B-splines, they would get the numeric search only.
- **An extrusion has no B-spline form.** It is unbounded in `v`, as a `Cylinder` is, so like the elementary surfaces it converts one patch at a time with `to_bspline_patch`.

Some profiles sweep an elementary surface:

- a line extruded across its own direction sweeps a plane;
- a circle extruded along its normal sweeps a cylinder;
- a line revolved in a plane through the axis sweeps a cylinder, a cone or a plane;
- a circle revolved in a plane through the axis sweeps a torus, or a sphere when its centre is on the axis.

The swept type keeps the caller's profile and parametrization, which is what a file or a modelling operation wrote. A query on such a surface runs on the equivalent elementary surface and maps the parameters back with an affine change in each direction: a shift, a scale and possibly a reversal.

## Considered options

- **Constructors that convert to `BSplineSurface`**: rejected. They lose the dimension-reducing projections and the exact intersections above, and an extrusion is unbounded, so it has no single B-spline form.
- **Construction that returns the elementary surface when one is swept**: rejected. The caller's parametrization would change under them, and a profile read from a file would not round-trip. The reduction is internal, behind the parameter maps.
- **A profile generic over its curve type (`ExtrusionSurface<C = Curve>`)**: rejected for now. `Surface` needs one concrete type, the Python classes need one, and evaluating the profile costs far more than the branch that picks its variant. Static dispatch can come later if a benchmark shows that branch matters.
- **A general sweep of a profile along a path curve**: not now. It has no exact representation, and extrusion and revolution are the two sweeps that STEP files and turned and prismatic parts use.
- **An offset surface type now**: deferred.
  - The offset of a freeform surface is not rational in general, so it has no exact B-spline conversion.
  - Its first derivative needs the base surface's second derivatives, and its second derivative needs the base's third, which surfaces do not evaluate today.
  - It folds into cusps and self-intersections wherever the distance exceeds a radius of curvature, and every numeric search over it has to cope with that.

  None of this is needed while no caller offsets a freeform surface. When one does, the offset of a revolution is the revolution of the profile's offset in its meridian plane, and the offset of an extrusion is the extrusion of the profile's offset across `d`. Both reduce to offset curves, which belong to the same decision.

## Consequences

Not every surface is elementary or a B-spline any more: an extrusion or revolution surface is a third kind, and the glossary says so.

`Surface` gains two variants that hold a `Curve`, so every match on `Surface` grows by two arms. Each new type also multiplies the pairs the intersection router has to handle. Through `ParametricSurface`, every pair with a swept surface gets the numeric search from the first version, so no pair is blocked; the exact pairs above are added one at a time, where a benchmark shows they pay.

A profile that sweeps no area is rejected at construction: a line extruded along its own direction, or a circle revolved about its own axis. A profile that crosses the revolution axis is accepted. It sweeps a surface with a singular point on the axis, as a cone has at its apex, and the normal rules of ADR 0010 apply there.
