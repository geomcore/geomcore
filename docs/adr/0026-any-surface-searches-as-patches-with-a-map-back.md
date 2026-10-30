# Any surface searches as patches with a map back

ADR 0007 converts an elementary shape to its exact rational B-spline so that a curve against a B-spline goes through the proven subdivision search. This decision applies the same conversion to surface–surface pairs, and puts it in one module, `intersect/numeric/surface_patches.rs`.

`SurfacePatches` holds a surface as Bézier-searchable B-spline patches and maps a point of a patch back to the surface's own parameters:

- a B-spline is its own patch, borrowed, with nothing to map;
- an elementary surface is converted as before (a cone holding its apex split there into its two nappes; a sphere or torus whole);
- an extrusion of a bounded profile is converted with the profile's B-spline across and the sweep along, over the part inside a cover box; a revolution of a bounded profile is converted whole.

An unbounded surface (a plane, cylinder, cone or extrusion) needs a cover box, which every hit lies inside: against a B-spline curve, the box of its control points grown by a margin; against another surface, the other's box grown the same way. The cover is a closure, called only when the surface needs one, so a bounded surface never computes it. Two unbounded surfaces have no cover and are not converted. The map back reads a straight direction as is and an angle (a circle or ellipse profile, the turn of a revolution) through the curve's own parameter, nearest the point's own position on the patch so that it stays on its own side of a seam, then settled into the surface's range. The hits and regions of a curve search are carried back through the same map (`carry_hits`).

Its consumers:

- **The curve searches** (`pairs_with_bspline.rs`): a B-spline curve against any surface that converts goes through the hull search. That adds the swept surfaces, which were sampled; a B-spline curve against an extrusion or a revolution of a bounded profile now reports complete or incomplete like the other pairs of ADR 0007. The conversion is exact, so the completeness it proves holds for the surface itself.
- **The patch-pair search** (`patch_pairs.rs`) takes two `SurfacePatches` and searches every pair of their patches. Two B-splines run as before, one patch each.
- **The seeding policy of marching** (ADR 0025) appends the patch-pair seeds after the iso-curve or grid seeds of every other pair that converts, so that a loop smaller than the spacing of the iso-curves and grids is found: a bump of a B-spline finer than its grid against a plane, cylinder or extrusion, a small cap a plane cuts off a torus. The policy leaves the search out in two cases where it found nothing new in the probes and cost the most:
  - both surfaces bounded and neither a B-spline (two spheres or tori, a sphere and a torus): both parameter grids already sample each surface at its own scale, and every small loop probed was found by them; the search cost up to 133 % more on `sphere_torus`;
  - every iso-curve or grid seed lies where the surfaces only touch: the search halves every pair along the touch down to its smallest pieces (623 % more on a plane tangent to a cylinder patch, 168 % on a plane tangent to a torus patch).

  The search still costs where it runs: on the benches, up to twice the seeding of a plane cutting a torus obliquely, 30 % on a cylinder against an extruded wall, 14 % to 29 % on the plane cuts of B-spline patches of swept surfaces. Those are the classes where it now finds the loops the tests show. A patch pair whose boxes do not overlap is dropped at once, so surfaces that miss stay cheap.

Marching stays sampled (ADR 0007): the patch-pair seeds add loops, but the iso-curve and grid seeds before them are still samples, and a patch-pair search does not prove that each curve it seeds was traced whole.

The curve side (`HullCurve`, `hull_curve`) stays in `pairs_with_bspline.rs`. It maps one parameter, restricts a trimmed curve to the part among the other's control points, and serves the curve–curve searches as well, which never meet a surface; folding it into `SurfacePatches` would make one module of two halves with no shared caller. It shares the angle mapping (`nearest`, `settle`) and the cover helpers (`grown`, `reach`), which move to `surface_patches.rs`, and its map back of hits becomes one method like the surface's.

## Considered options

- **Patch-pair seeds for every pair that converts**: rejected. On two bounded elementary surfaces and on surfaces that only touch, it found no loop the grids missed, at the costs above.
- **Converting the swept surfaces through the routing's operand view**: rejected. The searches are generic over `ParametricSurface`, and the swept surfaces join them through two hidden hooks, `as_extrusion` and `as_revolution`, like `as_bspline` and `as_elementary` (ADR 0022).
- **Converting a swept surface over an unbounded profile** (a line, parabola or hyperbola) over the cover too: left out. It needs the profile trimmed to the cover first, and no pair in the tests or fixtures has one; such a pair keeps its sampled search.
- **Computing the cover eagerly**: rejected. A B-spline curve against a B-spline surface, the most common pair, would compute and discard a box on every call.
