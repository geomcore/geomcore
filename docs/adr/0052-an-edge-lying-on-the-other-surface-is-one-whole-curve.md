# An edge lying on the other surface is one whole curve

Where an edge of one surface (the bound of a bounded parameter, an iso-curve) lies on the other surface, the march reports that edge as one curve, from end to end, and traces only the rest of the intersection. The edge is not cut where other curves of the intersection meet it, and it is not marked tangent unless the surfaces touch along all of it.

Two surfaces sharing a boundary meet along it: the plain and the rational surface on one control net share every edge whose poles are collinear, and two neighbouring patches of a model share the edge between them. The march traced such a curve along a parameter bound of both surfaces at once. A predicted step left the parameter domain by a rounding, and the trace ended there, so the curve came back in pieces wherever rounding took a step out, and their number changed with the model's size and placement (ADR 0051): the plain and rational sample surfaces met in 10 curves at the reference size, 4 at a scale of 2^-20 and 42 at a translation of 1e3, with seeds no trace could leave reported as unseparated.

An edge is an iso-curve of its own surface, so its points are found on it directly, without a correction step that can leave the domain. Before seeding, each edge of either surface is checked at 9 points, its middle first; a point outside the other surface's box needs no projection, so an edge off the other surface costs one projection at most. An edge lying on it is sampled along its own parameter, cut at the knots of a B-spline surface and halved until the fitted Hermite curve keeps within a quarter of confusion of it, as a traced step does, and its parameters on the other surface come from projection. It is then kept as a traced curve: a seed on it is found on it and not traced, and a trace reaching it ends on the bound as before. An edge of the second surface that is also one of the first is kept once.

## Considered options

- **Cut the edge where other curves meet it**: rejected. The points where a crossing curve meets the edge are where the trace of that curve ends, found to within its last step; cutting the edge there makes its pieces depend on that step. A caller that needs the pieces splits the edge at the ends of the other curves.
- **Report the shared boundary as an overlap region**: rejected. Where the surfaces cross along the edge, it is a curve of the intersection, not a region where they lie on each other; where they touch along all of it, the curve is reported tangent.
- **Let a trace run along a parameter bound**: rejected. Holding a step on the bound needs a correction with one parameter pinned for each side the curve lies on the bound of, and the step would still decide where the curve ends. The edge is known whole before any step.
- **Edges lying on the other surface in part**: not now. An edge that leaves the other surface partway is traced as before. It is reopened when a pair shows it.

## Consequences

A pair whose edges lie on each other's surfaces reports each shared edge as one curve, with the same count at every size and placement. Every pair of surfaces with a bounded parameter pays the check of its edges: a box test, and one projection for an edge inside the other surface's box.
