# Marching is a surface pair, a tracer and a seeding policy

The marching intersection (`intersect/numeric/marching.rs`) composes three crate-private modules under `intersect/numeric/marching/`, and runs the coincidence check (`coincidence.rs`) before them:

- `pair.rs`, the **surface pair**: the two surfaces and the point `x = [u1, v1, u2, v2]` on both. It converges a guess onto both surfaces (`meet`, with or without a plane through the predicted point; `meet_touching` where they touch; `land` on a parameter bound), clamps parameters into their bounds, predicts a step, and reads the local geometry: the direction where the surfaces cross (`tangent`) or touch (`contact_direction`), and whether they touch (`grazes`) or lie on each other (`osculates`). It knows nothing of steps or seeds. `swapped` and `swap` read the pair the other way round, once.
- `tracer.rs`, the **tracer**: the step loop. From a seed it traces both ways and returns a `Trace`: a curve with whether the surfaces touch along it, a touch at a point, or nothing. How the surfaces meet where a step is taken (crossing or touching) is the tracer's own per-step state, a private `Contact` value read at every step; it picks the direction and the corrector, and the pair and the caller never see it.
- `seeding.rs`, the **seeding policy**: which sources run on a pair (patch pairs on two B-splines, the iso-curves of an elementary side, the parameter grids), in what order, and when one is left out (the iso-curves where the grid seeds only graze; the grids once the patch-pair search exhausts its budget). It returns the seeds in the order they are traced. A patch-pair seed stays unconverged, with the box that skips it once an earlier curve crosses it, so that it is converged only when its turn comes, as before.

`marching_intersection` asks the policy for the seeds, hands each to the tracer, and fits and collects what it returns. Every result is bit for bit what it was before the split.

The split places the two roadmap items of marching. Tracing through a point where the surfaces touch changes only the tracer: a step that fails in one mode re-reads the pair (`grazes`) and switches its `Contact`, and the tracer reports the curve. Finding loops smaller than the iso-curve spacing changes only the seeding policy: a new source returns more seeds.

The surface pair stays generic over the two surface types, as the functions it replaces were, so the Newton loops are monomorphized for each pair of types and their evaluations inlined. Instruction counts on every marching bench moved by less than one percent.

## Considered options

- **One file of free functions** (before this decision): rejected. Its 28 generic functions each repeated the two surfaces and their bounds; the touching mode, decided once per seed, was threaded through four signatures; the seeding policy was spread over three functions; and the step loop could only be tested by running a whole march.
- **The touching mode on the pair, or passed by the caller**: rejected. Either way a mode change in the middle of a trace would cross a module seam. On the tracer it is one variable of the step loop.
- **The pair over `&dyn ParametricSurface`**: rejected. Every evaluation in the correctors would become an indirect call that cannot be inlined, in the innermost loops of the search.
- **Seeds as plain converged points**: rejected. A patch-pair seed is skipped when an earlier trace already crosses its patches, and is converged only when it is reached; converging every patch-pair seed up front would cost Newton runs the march skips today.
