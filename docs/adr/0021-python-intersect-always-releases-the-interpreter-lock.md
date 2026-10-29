# Python `intersect` always releases the interpreter lock

The Python `intersect` releases the interpreter lock around every core call, whatever the pair: 3D curve pairs, curve and surface pairs, surface pairs and 2D curve pairs. While one thread runs a numeric search, other Python threads keep running.

Before this decision the bindings released the lock only for the pairs they expected to search, through their own copy of the core routing: a B-spline or Bézier operand, or a curve other than a whole line against a torus. That copy had already drifted from the core. The 2D pairs were never in it, so a 2D B-spline intersection held the lock for the whole search, and every new variant (the Bézier shapes, ADR 0017) had to be added in two places. Whether a pair has a closed form is the router's decision (ADR 0011). The bindings should not have to know it.

Releasing the lock costs about 47 ns per call: in a tight Python loop, a closed-form `Line` × `Plane` call goes from about 290 ns to about 340 ns (release builds, minimum of eight interleaved rounds on a loaded machine). That is a 12 to 16% slowdown on the cheapest calls. A searching call takes from tens of microseconds to milliseconds, so there the cost does not show.

## Considered options

- **Keep the copied routing and add the 2D pairs**: rejected. It fixes today's bug and keeps the closed-form calls 47 ns faster. But the next routing change, a new exact pair or a new variant, has to be mirrored by hand, and a miss is silent: nothing fails, a thread just stops running in parallel.
- **Hold the lock when neither operand is a B-spline or Bézier shape**: rejected. It reads only the operands, not the routing, but elementary pairs search too: a circle against a torus searches for about 130 µs, and a sphere against a torus for about 1 ms, all while holding the lock.
- **Run the closed form under the lock, and release it only for the search**: deferred. It needs a core entry point that tries the exact solver and returns a family result, or says that the pair has none. That conversion from exact to family results is being reshaped in the core, so the bindings would be built on code that is about to change. Once the core has it, this option gets back the 47 ns without copying the routing.

## Consequences

Any thread can run `intersect` on any pair without blocking the others. A pair the core routes differently in the future needs no change in the bindings. Closed-form calls from Python cost about 47 ns more until the core offers a closed-form-first entry point.
