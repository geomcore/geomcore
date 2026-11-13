# A result pickles with its operands, in its own version

Every Python result pickles and copies: the intersections and their hits, overlaps, regions, completeness and curves, the projections, the curve–curve extrema, the surface curvatures, the approximations and the discretizations. A result that holds its operands (`first` and `second`, or `target`) pickles them with it, as the shapes they are. So a `multiprocessing` worker can return an intersection, and a bug report can attach one.

A result reduces to its class's private static `_from_state`, called with the geomcore version that pickled it, its operands, and its parts as plain tuples, lists, numbers and names. The curves and surfaces inside a result, such as an intersection curve, its curves in parameters, or an approximation's B-spline, are shapes and pickle as themselves. In Rust, the six results whose parts are private, and the two `non_exhaustive` discretizations, gain a `from_parts` constructor. It takes the parts their getters return and keeps them as they are, unchecked and unsorted, as a deserializer does.

A result's pickle loads only in the geomcore version that wrote it. Its state is positional and follows the Rust result's parts, which change from one release to the next (a hit gained its multiplicity, a completeness its reason). `_from_state` refuses a pickle from another version with a `ValueError` naming both versions, rather than misreading it. Shapes and values are unchanged: they reduce to their public constructors, and load in a later version as long as those constructors take the same arguments.

## Considered options

- **Recompute the result from its operands on load**: rejected. A worker would compute an intersection twice, once to return it and once to load it, and a cache would save nothing.
- **Serialize the Rust results with serde**: rejected. It adds a dependency to the library and derives on every curve and surface, for one consumer that pickles in Python, and a serde layout would not load across versions either.
- **Keep results unpicklable and document the arrays to take them apart into**: rejected. That is what every caller did by hand, and it loses the completeness, the regions and the curves in parameters.
- **Load pickles across versions**: not now. It would need a versioned state for every result and a migration for each change of its parts. It is reopened if a consumer stores results across releases.

## Consequences

`from_parts` lets a Rust caller build a result that no operation would return, such as an exact result that is incomplete. It is documented as the deserializer's constructor. Each result class has a hidden `_from_state`, which the stubs leave out like other private names. A pickle written by one release fails to load in the next with a message that names both versions.
