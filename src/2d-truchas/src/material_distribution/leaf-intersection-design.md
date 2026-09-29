# 2D volume-fraction leaf intersection

## Purpose

`t2d_vol_frac_init_procs` currently computes fractions by recursively subdividing
mixed triangles.  At the refinement limit it assigns a mixed triangle to the
region containing its centroid.  This note describes a higher-accuracy leaf
treatment modelled on the final-subtet treatment in mainline's
`compute_body_volumes_proc`.

The design retains recursive subdivision as the universal initialization
algorithm.  In particular, a region remains usable when it can answer only
the existing point-membership question.  No geometry-specific cases belong in
the subdivision algorithm.

## Existing ordered-region semantics

The input region list is an ordered decomposition.  At a point, the selected
region is the first specification that contains the point; later
specifications apply only to the portion not assigned by preceding ones.
This is the same precedence rule used by mainline VOF body initialization.

The computed volume fractions consequently describe this ordered assignment,
not the independent geometric overlap of the input regions.

## Optional implicit-boundary capability

Extend the abstract `region` interface with an optional capability to
evaluate a scalar implicit function, `f(x)`, for its geometric boundary:

* `f(x) < 0` denotes the enclosed side;
* `f(x) = 0` denotes the boundary; and
* `f(x) > 0` denotes the exterior side.

The sign convention must agree with `encloses`.  `f` need not be a Euclidean
signed distance: a consistently signed level-set function is sufficient.  A
true signed distance is convenient, but is not required by the intersection
calculation.

A region that cannot provide this information continues to use the current
centroid rule at a mixed leaf.  The precise Fortran API should be decided
when implemented.  Plausible forms are a deferred `implicit_value` binding
plus a `has_implicit_value` query, or an abstract binding returning an
availability flag.  It must be possible for `cell-set` regions and the
background region to decline the capability.

`region_func` will need an internal operation to obtain this information for
the selected region index.  The volume-fraction code should not select on
concrete region type.

## Mixed-leaf algorithm

At the recursion limit, retain the present vertex classification by the full
ordered region list.  If the triangle contains exactly two assigned region
indices and the earlier (higher-precedence) of those regions supplies an
implicit function:

1. Evaluate `f` at the three triangle vertices.
2. Locate zero crossings on the edges by linear interpolation of `f`.
3. The two crossings define a line segment, hence a line separating the two
   locally assigned regions.
4. Clip the triangle against that line and compute the two areas.
5. Assign the clipped area to the higher-precedence region and the remainder
   to the other region.

The clipping calculation is exact for the *linearly interpolated* implicit
function.  It is geometrically exact for a planar interface and is a local
linear approximation for a curved boundary.  It is the 2D counterpart of
mainline's construction of a plane from interpolated signed-distance edge
intersections followed by tet/plane intersection.

All other cases retain the centroid fallback: unavailable implicit function,
more than two sampled regions, or a degenerate/ambiguous edge-intersection
configuration.  This provides a conservative implementation path with no
change to the treatment of general user-defined regions.

## Interface detection remains separate

The procedure currently refines only when its sampled vertices have different
assigned regions.  Thus it can miss a feature wholly enclosed within a
triangle, just as mainline's vertex-based divide-and-conquer detection can
miss an interface that does not alter vertex body IDs.  Leaf intersection
improves the fraction after a mixed triangle has been detected; it does not
solve this detection limitation.

The existing user guidance still applies: geometric features should be
resolved relative to the local mesh scale.  If stronger detection guarantees
are wanted later, they require additional region information (for example,
intersection/bounding queries) or additional sampling, and should be designed
independently of the leaf treatment.

## Possible planar short-circuit

For a half-plane, `f` is affine, so the clipping result is exact at *any*
subdivision level.  It is tempting to stop recursion as soon as a mixed
triangle is associated with such a boundary.

That optimization is not automatically safe under ordered-region semantics.
Although two vertex labels may be present, a later region can occupy part of
the remainder without being sampled at those vertices; ordinary subdivision
could discover it.  A planar short-circuit therefore needs a proof that no
other region can affect the triangle's remainder, or an explicit restricted
case.  Examples of safe initial restrictions include a two-region
decomposition consisting of the planar region and background, or a suffix of
the ordered list known to be inactive on the triangle.

Accordingly, the first implementation should use exact clipping only at the
existing refinement leaf.  A later optimization may add a certified
short-circuit through a `region_func`-level query; it should not be introduced
as a concrete `t2d_half_plane_region` test in `t2d_vol_frac_init_procs`.

## Validation

Unit tests should cover a triangle cut by an affine implicit function at
several orientations, including cuts through a vertex and along an edge.
Compare the two fractions with analytic clipped areas and verify nonnegative
fractions summing to one.  Integration tests should retain the current
recursive path for regions without an implicit function and exercise ordered
precedence with at least two geometric specifications.
