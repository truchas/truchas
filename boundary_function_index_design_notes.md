# Boundary Function Index Design Notes

## Recommended Public Contract

The public index of a sparse boundary or interface function should contain
each managed mesh entity exactly once. Each returned result should correspond
one-to-one with that index.

Duplicate indices should not escape a concrete boundary-function object. This
restores the ordinary mathematical meaning of a sparse function: one value for
each entity in its domain.

For an interface function, each column of the public index represents one
managed matching entity pair and should likewise be unique.

## Additive Superposition

Some concrete types, such as oriented thermal flux, intentionally support
multiple contributions on the same entity. That overlap should remain
supported, but it should be an internal implementation detail.

Such a concrete implementation should:

1. Preserve a private grouped contribution list, which may contain repeated
   entity indices.
2. Construct a unique public `index`.
3. Construct a private `map` from each grouped contribution to its position in
   the public index.
4. Initialize each computed result array to zero and accumulate contributions
   through that mapping.
5. Accumulate derivative contributions in the same way when applicable.

For example, an implementation might privately maintain data equivalent to:

```fortran
integer, allocatable :: index(:)          ! public, unique entities
integer, allocatable :: map(:)            ! private result slot for each entry
integer, allocatable :: xgroup(:)         ! private specification segments
```

Multiple lasers still superimpose, but the caller receives one total
oriented-flux value per face. EM `nxH` uses the same pattern for additive edge
contributions. PEC likewise exposes the unique union of its zero-valued edges,
without requiring an accumulation map after construction.

This is preferable to exposing an `indices_unique` property. Such a property
would force every caller to branch on a representation detail and would leave
duplicate handling vulnerable to mistakes in new consumers.

## Input Overlap Policy

Input-policy decisions should be separate from the storage representation.
The factory should decide whether overlap is legal, while the concrete object
should reduce any legal superposition to one public value per entity.

Recommended defaults are:

- Ordinary HTC, radiation, scalar flux, and interface specifications reject
  overlap between specifications of the same category.
- Oriented flux explicitly permits same-category overlap because its intended
  model includes additive contributions from multiple sources.
- Different natural-flux categories remain independently additive.
- Dirichlet-like data reject ambiguous overlap.

Rejecting overlapping HTC specifications is useful input validation, not
excessive protection. Although two HTC terms are mathematically additive,
overlap is more likely to indicate intersecting face sets or duplicated input.
Silently accepting it could produce a plausible but unintended result.

If a genuine need for multiple same-category heat-transfer paths arises, that
behavior should be enabled explicitly and documented as superposition. Until
then, rejection is the safer default.

The builder interface should express the representation directly. The direct
form may return repeated entity indices when `no_overlap=.false.`, while a
separate unique form returns a unique entity array and a grouped map into it.

## Benefits

The unique public-index invariant provides several concrete benefits:

- Mask operations such as `mask(index) = .true.` are well-defined.
- Callers do not need to know which concrete implementations support
  superposition.
- Assignment versus accumulation semantics are no longer implicit.
- New consumers cannot accidentally mishandle duplicate entries.
- Cross-category accumulation remains visible at the physics-component level.
- An array containing one boundary-function object per input specification is
  unnecessary.

The base classes should document this invariant. Concrete finalization should
verify it with debug assertions. Face-to-node and face-to-edge implementations
should also deduplicate or reject overlap before publishing their indices.

## Tests

Focused tests should cover:

- Rejection of overlapping HTC specifications.
- Superposition of two overlapping oriented-flux specifications.
- Uniqueness of the resulting public oriented-flux index.
- Equivalence between the internally accumulated result and explicit summation
  of the individual contributions.
- The analogous additive edge case for EM `nxH`.
- Debug invariant checks for every migrated concrete boundary and interface
  function.
