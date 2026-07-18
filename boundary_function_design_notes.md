# Boundary Function Evaluation Design Notes

This note records the redesign motivated by species mass-transfer-coefficient
(MTC) boundary and interface functions, and the subsequent migration of a
selected family of thermal boundary functions.

## Original Issue

The species-only SD model has no temperature field, while the coupled HTSD
model has both concentration and temperature. Both models use the same
`species_component` type and therefore the same polymorphic MTC objects. The
old MTC interfaces required two mesh-wide state fields even for SD:

```fortran
call bc%compute(t, var1, var2)
```

A temporary implementation used a logical flag to make the lower-level MTC
object optionally ignore the first field. This avoided constructing a dummy
mesh-wide temperature array, but made a flag change the shape and meaning of
the user-function argument vector.

## Adopted Design

The sparse function object owns persistent topology, but evaluation results are
returned through allocatable arguments. The concrete implementation allocates
each result to the shape determined by its topology. The current abstract
classes are:

- `intfc_field_func`, with `index(:,:)`, for interface data evaluated with one
  mesh-wide field.
- `intfc_multifield_func`, with `index(:,:)`, for interface data that must work
  with either one or two mesh-wide fields.
- `bndry_field_func`, with `index(:)`, for boundary data evaluated with one
  mesh-wide field.

The interface MTC uses the multifield class so SD can supply concentration and
HTSD can supply concentration followed by temperature:

```fortran
call ic_mtc%compute_value(t, Cface, value)
call ic_mtc%compute_value(t, Cface, Tface, value)
```

The two-field call order is also the order in which the concrete implementation
slices the fields when constructing the anonymous argument array for the
user-specified scalar function. For the HTSD interface MTC that argument order
is `(C, T, t, x, y, z)`; for SD it is `(C, t, x, y, z)`.

The boundary MTC has analogous one-field and two-field interfaces. Its ambient
concentration depends only on `(t, x, y, z)`, while its MTC coefficient uses
the same concentration-first state ordering described above.

## Boundary Field Functions

The former `bndry_func2` class stored transient `value` and `deriv` arrays in
the object. It has been replaced for its thermal/species users by
`bndry_field_func`, whose interface returns callee-allocated results:

```fortran
call bc%compute_value(t, u, value)
call bc%compute_deriv(t, u, deriv)
```

The migrated concrete types are the external HTC, ambient radiation, oriented
flux, and evaporation heat-flux boundary functions. Existing physical behavior
was preserved except for the following corrections and derivative completion:

- Oriented-flux absorptivity is now evaluated with `u(index(j))`, the
  temperature on the face being processed. The old code incorrectly passed the
  entire mesh-local face array to a scalar function and effectively used an
  unrelated face temperature.
- The radiation derivative includes the temperature derivative of emissivity,
  evaluated by centered finite differences.
- The external HTC coefficient receives `(T, t, x, y, z)`, where `T` is the
  local boundary-face temperature. Its derivative includes the resulting
  product-rule term, with the coefficient derivative evaluated by centered
  finite differences. Ambient temperature remains a function of
  `(t, x, y, z)` only.
- The oriented-flux type can compute the temperature derivative of
  absorptivity, but it is intentionally omitted from the diffusion
  preconditioner. For incoming flux with absorptivity increasing with
  temperature, this tangent is negative and can destroy coercivity. The full
  nonlinear residual retains the physical thermal feedback.

Boundary and interface evaluations are entity-local. A value for one managed
face, or one managed pair of matching interface faces, depends only on the
supplied fields at those entities. There is no coupling to other managed
entities.

## Why Return Results

Returning allocatable results makes transient call state explicit, avoids
exposing mutable arrays through a polymorphic object, and keeps topology-based
result sizing out of the caller. It also lets one shared abstract type provide
the field signatures required by both SD and HTSD without constructing
mesh-wide packed state arrays.

The tradeoff is that caching is no longer implicit. Some other sparse boundary
functions cache constant or time-dependent results in the object, and complex
implementations may still benefit from private caching. Such a cache should be
an implementation detail with clearly defined lifetime and invalidation, not a
public result component.

## Remaining Questions

- Audit the other sparse boundary-function families before migrating them;
  caching may be materially valuable for some implementations, especially in
  electromagnetics.
- Audit which sparse functions permit duplicate indices. Residual assembly must
  use explicit loops rather than indexed assignment. The current diffusion
  matrix diagonal increment also uses an explicit loop and accumulates
  duplicate entries.
