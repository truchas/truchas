# Multimaterial Volume-Fraction Design Summary

The simulation owns the authoritative material-composition state and supplies
it to all coupled physics models.

## Authoritative composition

- The state stores material volume fractions, not material volumes or phase
  fractions.
- Materials are the conserved identities. Phase composition is derived by the
  material model from temperature and other state.
- The initial implementation uses a dense `real(r8)` array over simulation
  materials and owned cells:

  ```fortran
  volume_fraction(nmaterial, ncell_onP)
  ```

- This is a reference implementation behind an interface that permits later
  replacement by partition-local or block-sparse storage.
- Ghost material fractions are not generally stored or exchanged. Consumers
  derive the quantities they need on owned cells and communicate those
  quantities.

## Physics-specific views

Thermal and similar models evaluate effective properties on owned cells using
the material fractions and material-property models. They then gather the
resulting scalar property values to ghost cells.

Flow uses a separate compact dense view containing:

- each immiscible liquid phase;
- `VOID`, if present; and
- aggregate `SOLID`, if present.

This flow view has a fixed row ordering across ranks and is explicitly gathered
to ghosts. All liquid phases move with the same velocity but are reconstructed
and advected separately; each liquid phase is conserved during that advection
step. Liquid-phase fluxes are reduced through their fixed parent-material
mapping to update the authoritative material fractions.

Flow does not distinguish individual immobile materials; it sees only the
aggregate solid fraction. Phase transformation is handled separately from
advection by the material and thermodynamic models.

## Initialization

Initialization is separate from runtime composition storage:

1. Ordered cell-set and geometric regions are associated with material IDs.
2. Cut-cell calculations produce region fractions.
3. Region fractions associated with the same material are summed.
4. The result initializes the simulation's material volume fractions.

The existing `material_distribution` directory provides a useful prototype for cell
sets, boxes, half-planes, disks, background fill, and cut-cell subdivision, but
it must be audited and tested before adoption.

## Invariants and conservation

Composition invariants are established whenever the state is changed:

```text
0 <= alpha(m,c) <= 1
sum_m alpha(m,c) = 1
```

For now, fractions remain real-valued. Updates apply the cutoff and boundedness
policy and assign any small closure residual deterministically, preferably to
the largest component. Material conservation is measured using

```text
sum_c volume(c) * alpha(m,c).
```

Conceptually this state serves the role of a material mesh function, but it is
a new, narrower abstraction: composition storage and access only. Region
construction and property-evaluation machinery are separate concerns.
