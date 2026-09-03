# Non-isothermal Navier--Stokes plan

This is a working design and implementation plan.  It records the current
decisions and will be updated as the work progresses.

## Scope of the first milestone

Implement thermal transport coupled to incompressible Navier--Stokes. The
thermal state is advected by the flow and then advanced by the implicit
thermal transport solver. The current implementation includes phase-aware
flow mapping, phase-specific enthalpy advection, and a single phase-change
material with one liquid phase and an immobile solid phase.

## Ownership and time-step structure

`material_distribution` remains simulation-wide, authoritative state.  Flow
uses a reduced view consisting of its individual mobile liquid materials, an
optional VOID component, and one derived SOLID component.  A flow material
transport component owns only the mapping, scratch/trial storage, and
material-resolved flux volumes; it does not own an independent distribution.

The coupled solver owns the attempted-step sequence.  Given the committed
time `t_n` and requested endpoint `t_np1`, it derives
`dt = t_np1 - t_n` locally.  Endpoint times are authoritative; no interface
advances a persistent time with `t = t + dt`.

The volume tracker may construct and evolve trial phase fractions internally,
but its sole authoritative external result is material-resolved face flux
volumes. The coupled solver advances `material_distribution` from the
divergence of those fluxes; it never adopts a tracker trial VOF. After thermal
transport, flow phase fractions are derived anew from the updated material
distribution and temperature.

For a coupled attempt:

1. Flow material transport reconstructs/advects its reduced phase view from
   the accepted face velocity and records phase flux volumes.
2. The coupled solver applies the flux divergences to the authoritative
   material distribution, then computes the enthalpy advection increment from
   the same fluxes and committed thermal state.
3. It attempts thermal transport. On thermal failure, it restores the
   pre-attempt material distribution and retries using the smaller step size
   supplied by the thermal solver.
4. On thermal success, it advances NS momentum and pressure correction using
   flow phase fractions derived from the advanced material distribution and
   thermal state. A flow-solver failure
   is non-recoverable for this algorithm.
5. It accepts the resulting material, thermal, and flow states.

Thermal transport retains its existing pending/commit capability for later
downstream physics that can reject an otherwise successful thermal step.

## Implementation stages

- [x] **Endpoint-time flow API.**  Refactor the existing Stokes and
  Navier--Stokes solvers to accept `t_n` and `t_np1` and derive their time
  step internally.  Update simulations and unit tests accordingly.
- [x] **Flow material transport.**  Add a reusable flow-side component that
  initially supplies one-liquid flux volumes directly from face velocity, then
  becomes the owner of the volume-tracker interaction and reduced-view mapping.
- [x] **External thermal rate.**  Add the general cell-integrated external
  enthalpy-rate interface to `ht_2d_model`/`ht_2d_solver`; zero remains the
  default for standalone heat transport.
- [x] **Enthalpy advection.**  Add a focused, unit-tested helper that converts
  material-resolved flux volumes into a conservative cell enthalpy increment,
  including explicit inflow thermal data when needed. Material distribution
  supplies the state from which flow transport constructs those fluxes; the
  advector itself needs only their material-slot mapping.
- [x] **Coupled solver and sim.**  `ns_ht_2d_solver` owns retry and
  sequencing, while `ns_ht_2d_sim`/main own setup, output scheduling,
  logging, signal termination, timing summary, and combined VTKHDF output.
- [x] **Initial tests.**  Cover zero-velocity equivalence to thermal-only
  results, material-resolved enthalpy-advection conservation, a driven-flow
  thermal case, serial/parallel reader consistency, and basic log output.
- [x] **Phase-aware state and output.** Derive flow phase fractions from the
  material distribution and temperature, use phase enthalpies for advection,
  and write phase volume fractions to VTKHDF output.
- [x] **Single-fluid-phase solidification.** Advance material distribution
  from phase fluxes, retain the parent material exactly in the one-liquid-phase
  case, and test coupled solidification in parallel.
- [ ] **Multifluid normalization.** Reconcile moving phase fractions after
  flux-divergence updates so they fill the non-solid remainder exactly, without
  deleting small fragments or adopting a clipped tracker trial VOF. Measure and
  report the correction.
- [ ] **Later capability.** Add temperature/distribution-dependent flow
  properties and buoyancy to phase-change problems, then extend tests to
  moving interfaces and multiple materials.
