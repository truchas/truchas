# Non-isothermal Navier--Stokes plan

This is a working design and implementation plan.  It records the current
decisions and will be updated as the work progresses.

## Scope of the first milestone

Implement thermal transport coupled to incompressible Navier--Stokes for one
liquid material.  The thermal state is advected by the flow and then advanced
by the existing implicit thermal transport solver.  The initial milestone has
constant flow properties: thermal feedback to momentum, phase change, and VOF
interface reconstruction follow later.

## Ownership and time-step structure

`material_composition` remains simulation-wide, authoritative state.  Flow
uses a reduced view consisting of its individual mobile liquid materials, an
optional VOID component, and one derived SOLID component.  A flow material
transport component owns only the mapping, scratch/trial storage, and
material-resolved flux volumes; it does not own an independent composition.

The coupled solver owns the attempted-step sequence.  Given the committed
time `t_n` and requested endpoint `t_np1`, it derives
`dt = t_np1 - t_n` locally.  Endpoint times are authoritative; no interface
advances a persistent time with `t = t + dt`.

For a coupled attempt:

1. Flow material transport advects the mobile material fractions from the
   accepted face velocity and records material flux volumes.
2. The coupled solver computes the enthalpy advection increment from those
   fluxes and the committed thermal state, then attempts thermal transport.
3. On thermal failure, it restores the pre-advection mobile composition and
   retries the entire attempt using the smaller step size supplied by the
   thermal solver.
4. On thermal success, it advances NS momentum and pressure correction using
   the advanced material composition and thermal state.  A flow-solver failure
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
  including explicit inflow thermal data when needed. Material composition
  supplies the state from which flow transport constructs those fluxes; the
  advector itself needs only their material-slot mapping.
- [x] **Coupled solver and sim.**  `ns_ht_2d_solver` owns retry and
  sequencing, while `ns_ht_2d_sim`/main own setup, output scheduling,
  logging, signal termination, timing summary, and combined VTKHDF output.
- [x] **Initial tests.**  Cover zero-velocity equivalence to thermal-only
  results, material-resolved enthalpy-advection conservation, a driven-flow
  thermal case, serial/parallel reader consistency, and basic log output.
- [ ] **Later capability.**  Replace the one-liquid transport bridge with VOF,
  add temperature/composition-dependent flow properties and buoyancy, then
  extend tests to multi-material and moving-interface cases.
