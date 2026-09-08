# 2D Truchas TODOs

This is the cross-cutting capability list for the 2D flow and coupled
flow/thermal models.  More detailed subsystem items remain in
`physics/volume_tracking/TODO.md` and `simulations/flow_thermal/TODO.md`.

## Multimaterial and multiphysics capability

- Extend `flow_thermal` beyond the current single-phase material contract.  The
  material factory already has a parameter-list representation for
  multiphase materials, but the 2D flow path must first split each parent
  material fraction into temperature-dependent phase fractions.  The flow
  layout already distinguishes phase IDs from parent material IDs and treats
  SOLID as the residual after fluid and VOID fractions; the reverse mapping
  from that lumped slot remains to be designed.
- Add coupled passive-scalar/species transport.  Flow properties and buoyancy
  will eventually depend on temperature and concentrations; concentrations
  are not material identities.
- Decide whether a full variable-density formulation is needed beyond the
  current incompressible Boussinesq model.  The current constant density and
  constant thermal-expansion parameters are appropriate for Boussinesq flow.
- Add phase-change and other thermal/material coupling needed for parity with
  the 3D multiphysics capability.
- Add surface-tension and capillary forces for interface-driven flow.

## Volume tracking and interfaces

- Implement the VOID-collapse model required for the full free-surface flow
  capability.  The current 2D path supports VOID tracking and inactive-VOID
  flow cells, but does not yet implement the collapse treatment used by the
  mature 3D flow algorithm.
- Implement nested-dissection reconstruction for cells containing three or
  more materials; the current geometric tracker uses onion-skin ordering.
- Revisit conservation repair and small-volume threshold behavior as the
  number of mobile and immobile components grows.
- Complete the flow-level integration of VOID and immobile SOLID fractions,
  including their effects on transport and property evaluation.

## Spatial discretization and initialization

- Improve the flow operators for triangle and mixed triangle/quad meshes.
  Quad meshes are currently the supported and tested flow topology.
- Replace the current free-surface pressure-force split with a single
  face-consistent, rotationally accurate formulation.  The projection matrix
  and face-velocity correction use face densities and face pressure
  derivatives, but the current cell-centered pressure-force path uses the
  first-order `gradient_cf` reconstruction, followed by face-density scaling
  and interpolation.  This bypasses the improved rotationally invariant
  `gradient_cc` reconstruction whenever VOID is present.  See
  `physics/flow/pressure-correction-comparison.md` for the current analysis.
- Eliminate the global per-step switch between the ordinary cell-centered
  gradient path and the free-surface face-based path.  The switch can change
  when a small VOID fraction appears or disappears and therefore can cause a
  discrete change in the cell velocity correction and time-step trajectory.
  A unified operator should handle pure fluid, mixed fluid/VOID, and
  fluid/VOID interfaces continuously.
- Establish a geometric and resolution-aware basis for the minimum face
  fraction regularization.  Its relationship to the material cutoff, mesh
  resolution, and the admissible velocity error should be documented and
  tested rather than treating the current value as a purely empirical
  constant.
- Add focused regressions for mixed fluid/VOID domains that contain no pure
  VOID cells, including interface appearance/disappearance and perturbed
  meshes.  These should check pressure, velocity, projection convergence, and
  sensitivity to the face-fraction regularization.
- Improve the no-slip discretization at solid/fluid interfaces.  The current
  treatment places the effective wall at the solid-cell center rather than at
  the material interface; [the dormant solid-wall Poiseuille test](simulations/flow/test/ns/poiseuille_solid_wall.json)
  demonstrates the resulting velocity error.
- Continue improving accuracy on non-orthogonal meshes and reduce the finite-
  resolution hydrostatic well-balance residual.
- The cell-centered flow gradient currently uses a minimum-norm `DGELSY`
  fallback for rank-deficient stencils created by dynamically changing solid
  regions.  Replace this with a proper interface-aware reconstruction; the
  fallback is not a general treatment of solid/fluid interfaces.
- Revisit the initial pressure/velocity construction so that hydrostatic and
  boundary-compatible initial states are treated as accurately as possible.
- Handle pressure null spaces for multiple disconnected active-fluid
  components (and, more generally, multiply connected active domains).  The
  current dynamic pin selection supplies one reference per projection solve.

## Architecture and remaining legacy paths

- Retain one shared flow mechanics implementation for Stokes, isothermal NS,
  and non-isothermal NS while exposing the separate material-transport,
  momentum, and projection phases needed by coupled physics.
- Generalize the coupled solver structure as additional physics models are
  added, rather than introducing separate ad hoc combinations for each pair
  of models.
