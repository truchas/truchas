# 2D Truchas TODOs

This is the cross-cutting capability list for the 2D flow and coupled
flow/thermal models.  More detailed subsystem items remain in
`volume_tracking/TODO.md` and `nonisothermal/TODO.md`.

## Multimaterial and multiphysics capability

- Extend `ns_ht_2d` beyond the current single-phase material contract.  The
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

- Implement nested-dissection reconstruction for cells containing three or
  more materials; the current geometric tracker uses onion-skin ordering.
- Revisit conservation repair and small-volume threshold behavior as the
  number of mobile and immobile components grows.
- Complete the flow-level integration of VOID and immobile SOLID fractions,
  including their effects on transport and property evaluation.

## Spatial discretization and initialization

- Improve the flow operators for triangle and mixed triangle/quad meshes.
  Quad meshes are currently the supported and tested flow topology.
- Improve the no-slip discretization at solid/fluid interfaces.  The current
  treatment places the effective wall at the solid-cell center rather than at
  the material interface; [the dormant solid-wall Poiseuille test](flow/Test/ns/poiseuille_solid_wall.json)
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
