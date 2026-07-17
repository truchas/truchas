# Truchas Source Architecture

`src/truchas` builds mostly as one monolithic Fortran library, `truchas`, with
subsystem directories contributing sources through `add_subdirectory`. The
executable is a thin wrapper around `drivers/main.F90` and links that library.

At runtime, control is centralized in `drivers`. `main` calls `drivers:code`,
which initializes MPI communication and HYPRE, processes command-line input,
initializes logging and timers, reads the input deck, runs setup, enters the
time-cycle driver, then cleans up. The main time loop is a multiphysics
coordinator rather than a deeply abstract scheduler: `CYCLE_DRIVER` directly
sequences EM heat updates, volume tracking, flow, diffusion/heat/species, solid
mechanics, probes, output, event queues, and time-step synchronization.

Input is centrally orchestrated by `input_driver`. It reads namelist groups in
dependency order: functions, physics flags, materials, outputs, restart, meshes,
interfaces and bodies, numerics, then optional physics-specific namelists for
flow, solid mechanics, EM heating, diffusion/heat/species, microstructure, and
probes.

Setup creates the simulation state after input parsing. It initializes the mesh
manager, allocates base cell and material state through `zone_module` and legacy
material APIs, writes the primary mesh, applies restart time state, and calls
`INITIAL` for initial conditions.

## Major Subsystems

- `distributed_mesh`, `exodus`, `partitioning`, `communication`: mesh
  representation, mesh IO, METIS partitioning, and MPI communication helpers.
- `materials`: material database/model, material factories, phase change, and
  averaged material property evaluation.
- `thermal_species`: thermal/species transport models,
  thermal/species boundary conditions and sources, MFD discretization, and
  BDF2/IDA-style models.
- `enclosure_radiation`: view-factor enclosure radiation models and solvers.
- `physics/volume_tracking`: volume-of-fluid/interface tracking and geometric
  material advection.
- `flow`: incompressible flow state, boundary conditions, turbulence models,
  projection, and prediction algorithms.
- `solid_mechanics`: node-centered finite-volume thermoelastic/viscoplastic
  mechanics with contact boundary-condition handling.
- `electromagnetics`: induction/microwave heating, time/frequency-domain Maxwell
  solvers, EM property evaluation, and Joule/dielectric heat coupling.
- `solver` and `ode`: shared numerical infrastructure including HYPRE wrappers,
  sparse matrices, nonlinear solvers, vectors, MINRES/CG, NKA, and BDF2/DAE
  integration.
- `functions`, `mesh_functions`, `toolpath`: user-defined scalar/vector
  functions, mesh-attached boundary/interface functions, and AM/CNC-style
  toolpath support.
- `output`, `truchasio`: cycle/edit/probe output plus HDF5, Scorpio, and fVTKHDF
  facing IO.

## Design Pattern

The codebase is transitional. Newer subsystems use derived types, factories, and
parameter-list style construction, while the top-level driver still uses
singleton module state and broad `use` association. This is visible in the
materials driver, global physics flags, and central cycle driver.
