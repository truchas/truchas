# Heat/Species Transport Refactor Context

This note is a handoff aid for continuing the thermal/species physics refactor.

## Package Shape

`thermal_species_driver.F90` is the non-OO entry point used by the backbone cycle
driver. It owns the retained thermal/species parameter list and a single
`thermal_species_physics` instance. `thermal_species_physics_type.F90` owns physics-level
coordination: mesh/MMF setup, solver dispatch, advection source transfer,
moving view-factor event forwarding, and state-copy accessors for output.

The concrete solver solvers are:

- `ht_solver`: fixed-domain heat transfer/phase change.
- `sd_solver`: fixed-domain scalar/species diffusion.
- `htsd_solver`: coupled heat plus species diffusion with Soret/coupling terms.
- `fht_solver`: heat transfer with fluid flow and dynamically changing void;
  this is deliberately separate because it has its own moving-domain behavior
  and private ODE strategy.

The vector data structure has been adopted across `ht`, `sd`, `htsd`, and
`fht`. `fht` differs because enthalpy is derived/intermediate, not a state
component coupled algebraically to temperature.

## Input And Initialization

The intended pattern matches other Truchas drivers: input reading constructs a
parameter list, and later driver initialization consumes it. Diffusion no longer
keeps separate `input_*` staging variables. `input_driver` reads
`physics_module` declarations (`heat_transport`, `species_transport`,
`number_of_species`) and passes them to `thermal_species_read_namelists`.

`thermal_species_read_namelists` builds only the active solver subtree:

- heat only: `solvers/heat`
- species only: `solvers/species`
- coupled heat/species: `solvers/htsd/model/{heat,species}`

`thermal_species_init` then passes the retained parameter list to `thermal_species_physics%init`.
Capability queries such as `thermal_species_have_heat_transport`,
`thermal_species_have_species_transport`, and `thermal_species_num_species` are initialized-package
queries only. Species initial-condition allocation now happens after
`thermal_species_init`, so it uses initialized physics capabilities.

## Important Conventions

Boundary-condition factories are concrete types:
`thermal_bc_factory` and `species_bc_factory`. The old abstract/factory1 split
has been removed.

The heat-transfer mesh name is exposed through `thermal_species_mesh_name()` because
microstructure must use the same mesh as heat transfer. It is currently `MAIN`,
but Truchas supports multiple meshes and this should remain a deliberate
coupling point.

`thermal_species_physics` should not retain parameter-list pointers unless needed.
The driver owns the retained parameter list; physics/solver init should use
local sublist pointers or let subordinate objects retain their own required
references.

## Validation Workflow

Default development toolchain is NAG:

```sh
ml purge && ml nag mpich truchas-tpl
cmake --build build-nag-debug --parallel 32
ctest --test-dir build-nag-debug -j32 -L 'DIFFUSION|CONDUCTION' --output-on-failure
```

NAG compilation needs escalated privileges for license-server access. MPI tests
also need escalated privileges. Before pushing, also test representative GCC,
Intel, and LLVM toolchains for portability.

The focused diffusion/conduction label set exercises heat-only, species-only,
coupled, FHT, radiation, microstructure, restart, and species-advection cases.

## Recent Cleanup State

Committed cleanup includes removal of `diffusion_solver_data.F90` and removal
of stored `params`/`solver_params` pointers from `thermal_species_physics`.

Uncommitted work at the time this note was written includes:

- removing `diff_set_input` and all `input_*` diffusion-driver staging state;
- passing physics declarations directly into `thermal_species_read_namelists`;
- making heat/species model parameter-list population solver-specific;
- moving species initial-condition allocation until after `thermal_species_init`;
- cleaning stale BC factory header text and `refactor_plan.md`.

That uncommitted state passed the NAG build and the focused
`DIFFUSION|CONDUCTION` tests.

## Deferred Design Work

The main deferred design question is deeper `htsd` decomposition/reuse. Since
`htsd` is heat plus species plus coupling, decide later whether it should reuse
more `ht`/`sd` internals or remain mostly separate as the coupled component.
Avoid starting that as incidental cleanup; it is a larger model/solver factoring
task.
