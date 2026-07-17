# Heat/Species Transport Refactor Plan

This note records the intended component split before extracting a scalar
diffusion solver from the current `HTSD_solver`. It is a working guide for the
refactor, not user-facing documentation.

## Current Dispatch

`thermal_species_physics` owns physics-level state, material mesh-function exchange,
advection source transfer, event forwarding, and the non-OO entry points used by
the cycle driver. It selects one concrete solver:

- `ht_solver`: fixed-domain heat transfer and phase change with no species.
- `HTSD_solver`: fixed-domain species diffusion, with or without heat transfer.
- `fht_solver`: heat transfer with fluid flow and dynamically changing void.

## Target Split

Create `sd_solver` first and leave `HTSD_solver` unchanged initially. The new
solver should own scalar-only model construction, integration, preconditioning,
state access, restart, and source/advection updates for species diffusion. Once
that path is stable, `HTSD_solver` can be narrowed to the coupled heat/species
model and renamed or cleaned up later if useful.

The physics-level dispatch should then become:

- heat only: `ht_solver`
- scalar only: `sd_solver`
- coupled heat and scalar: `HTSD_solver`
- heat with moving void/fluid flow: `fht_solver`

Vector adoption is intentionally deferred until this component split is stable.

## Shared Code Boundaries

The split should preserve shared implementations where the behavior is common:

- MFD discretization and matrix/preconditioner infrastructure:
  `mfd_disc_type`, `mfd_diff_matrix_type`, `mfd_diff_precon_type`.
- Boundary-condition factories and policies:
  `thermal_bc_factory_type`, `species_bc_factory_type`.
- Source factories:
  `thermal_source_factory_type`, `species_source_factory_type`.
- Material/property views:
  `matl_mesh_func`, `prop_mesh_func_type`, and mesh interop helpers.
- Integrator-facing structure:
  IDAESOL model wrappers, norms, preconditioners, and restart/state-copy
  patterns should be separated only where scalar-only and coupled behavior
  genuinely diverge.

Avoid introducing the new vector representation while extracting `sd_solver`;
that should be a later, behavior-preserving pass.

## Regression Anchors

Existing tests cover the current solver modes and should be run around the
split:

- Scalar-only diffusion: `ds1`, `ds5`, `diffusion-mtc`, `species-src-1`.
- Scalar advection/diffusion: `species-adv-1`, `species-adv-2`,
  `species-adv-3`.
- Coupled heat/species: `ds3`, `ds6`.
- Heat-only fixed-domain path: `ds2`, `ds7`, `ds9`, `ds10`.
- Heat phase change: `ds4`, `ds11`, `enthalpy1`, `enthalpy3`, `enthalpy4`,
  `enthalpy5`.

Example focused command:

```sh
ctest --test-dir build-nag-debug -j32 -R \
  '^(ds[1-7]|ds9|ds10|ds11|diffusion-mtc|species-adv-[1-3]|species-src-1|enthalpy[1345])$' \
  --output-on-failure
```

`enthalpy2` is intentionally omitted from the NAG-focused command because it has
a known long-standing NAG numerical delta unrelated to this refactor.
