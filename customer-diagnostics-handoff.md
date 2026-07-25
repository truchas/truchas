# Customer Diagnostics Handoff

Context for a later session: this branch contains temporary diagnostics for a customer-reported apparent MPI hang during unstructured mesh instantiation.

## Branch and Commit

- Branch: `24.04-wnaa`
- Base tag: `24.04`
- Diagnostic commit: `2498cc00` (`Add mesh instantiation diagnostics`)
- Remote branch: `origin/24.04-wnaa`

## Customer Problem Summary

The customer reports an apparent MPI hang while instantiating an unstructured mesh in:

```text
src/truchas/distributed_mesh/unstr_mesh_factory.F90
```

The mesh has about 4M cells. Smaller meshes have run successfully, and the customer also reports that a larger 8M-cell mesh has run successfully. They specifically described this failing mesh as more geometrically complex.

The last visible terminal output was:

```text
partitioning the mesh nodes
```

That output may have been misleading because TLS output was not guaranteed to flush before these changes. Still, it suggests execution likely got past METIS cell partitioning. The physics is simple heat transfer with viewfactor radiation.

A preliminary viewfactor computation using the same mesh completed successfully after extracting a boundary surface mesh. That suggests basic mesh input and surface extraction can succeed, but it does not rule out pathological volume-mesh node support, face/ghost construction, or rank-0 memory pressure during distributed mesh instantiation.

## Diagnostic Changes Added

The diagnostic commit changes these files:

```text
src/truchas/distributed_mesh/unstr_mesh_factory.F90
src/truchas/output/cycle_output_module.F90
src/truchas/utilities/get_process_size.c
src/truchas/utilities/process_info_module.F90
src/truchas/utilities/truchas_logging_services.F90
suggested-run-config.md
```

Main diagnostics:

1. Added flushed phase markers through `new_unstr_mesh_aux`.
   These report progress after major post-METIS phases such as node partitioning, face numbering, face partitioning, ghost-cell identification, index-map initialization, permutation distribution, cell-node/cell-face initialization, and link initialization.

2. Added rank-0 memory stamps to the mesh phase markers.
   They print virtual, resident, and data sizes as:

   ```text
   memory: vsize=... MB, rsize=... MB, dsize=... MB
   ```

3. Added node-support diagnostics in `all_cell_neighbors`.
   These report the top nodes by number of incident cells and exact counts for candidate cell-pair expansion, cross-partition pairs, nodes crossing partitions, and max cross-partition pairs for any node.

4. Made TLS output flush after ordinary output paths.
   This should make the last printed line much more trustworthy if the run hangs, aborts, or is killed.

5. Fixed and verified the old memory helper.
   `get_process_size` now uses 64-bit values and was checked against `/proc/self/statm`.

## Build Used For Verification

The requested toolchain was:

```bash
ml purge
ml gcc/15.2 mpich truchas-tpl/v22-wnaa
```

Fresh build directory used:

```text
/tmp/truchas-24diag-gcc152-v22-wnaa
```

Build command:

```bash
cmake --build /tmp/truchas-24diag-gcc152-v22-wnaa --parallel 32
```

The build completed successfully.

## Suggested Customer Run Configuration

The note in `suggested-run-config.md` recommends using MPICH abort/error checking and rank-labeled output:

```bash
export MPICH_ABORT_ON_ERROR=1
export MPICH_ERROR_CHECKING=1
mpiexec -n 64 -l ./truchas input.inp > truchas.out 2> truchas.err
```

The actual MPI size and executable/input paths should match the customer's run.

## What To Look For In Returned Output

First, identify the last flushed mesh phase marker. If output stops after a specific marker, inspect the next phase in `new_unstr_mesh_aux`.

Second, check rank-0 memory growth in the phase markers. Large jumps in `rsize` or `dsize`, especially before the last marker, may indicate memory pressure or OOM rather than an MPI deadlock.

Third, inspect the node-support diagnostics. Important red flags are:

- One or more nodes shared by an unusually large number of cells.
- Very large total candidate cell-pair counts.
- Very large cross-partition pair counts.
- A large max cross-partition pair count for a single node.

These would support the hypothesis that geometric complexity is causing pathological work in `all_cell_neighbors` or subsequent ghost/link construction.

Fourth, check rank-labeled stderr for a first failing rank. If one rank aborts or is OOM-killed before a collective, other ranks may appear hung in MPI. The first rank-specific error is more important than the final apparent hang point.

## Likely Interpretation

The leading hypothesis is not the viewfactor physics itself. The more likely failure area is unstructured distributed mesh construction after cell partitioning, especially:

- node partitioning and node-support construction,
- face numbering/partitioning,
- ghost-cell identification from `all_cell_neighbors`,
- memory pressure on the I/O/rank-0 process,
- pathological high-valence mesh nodes caused by geometric complexity.

When the customer output is available, compare the final marker, memory progression, and node-support summary against the code path in `src/truchas/distributed_mesh/unstr_mesh_factory.F90`.
