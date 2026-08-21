# Volume-tracking modernization plan

The minimal path is additive: preserve the tracker algorithms and legacy gold files while replacing only the test-harness infrastructure.

## 1. Remove direct TLS usage

Start with the `volume_tracking` sources:

- Remove `TLS_initialize`, `TLS_set_verbosity`, and the `truchas_env` setup from the four VOF programs.
- Continue initializing `env%simlog`, since the 2D mesh factory already uses it.
- Replace informational TLS calls with either `env%simlog%info` or no output.
- Replace impossible-condition `TLS_panic` calls with `INSIST`/`ASSERT`.
- Replace the `TLS_fatal_if_any` call in flux renormalization with an explicit collective check followed by `INSIST`.
- Remove unused logging imports.

This eliminates direct TLS use from the volume-tracking code. `f90_assert` still currently uses TLS internally, so eliminating that lower-level dependency entirely would be a separate utility change.

Leave the legacy timer calls alone initially; modernizing them through `env%timer` is independent of getting the tests operational.

## 2. Add a VOF VTKHDF writer

Add a small `vof_2d_vtkhdf_writer_type.F90` module, patterned after the existing flow and thermal writers.

It should write:

- the complete local mesh, including ghost cells and adjacent nodes;
- `GlobalIds`, `PedigreeIds`, and `vtkGhostType` for cells and nodes;
- a temporal multi-component cell field such as `volume-fraction`;
- optionally the interface normal field later.

The writer should take `env%comm` directly rather than using the legacy global `parallel_communication` communicator. It should also accept an output filename instead of hard-coding `out.vtkhdf`, although `out.vtkhdf` can remain the default.

`trunc_vol` does not need VTKHDF output; it is a low-level scalar geometry test.

## 3. Wire output into the existing drivers

Initially write the initial and final VOF states, or every timestep if the cost is negligible. The least invasive location is the existing `vof_2d_test_driver`, with the executable opening and closing the writer.

Keep the current text-reference comparison in place during this transition. That gives an independent check that adding VTKHDF output has not changed the algorithm.

## 4. Move the regression tests to Python/VTKHDF

The existing `TruchasVTKHDFData.py` reader already does the important work:

- reads VTKHDF through VTK;
- removes ghost cells;
- orders cells by pedigree IDs;
- provides cell centers and fields.

Add a Python test harness for the five cases. It can initially read the existing `circle*.txt` files as expected data while obtaining the computed result from VTKHDF.

Each case should run in a temporary directory so serial and four-process output files do not collide. The tests should compare:

- final volume fractions against the current reference;
- serial versus four-process volume fractions;
- optionally serial versus four-process cell centers.

Once these tests pass, the Fortran-side text comparison can be removed and the `circle*.txt` files can either remain as historical references or be retired later.

## 5. Simplify CMake

Update `volume_tracking/CMakeLists.txt` to:

- add the VTKHDF writer to `volume_tracking_2d`;
- find the Python interpreter;
- replace the direct executable/reference tests with Python-driven tests;
- retain `truncated_volumes` as a direct one-process test.

This gives the first milestone:

```text
existing algorithms
        ↓
TLS-free test programs
        ↓
VTKHDF output with partition metadata
        ↓
Python reader and serial/parallel regression tests
```

The axisymmetric cases require no separate output architecture. They use the same writer and tracker; only the mesh measures and geometric calculations differ.
