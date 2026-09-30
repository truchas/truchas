#!/usr/bin/env python3

"""Check NaN preservation and range handling in a flow VTKHDF file."""

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOHDF import vtkHDFReader
from vtkmodules.vtkFiltersCore import vtkAppendFilter


def fail(message):
    raise RuntimeError(message)


def read_grid(reader):
    reader.Update()
    data = reader.GetOutput()
    if data.IsA("vtkPartitionedDataSet"):
        append = vtkAppendFilter()
        for i in range(data.GetNumberOfPartitions()):
            append.AddInputData(data.GetPartition(i))
        append.Update()
        return append.GetOutput()
    return data


def check_range(array, expected, label):
    for range_name in ("GetRange", "GetFiniteRange"):
        value = np.asarray(getattr(array, range_name)(), dtype=float)
        if value.shape != (2,) or not np.isfinite(value).all():
            fail(f"{label}: {range_name} returned {value}")
        if not np.allclose(value, expected, rtol=0.0, atol=0.0):
            fail(f"{label}: {range_name} returned {value}, expected {expected}")


def main():
    if len(sys.argv) not in (3, 4):
        print(f"usage: {sys.argv[0]} TEST_EXECUTABLE MPIEXEC [NPROC]", file=sys.stderr)
        return 2

    executable = Path(sys.argv[1]).resolve()
    mpiexec = shutil.which(sys.argv[2]) or sys.argv[2]
    nproc = sys.argv[3] if len(sys.argv) == 4 else "1"
    output_dir = Path(tempfile.mkdtemp(prefix="flow_vtkhdf_nan_"))
    (output_dir / "disabled").mkdir()
    (output_dir / "mainline").mkdir()
    result = subprocess.run(
        [mpiexec, "-n", nproc, str(executable)],
        cwd=output_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        print(result.stdout, end="")
        fail(f"NaN writer executable returned {result.returncode}")

    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        fail(f"NaN writer did not produce {output_file}")

    reader = vtkHDFReader()
    reader.SetFileName(str(output_file))
    grid = read_grid(reader)
    if grid is None or grid.GetNumberOfCells() < 8:
        fail("VTK did not read the expected eight-cell unstructured grid")
    ids = np.asarray(vtk_to_numpy(grid.GetCellData().GetArray("GlobalCellIds")))
    if not np.array_equal(np.unique(ids), np.arange(1, 9)):
        fail(f"unexpected global cell IDs: {ids}")
    active = ids % 4 != 0
    if int(nproc) > 1:
        ghosts = vtk_to_numpy(grid.GetCellData().GetArray("vtkGhostType"))
        if not np.any(ghosts):
            fail("parallel output did not contain ghost cells to check")

    pressure = grid.GetCellData().GetArray("pressure")
    velocity = grid.GetCellData().GetArray("velocity")
    if pressure is None or velocity is None:
        fail("VTK did not expose pressure and velocity cell data")

    pressure_values = np.asarray(vtk_to_numpy(pressure))
    velocity_values = np.asarray(vtk_to_numpy(velocity))
    if pressure_values.shape != ids.shape or not np.isnan(pressure_values[~active]).all():
        fail(f"pressure did not preserve the quiet NaN: {pressure_values}")
    if not np.allclose(pressure_values[active], 1.0):
        fail(f"pressure finite value changed: {pressure_values}")
    if velocity_values.shape != (len(ids), 3) or not np.isnan(velocity_values[~active]).all():
        fail(f"velocity did not preserve the quiet NaNs: {velocity_values}")
    if not np.allclose(velocity_values[active], [0.25, 0.0, 0.0]):
        fail(f"velocity finite value changed: {velocity_values}")

    check_range(pressure, [1.0, 1.0], "pressure")
    compliance = grid.GetCellData().GetArray("void_compliance")
    if compliance is None:
        fail("VTK did not expose void_compliance cell data")
    values = np.asarray(vtk_to_numpy(compliance))
    if (values.shape != ids.shape or not np.all(values[active] == 0.125)
            or not np.isnan(values[~active]).all()):
        fail(f"unexpected void compliance values: {values}")
    check_range(compliance, [0.125, 0.125], "void_compliance")
    target = grid.GetCellData().GetArray("void_target_divergence")
    if target is None:
        fail("VTK did not expose void_target_divergence cell data")
    expected = np.array([np.nan, -0.5, 0.25, 0.0])[ids % 4]
    if not np.allclose(vtk_to_numpy(target), expected, rtol=0, atol=0, equal_nan=True):
        fail("incorrect signed targets, equilibrium zeros, or inactive NaNs (including ghosts)")
    check_range(target, [-0.5, 0.25], "void_target_divergence")
    for component, expected in enumerate(([0.25, 0.25], [0.0, 0.0], [0.0, 0.0])):
        for range_name in ("GetRange", "GetFiniteRange"):
            value = np.asarray(getattr(velocity, range_name)(component), dtype=float)
            if value.shape != (2,) or not np.isfinite(value).all():
                fail(f"velocity component {component}: {range_name} returned {value}")
            if not np.allclose(value, expected, rtol=0.0, atol=0.0):
                fail(
                    f"velocity component {component}: {range_name} returned {value}, "
                    f"expected {expected}"
                )

    reader.SetStep(1)
    grid = read_grid(reader)
    compliance = grid.GetCellData().GetArray("void_compliance")
    if not np.isnan(vtk_to_numpy(compliance)).all():
        fail("absent compliance did not produce all NaNs at the next output time")
    target = grid.GetCellData().GetArray("void_target_divergence")
    if not np.isnan(vtk_to_numpy(target)).all():
        fail("absent target divergence did not produce all NaNs at the next output time")

    disabled_reader = vtkHDFReader()
    disabled_reader.SetFileName(str(output_dir / "disabled" / "out.vtkhdf"))
    disabled_grid = read_grid(disabled_reader)
    for name in ("void_compliance", "void_target_divergence"):
        if disabled_grid.GetCellData().GetArray(name) is not None:
            fail(f"disabled compliance still published {name}")
    if disabled_grid.GetCellData().GetArray("pressure") is None:
        fail("disabled compliance output lost pressure data")

    mainline_reader = vtkHDFReader()
    mainline_reader.SetFileName(str(output_dir / "mainline" / "out.vtkhdf"))
    mainline_grid = read_grid(mainline_reader)
    if mainline_grid.GetCellData().GetArray("void_compliance") is not None:
        fail("mainline output published a compliance coefficient")
    target = mainline_grid.GetCellData().GetArray("void_target_divergence")
    if target is None:
        fail("mainline output omitted the target divergence")
    mainline_ids = vtk_to_numpy(mainline_grid.GetCellData().GetArray("GlobalCellIds"))
    expected = np.array([np.nan, -0.5, -0.25, 0.0])[mainline_ids % 4]
    if not np.allclose(vtk_to_numpy(target), expected, rtol=0, atol=0, equal_nan=True):
        fail("incorrect mainline targets or inactive NaNs (including ghosts)")

    print("PASS: quiet NaNs survive flow VTKHDF output and are ignored by VTK ranges")
    print(f"      artifact: {output_file}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
