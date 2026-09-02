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


def fail(message):
    raise RuntimeError(message)


def check_range(array, expected, label):
    for range_name in ("GetRange", "GetFiniteRange"):
        value = np.asarray(getattr(array, range_name)(), dtype=float)
        if value.shape != (2,) or not np.isfinite(value).all():
            fail(f"{label}: {range_name} returned {value}")
        if not np.allclose(value, expected, rtol=0.0, atol=0.0):
            fail(f"{label}: {range_name} returned {value}, expected {expected}")


def main():
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} TEST_EXECUTABLE MPIEXEC", file=sys.stderr)
        return 2

    executable = Path(sys.argv[1]).resolve()
    mpiexec = shutil.which(sys.argv[2]) or sys.argv[2]
    output_dir = Path(tempfile.mkdtemp(prefix="flow_2d_vtkhdf_nan_"))
    result = subprocess.run(
        [mpiexec, "-n", "1", str(executable)],
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
    reader.Update()
    grid = reader.GetOutput()
    if grid is None or grid.GetNumberOfCells() != 2:
        fail("VTK did not read the expected two-cell unstructured grid")

    pressure = grid.GetCellData().GetArray("pressure")
    velocity = grid.GetCellData().GetArray("velocity")
    if pressure is None or velocity is None:
        fail("VTK did not expose pressure and velocity cell data")

    pressure_values = np.asarray(vtk_to_numpy(pressure))
    velocity_values = np.asarray(vtk_to_numpy(velocity))
    if pressure_values.shape != (2,) or not np.isnan(pressure_values[1]):
        fail(f"pressure did not preserve the quiet NaN: {pressure_values}")
    if not np.allclose(pressure_values[0], 1.0):
        fail(f"pressure finite value changed: {pressure_values}")
    if velocity_values.shape != (2, 3) or not np.isnan(velocity_values[1, :2]).all():
        fail(f"velocity did not preserve the quiet NaNs: {velocity_values}")
    if not np.allclose(velocity_values[0], [0.25, 0.0, 0.0]):
        fail(f"velocity finite value changed: {velocity_values}")

    check_range(pressure, [1.0, 1.0], "pressure")
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

    print("PASS: quiet NaNs survive flow VTKHDF output and are ignored by VTK ranges")
    print(f"      artifact: {output_file}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
