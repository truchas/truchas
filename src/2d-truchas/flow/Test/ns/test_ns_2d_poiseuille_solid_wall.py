#!/usr/bin/env python3

"""Analytic four-process test for viscous flow bounded by solid cells."""

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} TRUCHAS_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]
    output_dir = Path(tempfile.mkdtemp(prefix="ns_2d_poiseuille_solid_wall_4p_"))
    result = subprocess.run(
        [mpiexec, "-n", "4", str(executable), "--simulation", "flow", "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        print(result.stdout, end="")
        print(f"FAIL: ns_2d returned {result.returncode}")
        return 1

    data = TruchasVTKHDFData(output_dir / "out.vtkhdf")
    if data.num_steps != 2:
        print(f"FAIL: found {data.num_steps} output states, expected 2")
        return 1
    if abs(data.time(0)) > 1.0e-14 or abs(data.time(1) - 50.0) > 1.0e-12:
        print(f"FAIL: output times are {data.time(0):g}, {data.time(1):g}")
        return 1

    centers = data.cell_centers(1)
    fluid = data.field(1, "vf_fluid")
    wall = data.field(1, "vf_wall")
    fluid_cells = fluid > 0.5
    if not np.allclose(fluid + wall, 1.0, rtol=0.0, atol=1.0e-12):
        print("FAIL: material fractions do not sum to one")
        return 1
    if not np.all((centers[fluid_cells, 1] > 0.25) & (centers[fluid_cells, 1] < 0.75)):
        print("FAIL: fluid cells are not confined to the channel")
        return 1
    if not np.allclose(fluid[~fluid_cells], 0.0, rtol=0.0, atol=1.0e-12):
        print("FAIL: non-fluid cells contain fluid")
        return 1

    pressure = data.field(1, "pressure")[fluid_cells]
    velocity = data.field(1, "velocity")[fluid_cells]
    x = centers[fluid_cells, 0]
    y = centers[fluid_cells, 1]
    expected_pressure = 1.0 - x
    expected_velocity = 0.5 * (y - 0.25) * (0.75 - y)
    pressure_error = np.max(np.abs(pressure - expected_pressure))
    velocity_error = np.max(np.abs(velocity[:, 0] - expected_velocity))
    transverse_error = np.max(np.abs(velocity[:, 1]))
    if pressure_error > 1.0e-10 or velocity_error > 1.0e-3 or transverse_error > 1.0e-10:
        print(
            f"FAIL: channel errors are pressure={pressure_error:g}, "
            f"velocity={velocity_error:g}, transverse={transverse_error:g}"
        )
        return 1

    print("PASS: solid-wall Poiseuille flow matches the analytic solution")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
