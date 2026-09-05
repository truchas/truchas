#!/usr/bin/env python3

"""Four-process hydrostatic regression for mixed material-defined solid walls."""

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
    output_dir = Path(tempfile.mkdtemp(prefix="ns_2d_solid_wall_hydrostatic_mixed_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--simulation", "flow", "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        print(result.stdout, end="")
        print(f"FAIL: ns_2d returned {result.returncode}")
        return 1

    data = TruchasVTKHDFData(output_dir / "out.vtkhdf")
    if data.num_steps != 2:
        print(f"FAIL: found {data.num_steps} output states, expected 2")
        return 1

    initial_fluid = None
    initial_wall = None
    for step in range(data.num_steps):
        expected_time = 0.0 if step == 0 else 1.0
        if abs(data.time(step) - expected_time) > 1.0e-12:
            print(f"FAIL: output {step} has time {data.time(step):g}")
            return 1

        centers = data.cell_centers(step)
        fluid = data.field(step, "vf_fluid")
        wall = data.field(step, "vf_wall")
        mixed = (fluid > 1.0e-12) & (fluid < 1.0 - 1.0e-12)
        fluid_cells = fluid > 1.0e-12
        if not np.any(mixed):
            print(f"FAIL: step {step}: no mixed cells found")
            return 1
        if not np.allclose(fluid + wall, 1.0, rtol=0.0, atol=1.0e-12):
            print(f"FAIL: step {step}: material fractions do not sum to one")
            return 1
        if not np.allclose(fluid[~fluid_cells], 0.0, rtol=0.0, atol=1.0e-12):
            print(f"FAIL: step {step}: pure solid cells contain fluid")
            return 1

        if step == 0:
            initial_fluid = fluid.copy()
            initial_wall = wall.copy()
        elif not np.array_equal(fluid, initial_fluid) or not np.array_equal(wall, initial_wall):
            print("FAIL: material fractions changed in the hydrostatic state")
            return 1

        pressure = data.field(step, "pressure")
        velocity = data.field(step, "velocity")
        y = centers[fluid_cells, 1]
        pressure_error = np.max(np.abs(pressure[fluid_cells] - 2.0 * (1.0 - y)))
        velocity_error = np.max(np.abs(velocity[fluid_cells]))
        if pressure_error > 1.0e-10 or velocity_error > 1.0e-10:
            print(
                f"FAIL: step {step}: hydrostatic errors are "
                f"pressure={pressure_error:g}, velocity={velocity_error:g}"
            )
            return 1

    print("PASS: mixed solid-wall hydrostatic state is preserved")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
