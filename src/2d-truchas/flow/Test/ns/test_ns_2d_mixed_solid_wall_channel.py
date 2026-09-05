#!/usr/bin/env python3

"""Four-process regression for an inviscid channel with mixed solid cells."""

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
    output_dir = Path(tempfile.mkdtemp(prefix="ns_2d_mixed_solid_wall_channel_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--simulation", "flow", "--output-dir", ".", "--force", str(input_file)],
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
    if abs(data.time(0)) > 1.0e-14 or abs(data.time(1) - 0.1) > 1.0e-14:
        print(f"FAIL: output times are {data.time(0):g}, {data.time(1):g}")
        return 1

    initial_fluid = None
    initial_wall = None
    for step in range(data.num_steps):
        centers = data.cell_centers(step)
        fluid = data.field(step, "vf_fluid")
        wall = data.field(step, "vf_wall")
        mixed = (fluid > 1.0e-12) & (fluid < 1.0 - 1.0e-12)
        pure_fluid = fluid > 1.0 - 1.0e-12
        pure_solid = fluid < 1.0e-12

        if mixed.sum() != 32:
            print(f"FAIL: step {step}: found {mixed.sum()} mixed cells, expected 32")
            return 1
        if not np.allclose(fluid + wall, 1.0, rtol=0.0, atol=1.0e-12):
            print(f"FAIL: step {step}: material fractions do not sum to one")
            return 1
        if not np.all(centers[fluid > 0.0, 1] > 0.25) or not np.all(centers[fluid > 0.0, 1] < 0.75):
            print(f"FAIL: step {step}: fluid cells are not confined to the channel")
            return 1
        if not np.allclose(fluid[pure_solid], 0.0, rtol=0.0, atol=1.0e-12):
            print(f"FAIL: step {step}: pure solid cells contain fluid")
            return 1

        if step == 0:
            initial_fluid = fluid.copy()
            initial_wall = wall.copy()
        elif not np.array_equal(fluid, initial_fluid) or not np.array_equal(wall, initial_wall):
            print("FAIL: material fractions changed during stationary-wall transport")
            return 1

        pressure = data.field(step, "pressure")[fluid > 0.0]
        velocity = data.field(step, "velocity")[fluid > 0.0]
        x = centers[fluid > 0.0, 0]
        expected_pressure = 1.0 - x
        expected_velocity = np.zeros_like(x) if step == 0 else np.full_like(x, 0.1)
        pressure_error = np.max(np.abs(pressure - expected_pressure))
        velocity_error = np.max(np.abs(velocity[:, 0] - expected_velocity))
        transverse_error = np.max(np.abs(velocity[:, 1]))
        if pressure_error > 1.0e-10 or velocity_error > 1.3e-3 or transverse_error > 1.0e-10:
            print(
                f"FAIL: step {step}: channel errors are "
                f"pressure={pressure_error:g}, velocity={velocity_error:g}, "
                f"transverse={transverse_error:g}"
            )
            return 1

    print("PASS: inviscid channel with mixed SOLID walls is conserved and analytic")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
