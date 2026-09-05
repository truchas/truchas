#!/usr/bin/env python3

"""Four-process regression for a two-fluid isothermal hydrostatic column."""

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

    executable, input_file, mpiexec = map(Path, sys.argv[1:])
    output_dir = Path(tempfile.mkdtemp(prefix="ns_2d_multifluid_hydrostatic_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--simulation", "flow", "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"ns_2d returned {result.returncode} in {output_dir}")

    data = TruchasVTKHDFData(output_dir / "out.vtkhdf")
    if data.num_steps != 2:
        raise RuntimeError(f"expected two output states, found {data.num_steps}")

    for step, expected_time in enumerate((0.0, 0.1)):
        if abs(data.time(step) - expected_time) > 1.0e-12:
            raise RuntimeError(f"output {step} has time {data.time(step):g}")

        centers = data.cell_centers(step)
        y = centers[:, 1]
        light = data.field(step, "vf_light_liquid")
        heavy = data.field(step, "vf_heavy_liquid")
        pressure = data.field(step, "pressure")
        velocity = data.field(step, "velocity")

        expected_heavy = (y <= 0.5).astype(float)
        if np.max(np.abs(heavy - expected_heavy)) > 1.0e-12:
            raise RuntimeError(f"step {step}: incorrect heavy-fluid volume fractions")
        if np.max(np.abs(light + heavy - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: volume fractions do not sum to one")

        expected_pressure = np.where(y <= 0.5, 2.0*y, y + 0.5)
        pressure_error = np.max(np.abs(pressure - expected_pressure))
        velocity_error = np.max(np.abs(velocity[:, :2]))
        if pressure_error > 1.0e-10:
            raise RuntimeError(f"step {step}: pressure error={pressure_error:g}")
        if velocity_error > 1.0e-10:
            raise RuntimeError(f"step {step}: velocity error={velocity_error:g}")

    print("PASS: two-fluid isothermal hydrostatic pressure and zero velocity")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
