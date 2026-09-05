#!/usr/bin/env python3

"""Four-process regression for a two-fluid thermally buoyant column."""

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
    output_dir = Path(tempfile.mkdtemp(prefix="ns_ht_2d_multifluid_hydrostatic_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--simulation", "ns_ht_2d", str(input_file)],
        cwd=output_dir, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False)
    (output_dir / "stdout.txt").write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"ns_ht_2d returned {result.returncode} in {output_dir}\n{result.stdout}")

    data = TruchasVTKHDFData(output_dir / "input_hydrostatic" / "out.vtkhdf")
    expected_times = (0.0, 0.1)
    if data.num_steps != len(expected_times):
        raise RuntimeError(f"expected two output states, found {data.num_steps}")

    for step, expected_time in enumerate(expected_times):
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

        alpha = 0.9
        pressure_at_interface = 2.0*((1.0 - alpha)*0.5 + alpha*0.5**2/2.0)
        expected_pressure = np.where(
            y <= 0.5,
            2.0*((1.0 - alpha)*y + alpha*y**2/2.0),
            pressure_at_interface + (1.0 - alpha)*(y - 0.5) +
            alpha*(y**2 - 0.5**2)/2.0)
        pressure_error = np.max(np.abs(pressure - expected_pressure))
        temperature_error = np.max(np.abs(data.field(step, "T") - (1.0 - y)))
        if temperature_error > 1.0e-10:
            raise RuntimeError(f"step {step}: temperature error={temperature_error:g}")
        velocity_error = np.max(np.abs(velocity[:, :2]))
        if pressure_error > 1.0e-3:
            raise RuntimeError(f"step {step}: pressure error={pressure_error:g}")
        if velocity_error > 1.0e-10:
            raise RuntimeError(f"step {step}: velocity error={velocity_error:g}")

    print("PASS: two-fluid thermally buoyant hydrostatic pressure and zero velocity")
    print(f"      artifacts: {output_dir}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
