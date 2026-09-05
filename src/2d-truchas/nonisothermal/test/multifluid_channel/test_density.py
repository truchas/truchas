#!/usr/bin/env python3

"""Analytic regression for translating a discontinuous-density plug flow."""

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
    output_dir = Path(tempfile.mkdtemp(prefix="ns_ht_2d_multifluid_density_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--simulation", "flow_thermal", "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"ns_ht_2d returned {result.returncode}")

    data = TruchasVTKHDFData(output_dir / "out.vtkhdf")
    expected_times = (0.0, 0.0025, 0.005, 0.0075, 0.01)
    if data.num_steps != len(expected_times):
        raise RuntimeError("incorrect output times")

    dx = 1.0 / 16.0
    volume_fraction_tolerance = 5.0e-5
    for step, expected_time in enumerate(expected_times):
        if abs(data.time(step) - expected_time) > 1.0e-14:
            raise RuntimeError(f"output {step} has time {data.time(step):g}")
        light = data.field(step, "vf_light_liquid")
        heavy = data.field(step, "vf_heavy_liquid")
        centers = data.cell_centers(step)
        interface = 0.47 + expected_time
        expected_light = np.clip(
            (interface - (centers[:, 0] - 0.5*dx)) / dx, 0.0, 1.0
        )
        distribution_error = np.max(np.abs(light - expected_light))
        if distribution_error > volume_fraction_tolerance:
            raise RuntimeError(f"step {step}: distribution error={distribution_error:g}")
        if np.min(light) < -1.0e-12 or np.min(heavy) < -1.0e-12:
            raise RuntimeError(f"step {step}: negative volume fraction")
        if np.max(np.abs(light + heavy - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: volume fractions do not sum to one")

        velocity = data.field(step, "velocity")
        pressure = data.field(step, "pressure")
        temp = data.field(step, "T")
        enthalpy = data.field(step, "H")
        velocity_error = np.max(np.abs(velocity[:, :2] - [1.0, 0.0]))
        pressure_error = np.max(np.abs(pressure))
        enthalpy_error = np.max(np.abs(enthalpy - (light + 2.0*heavy)))
        temp_error = np.max(np.abs(temp - 1.0))
        if max(velocity_error, pressure_error, enthalpy_error, temp_error) > 2.0e-9:
            raise RuntimeError(
                f"step {step}: velocity={velocity_error:g}, pressure={pressure_error:g}, "
                f"enthalpy={enthalpy_error:g}, temperature={temp_error:g}"
            )

    print("PASS: variable-density plug flow preserves the analytic coupled state")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
