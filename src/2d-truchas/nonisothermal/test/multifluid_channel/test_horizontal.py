#!/usr/bin/env python3

"""Analytic four-process regression for two-fluid inviscid plug flow."""

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
    output_dir = Path(tempfile.mkdtemp(prefix="ns_ht_2d_multifluid_channel_4p_"))
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
    initial_lower = None
    for step, expected_time in enumerate(expected_times):
        if abs(data.time(step) - expected_time) > 1.0e-14:
            raise RuntimeError(f"output {step} has time {data.time(step):g}")
        lower = data.field(step, "vf_lower_liquid")
        upper = data.field(step, "vf_upper_liquid")
        centers = data.cell_centers(step)
        temp = data.field(step, "T")
        enthalpy = data.field(step, "H")
        velocity = data.field(step, "velocity")
        if np.max(np.abs(lower + upper - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: material fractions do not sum to one")
        if np.min(lower) < -1.0e-12 or np.min(upper) < -1.0e-12:
            raise RuntimeError(f"step {step}: negative volume fraction")
        if step == 0:
            initial_lower = lower.copy()
        distribution_error = np.max(np.abs(lower - initial_lower))
        temp_error = np.max(np.abs(temp - (1.0 + centers[:, 1])))
        enthalpy_error = np.max(
            np.abs(enthalpy - (lower + 4.0*upper)*(1.0 + centers[:, 1]))
        )
        velocity_error = np.max(np.abs(velocity[:, :2] - [1.0, 0.0]))
        if max(distribution_error, temp_error, enthalpy_error, velocity_error) > 2.0e-9:
            raise RuntimeError(
                f"step {step}: distribution={distribution_error:g}, "
                f"temperature={temp_error:g}, enthalpy={enthalpy_error:g}, "
                f"velocity={velocity_error:g}"
            )
        if step == 0 and not np.any((lower > 1.0e-8) & (lower < 1.0 - 1.0e-8)):
            raise RuntimeError("interface did not create mixed cells")

    print("PASS: two-fluid plug flow preserves the analytic distribution and thermal state")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
