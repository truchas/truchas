#!/usr/bin/env python3

"""Analytic four-process test for two-fluid Stokes Couette flow."""

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
        print(f"usage: {sys.argv[0]} FLOW_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]
    output_dir = Path(tempfile.mkdtemp(prefix="stokes_2d_couette_viscosity_4p_"))
    result = subprocess.run(
        [mpiexec, "-n", "4", str(executable), "--simulation", "ns_2d", "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        print(result.stdout, end="")
        print(f"FAIL: flow_2d returned {result.returncode}")
        return 1

    data = TruchasVTKHDFData(output_dir / "out.vtkhdf")
    if data.num_steps != 2:
        print(f"FAIL: found {data.num_steps} output states, expected 2")
        return 1

    for step, expected_time in enumerate((0.0, 0.1)):
        if abs(data.time(step) - expected_time) > 1.0e-12:
            print(f"FAIL: output {step} has time {data.time(step):g}, expected {expected_time:g}")
            return 1

        centers = data.cell_centers(step)
        velocity = data.field(step, "velocity")
        lower = data.field(step, "vf_lower_fluid")
        upper = data.field(step, "vf_upper_fluid")
        y = centers[:, 1]
        expected_velocity = np.zeros_like(velocity)
        expected_velocity[:, 0] = np.where(y <= 0.5, 4.0 * y / 3.0, 2.0 * y / 3.0 + 1.0 / 3.0)
        velocity_error = np.max(np.abs(velocity - expected_velocity))
        pressure_error = np.max(np.abs(data.field(step, "pressure")))
        fraction_error = np.max(np.abs(lower + upper - 1.0))
        if velocity_error > 1.0e-9 or pressure_error > 1.0e-9 or fraction_error > 1.0e-12:
            print(
                f"FAIL: errors at output {step} are "
                f"velocity={velocity_error:g}, pressure={pressure_error:g}, "
                f"volume-fraction={fraction_error:g}"
            )
            return 1

    print("PASS: two-fluid Stokes Couette flow matches the analytic profile")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
