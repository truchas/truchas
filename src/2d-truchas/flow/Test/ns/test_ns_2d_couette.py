#!/usr/bin/env python3

"""Analytic four-process test for mixed-BC single-fluid Couette flow."""

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def main():
    if len(sys.argv) not in (4, 5):
        print(f"usage: {sys.argv[0]} TRUCHAS_2D JSON_INPUT MPIEXEC [ROTATION_ANGLE]", file=sys.stderr)
        return 2

    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]
    angle = float(sys.argv[4]) if len(sys.argv) == 5 else 0.0
    output_dir = Path(tempfile.mkdtemp(prefix="ns_2d_couette_4p_"))
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
    expected_times = (0.0, 50.0)
    if data.num_steps != len(expected_times):
        print(f"FAIL: found {data.num_steps} output states, expected {len(expected_times)}")
        return 1

    for step, expected_time in enumerate(expected_times):
        if abs(data.time(step) - expected_time) > 1.0e-12:
            print(f"FAIL: output {step} has time {data.time(step):g}, expected {expected_time:g}")
            return 1

    theta = np.deg2rad(angle)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    for step in range(data.num_steps):
        centers = data.cell_centers(step)
        velocity = data.field(step, "velocity")
        expected_velocity = np.zeros_like(velocity)
        local_y = -sin_theta * centers[:, 0] + cos_theta * centers[:, 1]
        expected_velocity[:, 0] = cos_theta * local_y
        expected_velocity[:, 1] = sin_theta * local_y
        velocity_error = np.max(np.abs(velocity - expected_velocity))
        pressure_error = np.max(np.abs(data.field(step, "pressure")))
        pressure_tolerance = 2.0e-10 if step == 0 else 1.0e-10
        if velocity_error > 1.0e-10 or pressure_error > pressure_tolerance:
            print(
                f"FAIL: Couette errors at output {step} are "
                f"velocity={velocity_error:g}, pressure={pressure_error:g}"
            )
            return 1

    print("PASS: mixed-BC Couette flow matches the analytic profile")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
