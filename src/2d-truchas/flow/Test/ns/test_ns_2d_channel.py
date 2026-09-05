#!/usr/bin/env python3

"""End-to-end analytic test for the JSON-driven 2D NS application."""

import pathlib
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, mpiexec):
    output_dir = pathlib.Path(tempfile.mkdtemp(prefix="ns_2d_channel_4p_"))
    command = [str(executable), "--simulation", "flow", "--output-dir", ".", "--force", str(input_file)]
    command = [mpiexec, "-n", "4"] + command
    result = subprocess.run(
        command,
        cwd=output_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"ns_2d returned {result.returncode} with 4 processes")
    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"ns_2d did not produce {output_file}")
    return TruchasVTKHDFData(output_file), output_dir


def check_profile(data, case, angle=30.0, tolerance=5.0e-3):
    step = data.num_steps - 1
    if abs(data.time(step) - 50.0) > 1.0e-12:
        raise RuntimeError(f"{case}: final time is {data.time(step):g}, expected 50")
    centers = data.cell_centers(step)
    velocity = data.field(step, "velocity")
    angle = np.deg2rad(angle)
    tangent = np.array([np.cos(angle), np.sin(angle), 0.0])
    normal = np.array([-np.sin(angle), np.cos(angle), 0.0])
    normal_coordinate = centers @ normal
    expected = 0.5 * normal_coordinate * (1.0 - normal_coordinate)
    expected_velocity = expected[:, None] * tangent
    error = np.max(np.abs(velocity - expected_velocity))
    if error > tolerance:
        raise RuntimeError(f"{case}: Poiseuille profile error {error:g}")
    return centers, velocity


def main():
    if len(sys.argv) not in (4, 5):
        print(f"usage: {sys.argv[0]} TRUCHAS_2D JSON_INPUT MPIEXEC [ANGLE]", file=sys.stderr)
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]
    angle = float(sys.argv[4]) if len(sys.argv) == 5 else 30.0

    try:
        data, output_dir = run_case(executable, input_file, mpiexec)
        tolerance = 1.0e-2 if "noisy" in input_file.stem else 5.0e-3
        check_profile(data, "4-process", angle, tolerance)
    except RuntimeError as error:
        print(f"FAIL: {error}")
        return 1

    print("PASS: integrated NS pressure channel matches the Poiseuille profile")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
