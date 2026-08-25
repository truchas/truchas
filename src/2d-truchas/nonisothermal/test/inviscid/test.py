#!/usr/bin/env python3

"""Analytic four-process test for nonisothermal inviscid plug flow."""

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, mpiexec):
    output_dir = Path(tempfile.mkdtemp(prefix="ns_ht_2d_inviscid_4p_"))
    command = [
        mpiexec, "-n", "4", str(executable), "--output-dir", ".", "--force", str(input_file)
    ]
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
        raise RuntimeError(f"ns_ht_2d returned {result.returncode}")
    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"ns_ht_2d did not produce {output_file}")
    return TruchasVTKHDFData(output_file), output_dir


def check_solution(data, angle):
    if data.num_steps < 2:
        raise RuntimeError("expected initial and final output states")
    if abs(data.time(0)) > 1.0e-14:
        raise RuntimeError(f"initial output time is {data.time(0):g}, expected 0")
    if abs(data.time(data.num_steps - 1) - 0.01) > 1.0e-14:
        raise RuntimeError(
            f"final output time is {data.time(data.num_steps - 1):g}, expected 0.01"
        )

    theta = np.deg2rad(angle)
    analytic_velocity = np.array([np.cos(theta), np.sin(theta), 0.0])
    for step in range(data.num_steps):
        centers = data.cell_centers(step)
        pressure = data.field(step, "pressure")
        velocity = data.field(step, "velocity")
        temperature = data.field(step, "T")
        enthalpy = data.field(step, "H")
        expected_velocity = np.zeros_like(velocity)
        expected_velocity[:, :2] = analytic_velocity[:2]
        expected_temperature = 1.0 - np.sin(theta)*centers[:, 0] + np.cos(theta)*centers[:, 1]
        errors = {
            "velocity": np.max(np.abs(velocity - expected_velocity)),
            "pressure": np.max(np.abs(pressure)),
            "temperature": np.max(np.abs(temperature - expected_temperature)),
            "enthalpy": np.max(np.abs(enthalpy - expected_temperature)),
        }
        error = max(errors.values())
        if error > 1.0e-10:
            details = ", ".join(f"{name}={value:g}" for name, value in errors.items())
            raise RuntimeError(f"step {step}: analytic solution errors are {details}")


def main():
    if len(sys.argv) not in (4, 5):
        print(f"usage: {sys.argv[0]} NS_HT_2D JSON_INPUT MPIEXEC [ANGLE]", file=sys.stderr)
        return 2

    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]
    angle = float(sys.argv[4]) if len(sys.argv) == 5 else 0.0
    try:
        data, output_dir = run_case(executable, input_file, mpiexec)
        check_solution(data, angle)
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        return 1

    print("PASS: nonisothermal inviscid channel matches the analytic solution")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
