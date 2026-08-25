#!/usr/bin/env python3

"""Analytic and MPI-consistency test for viscous channel flow with heat transfer."""

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, mpiexec):
    output_dir = Path(tempfile.mkdtemp(prefix="ns_ht_2d_channel_4p_"))
    command = [str(executable), "--output-dir", ".", "--force", str(input_file)]
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
        raise RuntimeError(f"ns_ht_2d returned {result.returncode} with 4 processes")
    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"ns_ht_2d did not produce {output_file}")
    return TruchasVTKHDFData(output_file), output_dir


def check_solution(data, angle):
    if data.num_steps != 2:
        raise RuntimeError(f"expected initial and final output states, found {data.num_steps}")
    if abs(data.time(0)) > 1.0e-14 or abs(data.time(1) - 50.0) > 1.0e-12:
        raise RuntimeError("incorrect VTKHDF output times")

    theta = np.deg2rad(angle)
    tangent = np.array([np.cos(theta), np.sin(theta), 0.0])
    normal = np.array([-np.sin(theta), np.cos(theta), 0.0])
    for step in range(data.num_steps):
        centers = data.cell_centers(step)
        temperature = data.field(step, "T")
        enthalpy = data.field(step, "H")
        expected_temperature = 1.0 + centers @ normal
        temperature_error = np.max(np.abs(temperature - expected_temperature))
        enthalpy_error = np.max(np.abs(enthalpy - expected_temperature))
        if temperature_error > 1.0e-10 or enthalpy_error > 1.0e-10:
            raise RuntimeError(
                f"step {step}: linear thermal errors are "
                f"temperature={temperature_error:g}, enthalpy={enthalpy_error:g}"
            )
        if not np.all(np.isfinite(data.field(step, "pressure"))):
            raise RuntimeError(f"step {step}: non-finite pressure output")

    pressure = data.field(1, "pressure")
    velocity = data.field(1, "velocity")
    normal_coordinate = centers @ normal
    expected_speed = 0.5 * normal_coordinate * (1.0 - normal_coordinate)
    expected_velocity = expected_speed[:, None] * tangent
    velocity_error = np.max(np.abs(velocity - expected_velocity))
    if velocity_error > 5.0e-3:
        raise RuntimeError(f"Poiseuille velocity error is {velocity_error:g}")
    return centers, pressure, velocity, temperature, enthalpy


def main():
    if len(sys.argv) != 5:
        print(f"usage: {sys.argv[0]} NS_HT_2D JSON_INPUT MPIEXEC ANGLE", file=sys.stderr)
        return 2

    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]
    angle = float(sys.argv[4])
    try:
        data, output_dir = run_case(executable, input_file, mpiexec)
        check_solution(data, angle)
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        return 1

    print("PASS: viscous channel matches the Poiseuille/linear-temperature solution")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
