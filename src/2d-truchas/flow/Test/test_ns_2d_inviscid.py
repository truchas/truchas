#!/usr/bin/env python3

"""Analytic four-process test for 2D inviscid plug flow."""

import pathlib
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[4]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, nproc, mpiexec):
    output_dir = pathlib.Path(
        tempfile.mkdtemp(prefix=f"ns_2d_inviscid_channel_{nproc}p_")
    )
    command = [str(executable), "--output-dir", ".", "--force", str(input_file)]
    if nproc > 1:
        command = [mpiexec, "-n", str(nproc)] + command
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
        raise RuntimeError(f"ns_2d returned {result.returncode} with {nproc} processes")
    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"ns_2d did not produce {output_file}")
    return TruchasVTKHDFData(output_file), output_dir


def check_solution(data, analytic_velocity):
    if data.num_steps < 2:
        raise RuntimeError("expected initial and final output states")
    if abs(data.time(0)) > 1.0e-14:
        raise RuntimeError(f"initial output time is {data.time(0):g}, expected 0")
    if abs(data.time(data.num_steps - 1) - 0.01) > 1.0e-14:
        raise RuntimeError(
            f"final output time is {data.time(data.num_steps - 1):g}, expected 0.01"
        )

    for step in range(data.num_steps):
        pressure = data.field(step, "pressure")
        velocity = data.field(step, "velocity")
        expected = np.zeros_like(velocity)
        expected[:, :2] = analytic_velocity[:2]
        velocity_error = np.max(np.abs(velocity - expected))
        pressure_error = np.max(np.abs(pressure))
        if velocity_error > 1.0e-12 or pressure_error > 1.0e-12:
            raise RuntimeError(
                f"step {step}: plug-flow errors are "
                f"velocity={velocity_error:g}, pressure={pressure_error:g}"
            )


def main():
    if len(sys.argv) not in (4, 6):
        print(f"usage: {sys.argv[0]} NS_2D JSON_INPUT MPIEXEC [VX VY]", file=sys.stderr)
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]
    expected_velocity = np.array([1.0, 0.0, 0.0])
    if len(sys.argv) == 6:
        expected_velocity[:2] = [float(sys.argv[4]), float(sys.argv[5])]

    try:
        data, output_dir = run_case(executable, input_file, 4, mpiexec)
        check_solution(data, expected_velocity)
    except RuntimeError as error:
        print(f"FAIL: {error}")
        return 1

    print("PASS: inviscid channel matches the analytic plug-flow solution")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
