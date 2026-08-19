#!/usr/bin/env python3

"""End-to-end body-acceleration test for the JSON-driven 2D flow application."""

import pathlib
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[4]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, nproc, mpiexec):
    output_dir = pathlib.Path(tempfile.mkdtemp(prefix=f"flow_2d_gravity_{nproc}p_"))
    command = [str(executable), str(input_file)]
    if nproc > 1:
        command = [mpiexec, "-n", str(nproc)] + command
    result = subprocess.run(command, cwd=output_dir, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            check=False)
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"flow_2d returned {result.returncode} with {nproc} processes")
    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"flow_2d did not produce {output_file}")
    return TruchasVTKHDFData(output_file), output_dir


def check_body_balance(data, case):
    step = data.num_steps - 1
    if abs(data.time(step) - 50.0) > 1.0e-12:
        raise RuntimeError(f"{case}: final time is {data.time(step):g}, expected 50")
    centers = data.cell_centers(step)
    pressure = data.field(step, "pressure")
    velocity = data.field(step, "velocity")

    velocity_error = np.max(np.abs(velocity))
    if velocity_error > 5.0e-4:
        raise RuntimeError(f"{case}: pressure-balanced velocity error {velocity_error:g}")

    expected = 1.0 - centers[:, 0]
    pressure_error = np.max(np.abs((pressure - expected) - np.mean(pressure - expected)))
    if pressure_error > 1.0e-2:
        raise RuntimeError(f"{case}: balanced pressure error {pressure_error:g}")
    return centers, pressure, velocity


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} FLOW_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]

    try:
        serial, serial_dir = run_case(executable, input_file, 1, mpiexec)
        parallel, parallel_dir = run_case(executable, input_file, 4, mpiexec)
        serial_centers, serial_pressure, serial_velocity = check_body_balance(serial, "serial")
        parallel_centers, parallel_pressure, parallel_velocity = check_body_balance(parallel, "parallel")
    except RuntimeError as error:
        print(f"FAIL: {error}")
        return 1

    if not np.allclose(serial_centers, parallel_centers, rtol=0.0, atol=1.0e-14):
        print("FAIL: serial and parallel cell centers differ")
        return 1
    if not np.allclose(serial_pressure, parallel_pressure, rtol=0.0, atol=1.0e-10):
        print("FAIL: serial and parallel pressure differs")
        return 1
    if not np.allclose(serial_velocity, parallel_velocity, rtol=0.0, atol=1.0e-10):
        print("FAIL: serial and parallel velocity differs")
        return 1

    print("PASS: pressure gradient balances body acceleration")
    print(f"      serial artifacts: {serial_dir}")
    print(f"      parallel artifacts: {parallel_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
