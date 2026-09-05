#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, mpiexec):
    run_root = Path(tempfile.mkdtemp(prefix="ns_ht_2d_hydrostatic_4p_"))
    command = [mpiexec, "-n", "4", str(executable), "--simulation", "ns_ht_2d", str(input_file)]
    result = subprocess.run(
        command, cwd=run_root, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False,
    )
    (run_root / "stdout.txt").write_text(result.stdout)
    if result.returncode != 0:
        raise AssertionError(f"ns_ht_2d returned {result.returncode} in {run_root}\n{result.stdout}")
    return run_root, TruchasVTKHDFData(run_root / "input" / "out.vtkhdf")


def exact_pressure(y):
    rho = 1.0
    alpha = 0.5
    tref = 0.0
    acceleration = 1.0
    top_temperature = 1.0
    bottom_temperature = 0.0
    pressure_slope = 1.0 - alpha * (top_temperature - tref)
    return rho * acceleration * (
        pressure_slope * y
        + alpha * (top_temperature - bottom_temperature) * y**2 / 2.0
    )


def check_case(data):
    pressure_tol = 1.0e-3
    velocity_tol = 1.0e-4

    if data.num_steps != 2:
        raise AssertionError(f"expected two output states, found {data.num_steps}")
    for step, expected in enumerate((0.0, 1.0)):
        if abs(data.time(step) - expected) > 1.0e-10:
            raise AssertionError(f"output {step} has time {data.time(step):g}, expected {expected:g}")

    for step in (0, 1):
        centers = data.cell_centers(step)
        y = centers[:, 1]
        temperature = data.field(step, "T")
        pressure = data.field(step, "pressure")
        velocity = data.field(step, "velocity")
        if not (np.isfinite(temperature).all() and np.isfinite(pressure).all() and np.isfinite(velocity).all()):
            raise AssertionError("hydrostatic solution contains non-finite values")

        temperature_error = np.max(np.abs(temperature - (1.0 - y)))
        if temperature_error > 1.0e-7:
            raise AssertionError(f"temperature error at t={data.time(step):g} is {temperature_error:g}")

        pressure_error = np.max(np.abs(pressure - exact_pressure(y)))
        if pressure_error > pressure_tol:
            raise AssertionError(f"pressure error at t={data.time(step):g} is {pressure_error:g}")

        velocity_error = np.max(np.abs(velocity))
        if velocity_error > velocity_tol:
            raise AssertionError(f"velocity error at t={data.time(step):g} is {velocity_error:g}")

    return np.max(np.abs(data.field(1, "velocity")))


def test():
    if len(sys.argv) != 3:
        raise AssertionError(f"usage: {sys.argv[0]} TRUCHAS_2D MPIEXEC")
    executable = Path(sys.argv[1]).resolve()
    input_file = Path(__file__).with_name("input.json").resolve()
    run_root, data = run_case(executable, input_file, sys.argv[2])
    max_velocity = check_case(data)
    print(f"PASS: hydrostatic; final maximum velocity {max_velocity:g}; artifacts: {run_root}")


if __name__ == "__main__":
    try:
        test()
    except (AssertionError, OSError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
