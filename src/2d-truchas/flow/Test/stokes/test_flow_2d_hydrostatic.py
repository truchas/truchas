#!/usr/bin/env python3

"""End-to-end hydrostatic test for a closed 2D flow domain."""

import pathlib
import sys

import numpy as np

from test_flow_2d_gravity import run_case


def check_hydrostatic(data, case):
    pressure_tol = 1.0e-11
    step = data.num_steps - 1
    if abs(data.time(step) - 20.0) > 1.0e-12:
        raise RuntimeError(f"{case}: final time is {data.time(step):g}, expected 20")
    centers = data.cell_centers(step)
    pressure = data.field(step, "pressure")
    velocity = data.field(step, "velocity")
    velocity_error = np.max(np.abs(velocity))
    if velocity_error > 1.0e-12:
        raise RuntimeError(f"{case}: hydrostatic velocity error {velocity_error:g}")
    expected = -centers[:, 1]
    pressure_error = np.max(np.abs((pressure - expected) - np.mean(pressure - expected)))
    if pressure_error > pressure_tol:
        raise RuntimeError(f"{case}: hydrostatic pressure error {pressure_error:g}")
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
        serial_results = check_hydrostatic(serial, "serial")
        parallel_results = check_hydrostatic(parallel, "parallel")
    except RuntimeError as error:
        print(f"FAIL: {error}")
        return 1

    for serial_result, parallel_result, name in zip(serial_results, parallel_results,
                                                      ("cell centers", "pressure", "velocity")):
        tolerance = 1.0e-11 if name == "pressure" else 1.0e-12
        if not np.allclose(serial_result, parallel_result, rtol=0.0, atol=tolerance):
            print(f"FAIL: serial and parallel {name} differ")
            return 1

    print("PASS: closed-box hydrostatic state is preserved")
    print(f"      serial artifacts: {serial_dir}")
    print(f"      parallel artifacts: {parallel_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
