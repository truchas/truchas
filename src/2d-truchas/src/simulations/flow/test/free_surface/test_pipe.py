#!/usr/bin/env python3

"""Four-process analytic regressions for one-dimensional fluid/VOID pipes."""

import pathlib
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[7]
sys.path.insert(0, str(source_root / "src/2d-truchas/python"))
from TruchasVTKHDFData import TruchasVTKHDFData


def expected_water_volume_fraction(case, time):
    if case == "plug":
        lower, upper = 0.5 + time, 1.5 + time
    elif case == "fill":
        lower, upper = 0.0, 0.5 + time
    elif case == "drain":
        lower, upper = 0.0, 4.5 - time
    else:
        raise ValueError(f"unknown pipe case {case!r}")

    nx, ny = (25, 5) if case == "plug" else (15, 5)
    dx = 5.0 / nx
    values = [
        max(0.0, min(upper, (i + 1) * dx) - max(lower, i * dx)) / dx
        for i in range(nx)
    ]
    return np.tile(values, ny)


def run_case(executable, input_file, mpiexec):
    output_dir = pathlib.Path(tempfile.mkdtemp(prefix="flow_free_surface_pipe_4p_"))
    command = [
        mpiexec,
        "-n",
        "4",
        str(executable),
        "--simulation",
        "flow",
        "--output-dir",
        ".",
        "--force",
        str(input_file),
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
        raise RuntimeError(f"flow returned {result.returncode}")

    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"flow did not produce {output_file}")
    return TruchasVTKHDFData(output_file), output_dir


def check_solution(data, case):
    expected_times = (0.0, 1.0, 3.0)
    if data.num_steps != len(expected_times):
        raise RuntimeError(
            f"expected {len(expected_times)} output states, got {data.num_steps}"
        )

    for step, time in enumerate(expected_times):
        observed_time = data.time(step)
        if abs(observed_time - time) > 1.0e-14:
            raise RuntimeError(f"step {step}: time={observed_time:g}, expected {time:g}")

        water = data.field(step, "vf_water")
        water_error = np.sum(
            np.abs(water - expected_water_volume_fraction(case, time))
        )
        if water_error > 1.0e-2:
            raise RuntimeError(
                f"step {step}: water-volume-fraction L1 error={water_error:g}"
            )

        velocity = data.field(step, "velocity")
        active = np.isfinite(velocity[:, 0])
        expected_velocity = 1.0 if case != "drain" else -1.0
        velocity_error = np.max(
            np.abs(velocity[active, :2] - [expected_velocity, 0.0])
        )
        if velocity_error > 1.0e-12:
            raise RuntimeError(f"step {step}: velocity error={velocity_error:g}")

        pressure = data.field(step, "pressure")
        active = np.isfinite(pressure)
        pressure_error = np.max(np.abs(pressure[active]))
        if pressure_error > 1.0e-12:
            raise RuntimeError(f"step {step}: pressure error={pressure_error:g}")


def main():
    if len(sys.argv) != 5:
        print(
            f"usage: {sys.argv[0]} TRUCHAS_2D JSON_INPUT MPIEXEC PIPE_CASE",
            file=sys.stderr,
        )
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    case = sys.argv[4]
    try:
        data, output_dir = run_case(executable, input_file, sys.argv[3])
        check_solution(data, case)
    except (RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        return 1

    print(f"PASS: free-surface {case} pipe matches the analytic solution")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
