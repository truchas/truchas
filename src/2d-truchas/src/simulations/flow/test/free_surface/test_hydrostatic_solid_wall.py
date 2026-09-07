#!/usr/bin/env python3

"""Analytic hydrostatic regression for a water/VOID state beside SOLID."""

import pathlib
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[7]
sys.path.insert(0, str(source_root / "src/2d-truchas/python"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, mpiexec):
    output_dir = pathlib.Path(tempfile.mkdtemp(prefix="flow_free_surface_hydrostatic_solid_wall_4p_"))
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


def check_solution(data):
    if data.num_steps != 2:
        raise RuntimeError(f"expected initial and final states, got {data.num_steps}")

    for step, expected_time in enumerate((0.0, 1.0)):
        time = data.time(step)
        if abs(time - expected_time) > 1.0e-14:
            raise RuntimeError(f"step {step}: time={time:g}, expected {expected_time:g}")

        centers = data.cell_centers(step)
        water = data.field(step, "vf_water")
        wall = data.field(step, "vf_wall")
        expected_water = ((centers[:, 0] >= 0.25) & (centers[:, 1] < 0.5)).astype(float)
        expected_wall = (centers[:, 0] < 0.25).astype(float)
        if not np.allclose(water, expected_water, rtol=0.0, atol=1.0e-12):
            error = np.max(np.abs(water - expected_water))
            raise RuntimeError(f"step {step}: water-fraction error={error:g}")
        if not np.allclose(wall, expected_wall, rtol=0.0, atol=1.0e-12):
            error = np.max(np.abs(wall - expected_wall))
            raise RuntimeError(f"step {step}: wall-fraction error={error:g}")

        velocity = data.field(step, "velocity")
        fluid = water > 0.0
        velocity_error = np.max(np.abs(velocity[fluid, :2]))
        if velocity_error > 1.0e-11:
            raise RuntimeError(f"step {step}: hydrostatic velocity error={velocity_error:g}")

        pressure = data.field(step, "pressure")
        expected_pressure = 0.5 - centers[:, 1]
        pressure_error = np.max(np.abs(pressure[fluid] - expected_pressure[fluid]))
        if pressure_error > 1.0e-10:
            raise RuntimeError(f"step {step}: hydrostatic pressure error={pressure_error:g}")


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} TRUCHAS_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    try:
        data, output_dir = run_case(executable, input_file, sys.argv[3])
        check_solution(data)
    except (RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        return 1

    print("PASS: hydrostatic water/VOID state beside SOLID is preserved")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
