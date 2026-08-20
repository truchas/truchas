#!/usr/bin/env python3

"""End-to-end MPI-consistency test for the JSON-driven 2D NS application."""

import pathlib
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[4]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, nproc, mpiexec):
    output_dir = pathlib.Path(tempfile.mkdtemp(prefix=f"ns_2d_lid_{nproc}p_"))
    command = [str(executable), str(input_file)]
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


def final_solution(data, case):
    if data.num_steps != 3:
        raise RuntimeError(f"{case}: wrote {data.num_steps} states, expected initial plus two scheduled outputs")
    if abs(data.time(1) - 0.053) > 1.0e-12:
        raise RuntimeError(f"{case}: intermediate output time is {data.time(1):g}, expected 0.053")
    step = data.num_steps - 1
    if abs(data.time(step) - 0.1) > 1.0e-12:
        raise RuntimeError(f"{case}: final time is {data.time(step):g}, expected 0.1")
    centers = data.cell_centers(step)
    pressure = data.field(step, "pressure")
    velocity = data.field(step, "velocity")
    if not np.isfinite(pressure).all() or not np.isfinite(velocity).all():
        raise RuntimeError(f"{case}: non-finite solution value")
    if np.max(np.abs(velocity[:, :2])) < 1.0e-2:
        raise RuntimeError(f"{case}: lid motion did not reach the fluid")
    if np.max(np.abs(velocity[:, 1])) < 1.0e-4:
        raise RuntimeError(f"{case}: cavity flow lacks transverse motion")
    return centers, pressure, velocity


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} NS_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]

    try:
        serial, serial_dir = run_case(executable, input_file, 1, mpiexec)
        parallel, parallel_dir = run_case(executable, input_file, 4, mpiexec)
        serial_centers, serial_pressure, serial_velocity = final_solution(serial, "serial")
        parallel_centers, parallel_pressure, parallel_velocity = final_solution(parallel, "parallel")
    except RuntimeError as error:
        print(f"FAIL: {error}")
        return 1

    if not np.allclose(serial_centers, parallel_centers, rtol=0.0, atol=1.0e-14):
        print("FAIL: serial and parallel cell centers differ")
        return 1
    pressure_error = np.max(np.abs(serial_pressure - parallel_pressure))
    velocity_error = np.max(np.abs(serial_velocity - parallel_velocity))
    if pressure_error > 1.0e-10 or velocity_error > 1.0e-10:
        print(f"FAIL: serial/parallel errors are pressure={pressure_error:g}, velocity={velocity_error:g}")
        return 1

    print("PASS: integrated Navier--Stokes lid cavity is finite and MPI independent")
    print(f"      serial artifacts: {serial_dir}")
    print(f"      parallel artifacts: {parallel_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
