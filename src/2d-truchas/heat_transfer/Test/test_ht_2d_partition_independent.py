#!/usr/bin/env python3

"""Verify VTKHDF cell fields are independent of MPI mesh partitioning."""

import pathlib
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[4]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, nproc, mpiexec):
    output_dir = pathlib.Path(tempfile.mkdtemp(prefix=f"ht_2d_partition_{nproc}p_"))
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
        raise RuntimeError(f"ht_2d returned {result.returncode} with {nproc} processes")
    return TruchasVTKHDFData(output_dir / "out.vtkhdf")


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} HT_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]

    try:
        serial = run_case(executable, input_file, 1, mpiexec)
        parallel = run_case(executable, input_file, 4, mpiexec)
    except RuntimeError as error:
        print(f"FAIL: {error}")
        return 1

    serial_step = serial.num_steps - 1
    parallel_step = parallel.num_steps - 1
    serial_ids = serial.field(serial_step, "ExternalCellIds")
    parallel_ids = parallel.field(parallel_step, "ExternalCellIds")
    if not np.array_equal(serial_ids, parallel_ids):
        print("FAIL: external cell IDs differ between serial and 4-process output")
        return 1

    serial_node_ids = serial.field(serial_step, "ExternalNodeIds", "point")
    parallel_node_ids = parallel.field(parallel_step, "ExternalNodeIds", "point")
    if not np.array_equal(serial_node_ids, parallel_node_ids):
        print("FAIL: external node IDs differ between serial and 4-process output")
        return 1

    serial_centers = serial.cell_centers(serial_step)
    parallel_centers = parallel.cell_centers(parallel_step)
    if not np.allclose(serial_centers, parallel_centers, rtol=0.0, atol=1.0e-14):
        print("FAIL: cell centers differ between serial and 4-process output")
        return 1

    serial_temperature = serial.field(serial_step, "T")
    parallel_temperature = parallel.field(parallel_step, "T")
    error = np.max(np.abs(serial_temperature - parallel_temperature))
    if error > 1.0e-10:
        print(f"FAIL: serial and 4-process temperatures differ by {error:g}")
        return 1

    print("PASS: serial and 4-process cell output agree in external cell order")
    return 0


if __name__ == "__main__":
    sys.exit(main())
