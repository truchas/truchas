"""Shared execution and VTKHDF checks for the 2D heat-transport tests."""

from dataclasses import dataclass
import pathlib
import re
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[4]
sys.path.insert(0, str(source_root / "python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


@dataclass
class Run:
    output_dir: pathlib.Path
    data: TruchasVTKHDFData
    final_step: int
    final_time: float


def run_case(executable, input_file, nproc=1, mpiexec=None, expected_final_time=1.0e-2):
    """Run ht_2d and return its final VTKHDF result."""

    executable = pathlib.Path(executable).resolve()
    input_file = pathlib.Path(input_file).resolve()
    run_root = pathlib.Path(
        tempfile.mkdtemp(prefix=f"ht_2d_{input_file.stem}_{nproc}p_")
    )
    output_dir = run_root / input_file.stem
    command = [str(executable), str(input_file)]
    if nproc > 1:
        command = [mpiexec, "-n", str(nproc)] + command

    result = subprocess.run(
        command,
        cwd=run_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (run_root / "stdout.txt").write_text(result.stdout)
    if result.returncode != 0:
        raise AssertionError(
            f"ht_2d returned {result.returncode} in {run_root}\n{result.stdout}"
        )

    log_file = output_dir / "run.log"
    if not log_file.exists():
        raise AssertionError(f"ht_2d did not produce run.log in {output_dir}")
    log = log_file.read_text()
    if "unrecoverable integration failure" in log:
        raise AssertionError("unrecoverable integration failure")
    completed = re.findall(r"Completed integration to T =\s*([0-9.E+-]+)", log)
    if not completed:
        raise AssertionError("completion record not found in run.log")
    final_time = float(completed[-1])
    if abs(final_time - expected_final_time) > 1.0e-12:
        raise AssertionError(
            f"final time {final_time:g} differs from {expected_final_time:g}"
        )

    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise AssertionError(f"ht_2d did not produce out.vtkhdf in {output_dir}")
    data = TruchasVTKHDFData(output_file)
    final_step = data.num_steps - 1
    if abs(data.time(final_step) - final_time) > 1.0e-12:
        raise AssertionError("VTKHDF final time disagrees with run.log")
    if not np.all(np.isfinite(data.field(final_step, "T"))):
        raise AssertionError("non-finite temperature in VTKHDF output")

    return Run(output_dir, data, final_step, final_time)


def check_close(actual, expected, tolerance, label):
    error = np.max(np.abs(actual - expected))
    if error > tolerance:
        raise AssertionError(f"{label} error {error:g}")


def run_from_argv(input_file, expected_final_time=1.0e-2):
    if len(sys.argv) not in (2, 4):
        raise AssertionError(f"usage: {sys.argv[0]} HT_2D [NPROC MPIEXEC]")
    nproc = int(sys.argv[2]) if len(sys.argv) == 4 else 1
    mpiexec = sys.argv[3] if len(sys.argv) == 4 else None
    return run_case(
        sys.argv[1], input_file, nproc, mpiexec, expected_final_time
    )


def partition_from_argv(input_file, volume_fraction_names=None):
    if len(sys.argv) != 3:
        raise AssertionError(f"usage: {sys.argv[0]} HT_2D MPIEXEC")
    check_partition_independence(
        sys.argv[1], input_file, sys.argv[2], volume_fraction_names
    )


def check_partition_independence(
    executable, input_file, mpiexec, volume_fraction_names=None
):
    """Compare serial and four-process output in external cell order."""

    serial = run_case(executable, input_file, 1, mpiexec)
    parallel = run_case(executable, input_file, 4, mpiexec)

    serial_ids = serial.data.field(serial.final_step, "ExternalCellIds")
    parallel_ids = parallel.data.field(parallel.final_step, "ExternalCellIds")
    if not np.array_equal(serial_ids, parallel_ids):
        raise AssertionError("external cell IDs differ between serial and 4-process output")

    serial_node_ids = serial.data.field(serial.final_step, "ExternalNodeIds", "point")
    parallel_node_ids = parallel.data.field(parallel.final_step, "ExternalNodeIds", "point")
    if not np.array_equal(serial_node_ids, parallel_node_ids):
        raise AssertionError("external node IDs differ between serial and 4-process output")

    check_close(
        serial.data.cell_centers(serial.final_step),
        parallel.data.cell_centers(parallel.final_step),
        1.0e-14,
        "cell-center",
    )
    check_close(
        serial.data.field(serial.final_step, "T"),
        parallel.data.field(parallel.final_step, "T"),
        1.0e-10,
        "serial/4-process temperature",
    )
    check_close(
        serial.data.field(serial.final_step, "H"),
        parallel.data.field(parallel.final_step, "H"),
        1.0e-10,
        "serial/4-process enthalpy",
    )
    for name in volume_fraction_names or ():
        check_close(
            serial.data.field(serial.final_step, name),
            parallel.data.field(parallel.final_step, name),
            1.0e-14,
            f"serial/4-process {name}",
        )


def finish(label, run=None):
    if run is None:
        print(f"PASS: {label}")
    else:
        print(f"PASS: {label} reached final time {run.final_time:g}")
        print(f"      artifacts: {run.output_dir}")


def execute(test):
    try:
        test()
    except (AssertionError, OSError, ValueError) as error:
        print(f"FAIL: {error}")
        return 1
    return 0
