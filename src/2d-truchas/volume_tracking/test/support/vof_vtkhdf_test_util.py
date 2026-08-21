"""Shared VTKHDF helpers for the two-dimensional VOF integration tests."""

import pathlib
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_test(label, final_time):
    if len(sys.argv) != 5:
        raise RuntimeError(
            f"usage: {sys.argv[0]} VOF_EXECUTABLE INPUT_FILE REFERENCE_FILE MPIEXEC"
        )

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    reference_file = pathlib.Path(sys.argv[3]).resolve()
    mpiexec = sys.argv[4]
    output_dir = pathlib.Path(tempfile.mkdtemp(prefix=f"vof_{label}_4p_"))
    command = [mpiexec, "-n", "4", str(executable), str(input_file), str(reference_file)]
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
        raise RuntimeError(f"{label}: VOF executable returned {result.returncode}")

    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"{label}: VOF executable did not produce {output_file}")
    data = TruchasVTKHDFData(output_file)
    if data.num_steps != 2:
        raise RuntimeError(f"{label}: expected two VTKHDF steps, found {data.num_steps}")
    if abs(data.time(0)) > 1.0e-12:
        raise RuntimeError(f"{label}: initial time is {data.time(0):g}, expected 0")
    if abs(data.time(1) - final_time) > 1.0e-12:
        raise RuntimeError(
            f"{label}: final time is {data.time(1):g}, expected {final_time:g}"
        )

    reference = _read_reference(reference_file)
    expected_ids = np.arange(1, reference.shape[0] + 1)
    if not np.array_equal(reference[:, 0].astype(int), expected_ids):
        raise RuntimeError(f"{label}: text reference is not in external cell-ID order")

    vfrac_1 = data.field(1, "volume-fraction-1")
    vfrac_2 = data.field(1, "volume-fraction-2")
    if vfrac_1.shape != reference[:, 1].shape:
        raise RuntimeError(
            f"{label}: VTKHDF has {vfrac_1.size} cells, reference has {reference.shape[0]}"
        )
    error = np.max(np.abs(vfrac_1 - reference[:, 1]))
    if error > 1.0e-7:
        raise RuntimeError(f"{label}: VTKHDF/text volume-fraction error is {error:g}")
    sum_error = np.max(np.abs(vfrac_1 + vfrac_2 - 1.0))
    if sum_error > 1.0e-10:
        raise RuntimeError(f"{label}: volume fractions do not sum to one ({sum_error:g})")

    print(f"PASS: {label} VTKHDF output matches its text reference")
    print(f"      artifact: {output_dir}")


def _read_reference(filename):
    """Read the legacy Fortran reference format, including exponentless values."""
    values = []
    with open(filename, encoding="ascii") as reference_file:
        next(reference_file)
        for line in reference_file:
            fields = line.split()
            if len(fields) != 2:
                raise RuntimeError(f"invalid text reference line: {line.rstrip()}")
            value = fields[1]
            if "E" not in value.upper() and "D" not in value.upper():
                for index, character in enumerate(value[1:], start=1):
                    if character in "+-":
                        value = value[:index] + "E" + value[index:]
                        break
            values.append((int(fields[0]), float(value.replace("D", "E"))))
    return np.asarray(values)
