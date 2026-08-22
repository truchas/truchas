"""Shared VTKHDF helpers for the two-dimensional VOF integration tests."""

import pathlib
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_test(label, final_time, field_tolerance=1.0e-7):
    if len(sys.argv) != 5:
        raise RuntimeError(
            f"usage: {sys.argv[0]} VOF_EXECUTABLE INPUT_FILE GOLD_FILE MPIEXEC"
        )

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    gold_file = pathlib.Path(sys.argv[3]).resolve()
    mpiexec = sys.argv[4]
    output_dir = pathlib.Path(tempfile.mkdtemp(prefix=f"vof_{label}_4p_"))
    result = subprocess.run(
        [mpiexec, "-n", "4", str(executable), str(input_file)],
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
    gold = TruchasVTKHDFData(gold_file)
    if data.num_steps != 2 or gold.num_steps != 2:
        raise RuntimeError(
            f"{label}: expected two VTKHDF steps, found {data.num_steps} and {gold.num_steps}"
        )
    if abs(gold.time(0)) > 1.0e-12 or abs(data.time(0) - gold.time(0)) > 1.0e-12:
        raise RuntimeError(f"{label}: initial time does not match gold output")
    if abs(gold.time(1) - final_time) > 1.0e-12:
        raise RuntimeError(
            f"{label}: gold final time is {gold.time(1):g}, expected {final_time:g}"
        )
    if abs(data.time(1) - gold.time(1)) > 1.0e-12:
        raise RuntimeError(f"{label}: final time does not match gold output")

    for step in range(2):
        fields = [data.field(step, f"volume-fraction-{material}") for material in (1, 2)]
        gold_fields = [gold.field(step, f"volume-fraction-{material}") for material in (1, 2)]
        for material, (values, expected) in enumerate(zip(fields, gold_fields), start=1):
            if values.shape != expected.shape:
                raise RuntimeError(
                    f"{label}: material {material} has {values.size} cells, "
                    f"gold has {expected.size}"
                )
            error = np.max(np.abs(values - expected))
            if error > field_tolerance:
                raise RuntimeError(
                    f"{label}: material {material}, step {step} error against gold is {error:g}"
                )
        sum_error = np.max(np.abs(fields[0] + fields[1] - 1.0))
        if sum_error > 1.0e-10:
            raise RuntimeError(
                f"{label}: step {step} volume fractions do not sum to one ({sum_error:g})"
            )

    print(f"PASS: {label} VTKHDF output matches serial gold data")
    print(f"      artifact: {output_dir}")
