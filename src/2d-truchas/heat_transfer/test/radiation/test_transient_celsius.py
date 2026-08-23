#!/usr/bin/env python3

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import check_close, execute, finish, run_from_argv
from TruchasVTKHDFData import TruchasVTKHDFData


def test():
    run = run_from_argv(
        Path(__file__).with_name("transient-celsius.json"),
        expected_final_time=1.0e-3,
    )
    reference = TruchasVTKHDFData(
        Path(__file__).with_name("transient-kelvin") / "out.vtkhdf"
    )
    if reference.num_steps != 3:
        raise AssertionError("Kelvin reference is missing a transient output")
    if abs(run.data.time(run.final_step) - reference.time(2)) > 1.0e-14:
        raise AssertionError("Celsius/Kelvin final time mismatch")
    for name in ("T", "H"):
        check_close(
            run.data.field(run.final_step, name) + 273.15,
            reference.field(2, name),
            1.0e-4,
            f"Celsius/Kelvin shifted {name}",
        )
    finish("transient Celsius radiation", run)


if __name__ == "__main__":
    sys.exit(execute(test))
