#!/usr/bin/env python3

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import check_close, execute, finish, run_from_argv
from TruchasVTKHDFData import TruchasVTKHDFData


def check_kelvin(input_name, label):
    run = run_from_argv(
        Path(__file__).with_name(input_name),
        expected_final_time=1.0e-3,
    )
    reference = TruchasVTKHDFData(
        Path(__file__).with_name("transient-kelvin") / "out.vtkhdf"
    )
    if reference.num_steps != 3 or run.data.num_steps != 3:
        raise AssertionError("expected initial and two transient output states")
    for step in (1, 2):
        if abs(run.data.time(step) - reference.time(step)) > 1.0e-14:
            raise AssertionError(f"reference time mismatch at step {step}")
        for name in ("T", "H"):
            check_close(
                run.data.field(step, name),
                reference.field(step, name),
                1.0e-4,
                f"Kelvin reference {name} at step {step}",
            )
    finish(label, run)


def test():
    check_kelvin("transient-kelvin.json", "transient Kelvin radiation")


if __name__ == "__main__":
    sys.exit(execute(test))
