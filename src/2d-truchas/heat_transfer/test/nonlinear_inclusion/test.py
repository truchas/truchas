#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import (
    TruchasVTKHDFData,
    check_close,
    execute,
    finish,
    run_from_argv,
)


def check_reference(actual, reference, tolerance, label):
    if actual.shape != reference.shape:
        raise AssertionError(
            f"{label} shapes differ: {actual.shape} != {reference.shape}"
        )
    error = np.max(np.abs((actual - reference) / reference))
    if not np.isfinite(error) or error > tolerance:
        raise AssertionError(f"{label} relative error {error:g}")


def test():
    run = run_from_argv(
        Path(__file__).with_name("input.json"), expected_final_time=0.1
    )
    reference = TruchasVTKHDFData(
        Path(__file__).with_name("reference") / "out.vtkhdf"
    )
    reference_step = reference.num_steps - 1
    if abs(reference.time(reference_step) - run.final_time) > 1.0e-12:
        raise AssertionError("reference and candidate final times differ")

    check_close(
        run.data.cell_centers(run.final_step),
        reference.cell_centers(reference_step),
        1.0e-14,
        "mesh cell centers",
    )
    for name in ("T", "H"):
        check_reference(
            run.data.field(run.final_step, name),
            reference.field(reference_step, name),
            1.0e-5,
            name,
        )
    finish("nonlinear material inclusion", run)


if __name__ == "__main__":
    sys.exit(execute(test))
