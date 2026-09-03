#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import TruchasVTKHDFData, execute, finish, run_from_argv


def check_reference(actual, reference, tolerance, label):
    if actual.shape != reference.shape:
        raise AssertionError(
            f"{label} shapes differ: {actual.shape} != {reference.shape}"
        )
    error = np.max(np.abs(actual - reference))
    scale = max(1.0, np.max(np.abs(reference)))
    if not np.isfinite(error) or error > tolerance * scale:
        raise AssertionError(f"{label} absolute error {error:g}")


def test():
    run = run_from_argv(
        Path(__file__).with_name("input.json"), expected_final_time=1000.0
    )
    reference = TruchasVTKHDFData(
        Path(__file__).with_name("reference") / "out.vtkhdf"
    )
    if reference.num_steps != run.data.num_steps:
        raise AssertionError(
            f"reference and candidate have different numbers of outputs: "
            f"{reference.num_steps} != {run.data.num_steps}"
        )

    fields = (
        "T",
        "H",
        "vf_phase_change_material_solid",
        "vf_phase_change_material_liquid",
    )
    for step in range(run.data.num_steps):
        if abs(run.data.time(step) - reference.time(step)) > 1.0e-12:
            raise AssertionError(f"reference and candidate times differ at output {step}")
        for name in fields:
            check_reference(
                run.data.field(step, name),
                reference.field(step, name),
                1.0e-8,
                f"{name} at t={run.data.time(step):g}",
            )
    finish("multiphase solidification reference comparison", run)


if __name__ == "__main__":
    sys.exit(execute(test))
