#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import check_close, execute, finish, run_from_argv


def test():
    run = run_from_argv(Path(__file__).with_name("input.json"))
    low_fraction = run.data.field(run.final_step, "vf_low_capacity")
    high_fraction = run.data.field(run.final_step, "vf_high_capacity")
    check_close(low_fraction + high_fraction, 1.0, 1.0e-12, "volume-fraction sum")
    x = run.data.cell_centers(run.final_step)[:, 0]
    check_close(
        low_fraction,
        np.where(x <= 0.5, 1.0, 0.0),
        1.0e-12,
        "low-capacity volume fraction",
    )
    expected = np.where(x <= 0.5, 2.01, 2.0 + 1.0e-2 / 6.0)
    check_close(run.data.field(run.final_step, "T"), expected, 1.0e-9, "mixed material source")
    finish("multimaterial source", run)


if __name__ == "__main__":
    sys.exit(execute(test))
