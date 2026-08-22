#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import check_close, execute, finish, run_from_argv


def test():
    run = run_from_argv(Path(__file__).with_name("input.json"))
    x = run.data.cell_centers(run.final_step)[:, 0]
    low_fraction = np.clip((0.4375 - (x - 0.0625)) / 0.125, 0.0, 1.0)
    heat_capacity = low_fraction + 3.0 * (1.0 - low_fraction)
    expected = 2.0 + 1.0e-2 / heat_capacity
    check_close(run.data.field(run.final_step, "T"), expected, 1.0e-9, "mixed-cell source")
    finish("multimaterial mixed source", run)


if __name__ == "__main__":
    sys.exit(execute(test))
