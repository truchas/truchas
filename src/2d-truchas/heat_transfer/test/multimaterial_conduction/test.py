#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import check_close, execute, finish, run_from_argv


def test():
    run = run_from_argv(Path(__file__).with_name("input.json"))
    flux = 1.0 / (0.5 / 1.0 + 0.5 / 10.0)
    x = run.data.cell_centers(run.final_step)[:, 0]
    expected = np.where(x <= 0.5, 1.0 - flux * x,
                        1.0 - 0.5 * flux - flux / 10.0 * (x - 0.5))
    check_close(run.data.field(run.final_step, "T"), expected, 1.0e-9, "mixed material steady profile")
    finish("multimaterial conduction", run)


if __name__ == "__main__":
    sys.exit(execute(test))
