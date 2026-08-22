#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import check_close, execute, finish, run_from_argv


def test():
    run = run_from_argv(Path(__file__).with_name("input.json"))
    expected = np.sqrt(1.0 + 2.0e-2) - 1.0
    check_close(run.data.field(run.final_step, "T"), expected, 1.0e-8, "nonlinear temperature")
    finish("nonlinear", run)


if __name__ == "__main__":
    sys.exit(execute(test))
