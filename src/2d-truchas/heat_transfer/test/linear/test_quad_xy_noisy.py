#!/usr/bin/env python3

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import check_close, execute, finish, run_from_argv


def test():
    run = run_from_argv(Path(__file__).with_name("input_quad_xy_noisy.json"))
    centers = run.data.cell_centers(run.final_step)
    expected = 1.0 + centers[:, 0] + 2.0 * centers[:, 1]
    check_close(run.data.field(run.final_step, "T"), expected, 1.0e-10, "noisy quad temperature")
    check_close(run.data.field(run.final_step, "H"), 2.0 * expected, 1.0e-10, "noisy quad enthalpy")
    finish("noisy quad xy linear", run)


if __name__ == "__main__":
    sys.exit(execute(test))
