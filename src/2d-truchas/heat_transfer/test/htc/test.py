#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import check_close, execute, finish, run_from_argv


def test():
    run = run_from_argv(Path(__file__).with_name("input.json"))
    temperature = run.data.field(run.final_step, "T")
    enthalpy = run.data.field(run.final_step, "H")
    if not np.all(temperature < 1.0):
        raise AssertionError("HTC boundary did not cool every cell")
    if not np.all(temperature > 0.0):
        raise AssertionError("HTC boundary produced a non-positive temperature")
    check_close(enthalpy, temperature, 1.0e-10, "unit enthalpy")
    finish("external HTC", run)


if __name__ == "__main__":
    sys.exit(execute(test))
