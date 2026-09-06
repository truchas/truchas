#!/usr/bin/env python3

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "support"))
from vof_vtkhdf_test_util import run_test


if __name__ == "__main__":
    try:
        run_test("axial-flow", 0.1, field_tolerance=1.0e-6)
    except RuntimeError as error:
        print(f"FAIL: {error}")
        sys.exit(1)
