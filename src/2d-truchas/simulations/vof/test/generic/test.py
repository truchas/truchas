#!/usr/bin/env python3
"""Regression test for the generic two-dimensional VOF simulation."""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "support"))
from vof_vtkhdf_test_util import run_test


try:
    run_test("generic-vof", 0.002)
except Exception as error:
    print(f"FAIL: {error}")
    raise SystemExit(1)
