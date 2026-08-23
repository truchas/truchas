#!/usr/bin/env python3

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import execute
from test_transient_kelvin import check_kelvin


def test():
    check_kelvin(
        "transient-kelvin-scaled.json",
        "transient Kelvin radiation with scaled constants",
    )


if __name__ == "__main__":
    sys.exit(execute(test))
