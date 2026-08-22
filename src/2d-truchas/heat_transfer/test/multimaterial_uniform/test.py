#!/usr/bin/env python3

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import execute, finish, run_from_argv


def test():
    run = run_from_argv(Path(__file__).with_name("input.json"))
    finish("multimaterial uniform", run)


if __name__ == "__main__":
    sys.exit(execute(test))
