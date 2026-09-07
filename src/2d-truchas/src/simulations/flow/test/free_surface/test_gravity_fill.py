#!/usr/bin/env python3

import pathlib
import sys

source_root = pathlib.Path(__file__).resolve().parents[7]
sys.path.insert(0, str(source_root / "src/2d-truchas/python"))
from gravity_pipe_test import run_test


if len(sys.argv) != 4:
    print(f"usage: {sys.argv[0]} TRUCHAS_2D JSON_INPUT MPIEXEC", file=sys.stderr)
    raise SystemExit(2)

run_test(pathlib.Path(sys.argv[1]).resolve(), pathlib.Path(sys.argv[2]).resolve(), sys.argv[3], 0.5, 1.0, 2.5e-3, "fill")
