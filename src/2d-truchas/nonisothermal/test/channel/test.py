#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(source_root / "python" / "truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, nproc=1, mpiexec=None):
    run_root = Path(tempfile.mkdtemp(prefix=f"ns_ht_2d_channel_{nproc}p_"))
    command = [str(executable), str(input_file)]
    if nproc > 1:
        command = [mpiexec, "-n", str(nproc)] + command
    result = subprocess.run(
        command, cwd=run_root, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False,
    )
    (run_root / "stdout.txt").write_text(result.stdout)
    if result.returncode != 0:
        raise AssertionError(f"ns_ht_2d returned {result.returncode} in {run_root}\n{result.stdout}")
    return run_root, TruchasVTKHDFData(run_root / "input" / "out.vtkhdf")


def check_channel(run_root, data):
    log = (run_root / "input" / "run.log").read_text()
    for text in ("Beginning integration", "Completed integration", "Timing Summary:"):
        if text not in log:
            raise AssertionError(f"missing {text!r} in simulation log")
    if data.num_steps != 2:
        raise AssertionError(f"expected two output steps, found {data.num_steps}")
    if abs(data.time(0)) > 1.0e-14 or abs(data.time(1) - 1.0e-2) > 1.0e-12:
        raise AssertionError("incorrect VTKHDF output times")
    temp0 = data.field(0, "T")
    temp1 = data.field(1, "T")
    if not np.allclose(temp0, 2.0, rtol=0.0, atol=1.0e-12):
        raise AssertionError("incorrect initial temperature")
    if not np.all(np.isfinite(temp1)) or temp1.max() <= 2.01 or temp1.min() < 2.0:
        raise AssertionError("warm inflow did not produce the expected bounded thermal response")
    for name in ("H", "pressure", "velocity"):
        if not np.all(np.isfinite(data.field(1, name))):
            raise AssertionError(f"non-finite {name} output")


def test():
    executable = Path(sys.argv[1]).resolve()
    input_file = Path(__file__).with_name("input.json").resolve()
    run_root, serial = run_case(executable, input_file)
    check_channel(run_root, serial)
    if len(sys.argv) == 3:
        parallel_root, parallel = run_case(executable, input_file, 4, sys.argv[2])
        check_channel(parallel_root, parallel)
        for name in ("T", "H", "pressure", "velocity"):
            if not np.allclose(serial.field(1, name), parallel.field(1, name), rtol=0.0, atol=1.0e-9):
                raise AssertionError(f"serial/parallel {name} differs")
        print(f"PASS: channel serial/parallel; artifacts: {run_root}, {parallel_root}")
    else:
        print(f"PASS: channel warm-inflow; artifacts: {run_root}")


if __name__ == "__main__":
    try:
        test()
    except (AssertionError, OSError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
