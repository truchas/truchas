#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, mpiexec):
    run_root = Path(tempfile.mkdtemp(prefix="ns_ht_2d_nat_conv_4p_"))
    command = [mpiexec, "-n", "4", str(executable), str(input_file)]
    result = subprocess.run(
        command, cwd=run_root, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False,
    )
    (run_root / "stdout.txt").write_text(result.stdout)
    if result.returncode != 0:
        raise AssertionError(f"ns_ht_2d returned {result.returncode} in {run_root}\n{result.stdout}")
    output_file = run_root / input_file.stem / "out.vtkhdf"
    return run_root, TruchasVTKHDFData(output_file)


def nearest_field(data, step, x, y, component):
    centers = data.cell_centers(step)
    index = np.argmin((centers[:, 0] - x)**2 + (centers[:, 1] - y)**2)
    return data.field(step, "velocity")[index, component]


def check_reference(data, reference, bdf1):
    if reference.num_steps != 3:
        raise AssertionError(f"expected three reference states, found {reference.num_steps}")
    for step, expected in enumerate((0.0, 1000.0, 30000.0)):
        if abs(reference.time(step) - expected) > 1.0e-10:
            raise AssertionError(
                f"reference output {step} has time {reference.time(step):g}, "
                f"expected {expected:g}"
            )
        if abs(data.time(step) - expected) > 1.0e-10:
            raise AssertionError(
                f"output {step} has time {data.time(step):g}, expected {expected:g}"
            )

    if bdf1:
        tolerances = {
            "T": (3.0e-4, 1.0e-4),
            "velocity": (1.0e-8, 2.0e-8),
        }
    else:
        tolerances = {
            "T": (2.0e-6, 5.0e-5),
            "velocity": (2.0e-10, 2.0e-8),
        }
    for step in (1, 2):
        for name, step_tolerances in tolerances.items():
            error = np.max(np.abs(data.field(step, name) - reference.field(step, name)))
            tolerance = step_tolerances[step - 1]
            if error > tolerance:
                raise AssertionError(
                    f"{name} error at t={data.time(step):g} is {error:g}, "
                    f"expected no more than {tolerance:g}"
                )


def check_case(run_root, data, reference, bdf1):
    if data.num_steps != 3:
        raise AssertionError(f"expected three output states, found {data.num_steps}")
    expected_times = [0.0, 1000.0, 30000.0]
    for step, expected in enumerate(expected_times):
        if abs(data.time(step) - expected) > 1.0e-10:
            raise AssertionError(f"output {step} has time {data.time(step):g}, expected {expected:g}")

    check_reference(data, reference, bdf1)

    initial_centers = data.cell_centers(0)
    initial_temperature = data.field(0, "T")
    if not np.allclose(initial_temperature, 2.0 - initial_centers[:, 0], atol=1.0e-12, rtol=0.0):
        raise AssertionError("initial temperature is not the prescribed linear profile")

    final_step = data.num_steps - 1
    temperature = data.field(final_step, "T")
    pressure = data.field(final_step, "pressure")
    velocity = data.field(final_step, "velocity")
    if not (np.isfinite(temperature).all() and np.isfinite(pressure).all() and np.isfinite(velocity).all()):
        raise AssertionError("natural-convection solution contains non-finite values")
    if temperature.min() < 1.0 - 1.0e-8 or temperature.max() > 2.0 + 1.0e-8:
        raise AssertionError("temperature escaped the boundary-value range")

    horizontal = abs(nearest_field(data, final_step, 0.5, 0.813, 0))
    vertical = abs(nearest_field(data, final_step, 0.178, 0.5, 1))
    expected_horizontal = 7.585e-5
    expected_vertical = 7.685e-5
    tolerance = 0.02
    if abs(horizontal - expected_horizontal) / expected_horizontal > tolerance:
        raise AssertionError(f"horizontal centerline velocity {horizontal:g} differs from benchmark")
    if abs(vertical - expected_vertical) / expected_vertical > tolerance:
        raise AssertionError(f"vertical centerline velocity {vertical:g} differs from benchmark")
    return horizontal, vertical


def test():
    if len(sys.argv) not in (3, 4):
        raise AssertionError(f"usage: {sys.argv[0]} NS_HT_2D MPIEXEC [INPUT_FILE]")
    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[3]).resolve() if len(sys.argv) == 4 \
    else Path(__file__).with_name("bdf2.json").resolve()
    bdf1 = input_file.name == "bdf1.json"
    reference = TruchasVTKHDFData(
        Path(__file__).with_name("reference") / "out.vtkhdf"
    )
    run_root, data = run_case(executable, input_file, sys.argv[2])
    horizontal, vertical = check_case(run_root, data, reference, bdf1)
    scheme = "BDF1" if bdf1 else "BDF2"
    print(f"PASS: natural convection ({scheme}); centerline velocities {horizontal:g}, {vertical:g}; artifacts: {run_root}")


if __name__ == "__main__":
    try:
        test()
    except (AssertionError, OSError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
