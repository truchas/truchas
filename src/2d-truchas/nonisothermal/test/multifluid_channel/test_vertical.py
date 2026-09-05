#!/usr/bin/env python3

"""Analytic four-process regression for a translating two-fluid interface."""

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} TRUCHAS_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable, input_file, mpiexec = map(Path, sys.argv[1:])
    output_dir = Path(tempfile.mkdtemp(prefix="ns_ht_2d_multifluid_vertical_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--simulation", "flow_thermal", "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"ns_ht_2d returned {result.returncode}")

    data = TruchasVTKHDFData(output_dir / "out.vtkhdf")
    expected_times = (0.0, 0.0025, 0.005, 0.0075, 0.01)
    if data.num_steps != len(expected_times):
        raise RuntimeError("incorrect output times")
    dx = 1.0 / 16.0
    volume_fraction_tolerance = 5.0e-5
    for step, expected_time in enumerate(expected_times):
        if abs(data.time(step) - expected_time) > 1.0e-14:
            raise RuntimeError(f"output {step} has time {data.time(step):g}")
        inlet = data.field(step, "vf_inlet_liquid")
        outlet = data.field(step, "vf_outlet_liquid")
        centers = data.cell_centers(step)
        interface = 0.47 + expected_time
        expected_inlet = np.clip(
            (interface - (centers[:, 0] - 0.5*dx)) / dx, 0.0, 1.0
        )
        if np.min(inlet) < -1.0e-12 or np.min(outlet) < -1.0e-12:
            raise RuntimeError(f"step {step}: negative volume fraction")
        if np.max(np.abs(inlet + outlet - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: volume fractions do not sum to one")
        distribution_error = np.max(np.abs(inlet - expected_inlet))
        temp = data.field(step, "T")
        enthalpy = data.field(step, "H")
        temp_error = np.max(np.abs(temp - (1.0 + centers[:, 1])))
        enthalpy_error = np.max(
            np.abs(enthalpy - (inlet + 4.0*outlet)*(1.0 + centers[:, 1]))
        )
        velocity = data.field(step, "velocity")
        velocity_error = np.max(np.abs(velocity[:, :2] - [1.0, 0.0]))
        # The geometric initializer and tracker use a finite refinement level,
        # so the axis-aligned area fraction is only resolved to this scale.
        if (distribution_error > volume_fraction_tolerance or
                max(temp_error, enthalpy_error, velocity_error) > 2.0e-9):
            raise RuntimeError(
                f"step {step}: distribution={distribution_error:g}, "
                f"temperature={temp_error:g}, enthalpy={enthalpy_error:g}, "
                f"velocity={velocity_error:g}"
            )
        if step == 0 and not np.any((inlet > 1.0e-8) & (inlet < 1.0 - 1.0e-8)):
            raise RuntimeError("initial interface did not create mixed cells")

    print("PASS: vertical two-fluid interface follows the analytic translation")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
