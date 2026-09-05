#!/usr/bin/env python3

"""Reference regression for a two-fluid rising bubble."""

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[6]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, mpiexec):
    output_dir = Path(tempfile.mkdtemp(prefix="ns_2d_bubble_rise_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--simulation", "flow", "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False)
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"ns_2d returned {result.returncode} in {output_dir}")
    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"ns_2d did not produce {output_file}")
    return TruchasVTKHDFData(output_file), output_dir


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} TRUCHAS_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[2]).resolve()
    mpiexec = Path(sys.argv[3])
    reference_file = Path(__file__).with_name("reference") / "out.vtkhdf"
    data, output_dir = run_case(executable, input_file, mpiexec)
    reference = TruchasVTKHDFData(reference_file)

    expected_times = np.arange(0.0, 8.0 + 1.0, 1.0)
    if data.num_steps != len(expected_times) or reference.num_steps != len(expected_times):
        raise RuntimeError("expected initial output plus eight scheduled outputs")

    initial_vf = data.field(0, "vf_light_liquid")
    final_vf = data.field(data.num_steps - 1, "vf_light_liquid")
    if np.max(np.abs(final_vf - initial_vf)) < 1.0e-3:
        raise RuntimeError("bubble did not move or deform")

    for step, expected_time in enumerate(expected_times):
        if abs(data.time(step) - expected_time) > 1.0e-12:
            raise RuntimeError(f"output {step} has time {data.time(step):g}")
        if abs(reference.time(step) - expected_time) > 1.0e-12:
            raise RuntimeError(f"reference output {step} has time {reference.time(step):g}")

        light = data.field(step, "vf_light_liquid")
        reference_light = reference.field(step, "vf_light_liquid")
        if not np.isfinite(light).all():
            raise RuntimeError(f"step {step}: non-finite volume fraction")
        relative_l1_error = np.sum(np.abs(light - reference_light)) / np.sum(np.abs(reference_light))
        if relative_l1_error > 2.0e-6:
            raise RuntimeError(f"step {step}: relative light-fluid VF L1 error={relative_l1_error:g}")

        heavy = data.field(step, "vf_heavy_liquid")
        if np.min(light) < -1.0e-12 or np.min(heavy) < -1.0e-12:
            raise RuntimeError(f"step {step}: negative volume fraction")
        if np.max(np.abs(light + heavy - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: volume fractions do not sum to one")

        velocity = data.field(step, "velocity")
        pressure = data.field(step, "pressure")
        if not np.isfinite(velocity).all() or not np.isfinite(pressure).all():
            raise RuntimeError(f"step {step}: non-finite flow solution")

    if np.max(np.abs(data.field(data.num_steps - 1, "velocity")[:, :2])) < 1.0e-2:
        raise RuntimeError("bubble motion did not generate a flow field")

    print("PASS: rising-bubble volume fraction agrees with four-process reference")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
