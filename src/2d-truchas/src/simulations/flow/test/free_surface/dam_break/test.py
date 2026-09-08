#!/usr/bin/env python3

"""Four-process reference regression for the small broken-dam problem."""

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[8]
sys.path.insert(0, str(source_root / "src/2d-truchas/python"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, mpiexec):
    output_dir = Path(tempfile.mkdtemp(prefix="flow_free_surface_dam_break_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--simulation", "flow",
         "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"flow returned {result.returncode} in {output_dir}")
    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"flow did not produce {output_file}")
    return TruchasVTKHDFData(output_file), output_dir


def compare_field(actual, reference, name, step):
    if name in ("velocity", "pressure"):
        if actual.ndim == 1:
            actual_active = np.isfinite(actual)
            reference_active = np.isfinite(reference)
        else:
            actual_active = np.isfinite(actual).all(axis=1)
            reference_active = np.isfinite(reference).all(axis=1)
        if not np.array_equal(actual_active, reference_active):
            raise RuntimeError(f"step {step}: {name} active-cell mask differs")
        actual = actual[actual_active]
        reference = reference[reference_active]
    elif not np.isfinite(actual).all():
        raise RuntimeError(f"step {step}: non-finite {name}")

    denominator = max(np.linalg.norm(reference), 1.0e-14)
    error = np.linalg.norm(actual - reference) / denominator
    if error > 2.0e-6:
        raise RuntimeError(f"step {step}: {name} relative L2 error={error:g}")


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
    expected_times = np.array([0.0, 0.1, 0.5, 1.0])
    if data.num_steps != len(expected_times) or reference.num_steps != len(expected_times):
        raise RuntimeError("expected initial output plus three scheduled outputs")

    for step, expected_time in enumerate(expected_times):
        if abs(data.time(step) - expected_time) > 1.0e-12:
            raise RuntimeError(f"output {step} has time {data.time(step):g}")
        if abs(reference.time(step) - expected_time) > 1.0e-12:
            raise RuntimeError(f"reference output {step} has time {reference.time(step):g}")

        for name in ("vf_water", "vf_VOID", "velocity", "pressure"):
            compare_field(data.field(step, name), reference.field(step, name), name, step)

        water = data.field(step, "vf_water")
        void = data.field(step, "vf_VOID")
        if np.min(water) < -1.0e-12 or np.min(void) < -1.0e-12:
            raise RuntimeError(f"step {step}: negative volume fraction")
        if np.max(np.abs(water + void - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: volume fractions do not sum to one")

    print("PASS: small broken-dam flow agrees with four-process reference")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
