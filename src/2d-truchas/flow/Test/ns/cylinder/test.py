#!/usr/bin/env python3

"""Reference regression for laminar flow past a solid cylinder."""

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[6]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def run_case(executable, input_file, mpiexec):
    output_dir = Path(tempfile.mkdtemp(prefix="ns_2d_cylinder_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False)
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"ns_2d returned {result.returncode} in {output_dir}")
    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"ns_2d did not produce {output_file}")
    return TruchasVTKHDFData(output_file), output_dir


def relative_l2_error(actual, reference):
    return np.linalg.norm(actual - reference) / np.linalg.norm(reference)


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} NS_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[2]).resolve()
    mpiexec = Path(sys.argv[3])
    reference_file = Path(__file__).with_name("reference") / "out.vtkhdf"
    data, output_dir = run_case(executable, input_file, mpiexec)
    reference = TruchasVTKHDFData(reference_file)

    expected_times = np.array([0.0, 0.25, 0.5, 1.0, 2.0])
    if data.num_steps != len(expected_times) or reference.num_steps != len(expected_times):
        raise RuntimeError("expected initial output plus four scheduled outputs")

    for step, expected_time in enumerate(expected_times):
        if abs(data.time(step) - expected_time) > 1.0e-12:
            raise RuntimeError(f"output {step} has time {data.time(step):g}")
        if abs(reference.time(step) - expected_time) > 1.0e-12:
            raise RuntimeError(f"reference output {step} has time {reference.time(step):g}")

        for name in ("vf_cylinder", "vf_fluid", "velocity", "pressure"):
            actual = data.field(step, name)
            expected = reference.field(step, name)
            error = relative_l2_error(actual, expected)
            if error > 2.0e-6:
                raise RuntimeError(f"step {step}: {name} relative L2 error={error:g}")
            if not np.isfinite(actual).all():
                raise RuntimeError(f"step {step}: non-finite {name}")

        cylinder = data.field(step, "vf_cylinder")
        fluid = data.field(step, "vf_fluid")
        if np.min(cylinder) < -1.0e-12 or np.min(fluid) < -1.0e-12:
            raise RuntimeError(f"step {step}: negative material volume fraction")
        if np.max(np.abs(cylinder + fluid - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: material fractions do not sum to one")

    velocity = data.field(data.num_steps - 1, "velocity")[:, :2]
    if np.min(velocity[:, 0]) > -1.0e-3 or np.max(np.abs(velocity[:, 1])) < 1.0e-2:
        raise RuntimeError("final flow field does not show a nontrivial cylinder wake")

    print("PASS: cylinder flow agrees with four-process reference")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
