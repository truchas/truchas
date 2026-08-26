#!/usr/bin/env python3

"""Four-process regression for two-fluid inviscid plug flow."""

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
        print(f"usage: {sys.argv[0]} NS_HT_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable, input_file, mpiexec = map(Path, sys.argv[1:])
    output_dir = Path(tempfile.mkdtemp(prefix="ns_ht_2d_multifluid_channel_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"ns_ht_2d returned {result.returncode}")

    data = TruchasVTKHDFData(output_dir / "out.vtkhdf")
    if data.num_steps != 2 or abs(data.time(1) - 0.01) > 1.0e-14:
        raise RuntimeError("incorrect output times")
    lower0, upper0 = (data.field(0, name) for name in ("vf_lower_liquid", "vf_upper_liquid"))
    lower1, upper1 = (data.field(1, name) for name in ("vf_lower_liquid", "vf_upper_liquid"))
    if np.max(np.abs(lower0 + upper0 - 1.0)) > 1.0e-12:
        raise RuntimeError("initial material fractions do not sum to one")
    if not np.any((lower0 > 1.0e-8) & (lower0 < 1.0 - 1.0e-8)):
        raise RuntimeError("interface did not create mixed cells")
    vf_error = max(np.max(np.abs(lower1 - lower0)), np.max(np.abs(upper1 - upper0)))
    if vf_error > 1.0e-10:
        raise RuntimeError(f"plug flow changed the horizontal material interface: {vf_error:g}")

    for step, lower, upper in ((0, lower0, upper0), (1, lower1, upper1)):
        centers = data.cell_centers(step)
        temp = data.field(step, "T")
        enthalpy = data.field(step, "H")
        temp_error = np.max(np.abs(temp - (1.0 + centers[:, 1])))
        enthalpy_error = np.max(np.abs(enthalpy - (lower + 4.0*upper)*(1.0 + centers[:, 1])))
        if max(temp_error, enthalpy_error) > 2.0e-9:
            raise RuntimeError(f"step {step}: temperature={temp_error:g}, enthalpy={enthalpy_error:g}")

    print("PASS: two-fluid split-cell plug flow preserves composition and thermal state")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
