#!/usr/bin/env python3

"""Short four-process smoke test for a translating two-fluid interface."""

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
    output_dir = Path(tempfile.mkdtemp(prefix="ns_ht_2d_multifluid_vertical_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"ns_ht_2d returned {result.returncode}")

    data = TruchasVTKHDFData(output_dir / "out.vtkhdf")
    if data.num_steps != 2 or abs(data.time(1) - 0.01) > 1.0e-14:
        raise RuntimeError("incorrect output times")
    inlet0, outlet0 = (data.field(0, name) for name in ("vf_inlet_liquid", "vf_outlet_liquid"))
    inlet1, outlet1 = (data.field(1, name) for name in ("vf_inlet_liquid", "vf_outlet_liquid"))
    for step, inlet, outlet in ((0, inlet0, outlet0), (1, inlet1, outlet1)):
        if np.min(inlet) < -1.0e-12 or np.min(outlet) < -1.0e-12:
            raise RuntimeError(f"step {step}: negative volume fraction")
        if np.max(np.abs(inlet + outlet - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: volume fractions do not sum to one")
    if not np.any((inlet0 > 1.0e-8) & (inlet0 < 1.0 - 1.0e-8)):
        raise RuntimeError("initial interface did not create mixed cells")
    mean_shift = np.mean(inlet1 - inlet0)
    if abs(mean_shift - 0.01) > 2.0e-8:
        raise RuntimeError(f"inlet-material volume change is {mean_shift:g}, expected 0.01")

    for step, inlet, outlet in ((0, inlet0, outlet0), (1, inlet1, outlet1)):
        centers = data.cell_centers(step)
        temp = data.field(step, "T")
        enthalpy = data.field(step, "H")
        temp_error = np.max(np.abs(temp - (1.0 + centers[:, 1])))
        enthalpy_error = np.max(np.abs(enthalpy - (inlet + 4.0*outlet)*(1.0 + centers[:, 1])))
        if max(temp_error, enthalpy_error) > 1.0e-10:
            raise RuntimeError(f"step {step}: temperature={temp_error:g}, enthalpy={enthalpy_error:g}")

    print("PASS: vertical two-fluid interface advances through the coupled step")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
