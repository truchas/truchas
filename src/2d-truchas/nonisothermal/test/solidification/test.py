#!/usr/bin/env python3

"""Four-process regression for phase-aware coupled solidification output."""

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

    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]
    output_dir = Path(tempfile.mkdtemp(prefix="ns_ht_2d_solidification_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--simulation", "ns_ht_2d", "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"ns_ht_2d returned {result.returncode}")

    data = TruchasVTKHDFData(output_dir / "out.vtkhdf")
    expected_times = (0.0, 0.25, 0.5)
    if data.num_steps != len(expected_times):
        raise RuntimeError(f"expected {len(expected_times)} output states, found {data.num_steps}")

    for step, expected_time in enumerate(expected_times):
        if abs(data.time(step) - expected_time) > 1.0e-14:
            raise RuntimeError(f"output {step} has time {data.time(step):g}")
        solid = data.field(step, "vf_phase_change_material_solid")
        liquid = data.field(step, "vf_phase_change_material_liquid")
        temperature = data.field(step, "T")
        enthalpy = data.field(step, "H")
        if not np.isfinite(temperature).all() or not np.isfinite(enthalpy).all():
            raise RuntimeError(f"step {step}: non-finite thermal state")
        if np.min(solid) < -1.0e-12 or np.min(liquid) < -1.0e-12:
            raise RuntimeError(f"step {step}: negative phase volume fraction")
        if np.max(np.abs(solid + liquid - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: phase fractions do not partition the material")
        if step == 0:
            if np.max(np.abs(solid)) > 1.0e-14 or np.max(np.abs(liquid - 1.0)) > 1.0e-14:
                raise RuntimeError("initial material is not entirely liquid")
        elif np.max(solid) <= 1.0e-4 or np.min(liquid) >= 1.0 - 1.0e-4:
            raise RuntimeError(f"step {step}: cooling did not produce a mixed solid/liquid region")

    print("PASS: coupled solidification preserves material volume while producing phase fractions")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
