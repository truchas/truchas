#!/usr/bin/env python3

"""Four-process regression for buoyancy-driven solidification."""

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

    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[2]).resolve()
    mpiexec = sys.argv[3]
    run_root = Path(tempfile.mkdtemp(prefix="ns_ht_2d_solidification_flow_4p_"))
    result = subprocess.run(
        [mpiexec, "-n", "4", str(executable), "--output-dir", str(run_root),
         "--force", str(input_file)],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    (run_root / "stdout.txt").write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"ns_ht_2d returned {result.returncode} in {run_root}\n{result.stdout}")

    data = TruchasVTKHDFData(run_root / "out.vtkhdf")
    expected_times = (0.0, 1.0, 5.0, 20.0, 100.0, 200.0)
    if data.num_steps != len(expected_times):
        raise RuntimeError(f"expected {len(expected_times)} output states, found {data.num_steps}")

    final_speed = 0.0
    final_solid = 0.0
    final_mixed = 0
    final_pure_solid = 0
    for step, expected_time in enumerate(expected_times):
        if abs(data.time(step) - expected_time) > 1.0e-12:
            raise RuntimeError(f"output {step} has time {data.time(step):g}")

        solid = data.field(step, "vf_phase_change_material_solid")
        liquid = data.field(step, "vf_phase_change_material_liquid")
        temperature = data.field(step, "T")
        enthalpy = data.field(step, "H")
        pressure = data.field(step, "pressure")
        velocity = data.field(step, "velocity")
        if np.min(solid) < -1.0e-12 or np.min(liquid) < -1.0e-12:
            raise RuntimeError(f"step {step}: negative phase volume fraction")
        if np.max(np.abs(solid + liquid - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: phase fractions do not partition the material")

        active = liquid > 1.0e-2
        if not all(np.isfinite(field[active]).all() for field in
                   (temperature, enthalpy, pressure, velocity)):
            raise RuntimeError(f"step {step}: non-finite active-cell solution value")

        if step == 0:
            if np.max(np.abs(solid)) > 1.0e-14 or np.max(np.abs(liquid - 1.0)) > 1.0e-14:
                raise RuntimeError("initial material is not entirely liquid")
        if step == len(expected_times) - 1:
            final_solid = np.max(solid)
            final_mixed = np.count_nonzero((solid > 1.0e-6) & (liquid > 1.0e-6))
            final_pure_solid = np.count_nonzero(liquid <= 1.0e-2)
            final_speed = np.max(np.linalg.norm(velocity[active], axis=1))

    if final_solid <= 1.0e-2 or final_mixed == 0 or final_pure_solid == 0:
        raise RuntimeError("flow case did not produce solid/mixed cells")
    if final_speed <= 1.0e-3:
        raise RuntimeError("flow case did not produce buoyancy-driven motion")

    print("PASS: buoyancy-driven flow couples to solidification")
    print(f"      final maximum speed: {final_speed:g}")
    print(f"      final maximum solid fraction: {final_solid:g}")
    print(f"      final pure solid cells: {final_pure_solid}")
    print(f"      final mixed cells: {final_mixed}")
    print(f"      artifacts: {run_root}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
