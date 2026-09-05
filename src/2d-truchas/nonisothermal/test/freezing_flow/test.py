#!/usr/bin/env python3

"""Four-process reference regression for inviscid freezing plug flow."""

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

source_root = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


EXPECTED_TIMES = (0.0, 0.5, 3.0)
PHASE_FIELDS = (
    "vf_freezing_material_liquid",
    "vf_freezing_material_solid",
)


def run_case(executable, input_file, mpiexec):
    output_dir = Path(tempfile.mkdtemp(prefix="ns_ht_2d_freezing_flow_4p_"))
    result = subprocess.run(
        [str(mpiexec), "-n", "4", str(executable), "--simulation", "flow_thermal", "--output-dir", ".", "--force",
         str(input_file)],
        cwd=output_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (output_dir / "stdout.txt").write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(
            f"ns_ht_2d returned {result.returncode} in {output_dir}\n"
            f"{result.stdout}"
        )
    output_file = output_dir / "out.vtkhdf"
    if not output_file.exists():
        raise RuntimeError(f"ns_ht_2d did not produce {output_file}")
    return TruchasVTKHDFData(output_file), output_dir


def compare_field(actual, reference, name, step, tolerance):
    if actual.shape != reference.shape:
        raise RuntimeError(
            f"{name} shape differs at t={EXPECTED_TIMES[step]:g}: "
            f"{actual.shape} != {reference.shape}"
        )

    actual_finite = np.isfinite(actual)
    reference_finite = np.isfinite(reference)
    if not np.array_equal(actual_finite, reference_finite):
        raise RuntimeError(
            f"{name} finite-value mask differs at t={EXPECTED_TIMES[step]:g}"
        )

    values = np.abs(actual[reference_finite] - reference[reference_finite])
    error = np.max(values) if values.size else 0.0
    if error > tolerance:
        raise RuntimeError(
            f"{name} error at t={EXPECTED_TIMES[step]:g} is {error:g}; "
            f"expected no more than {tolerance:g}"
        )


def check_case(data, reference):
    if data.num_steps != len(EXPECTED_TIMES):
        raise RuntimeError(
            f"expected {len(EXPECTED_TIMES)} output states, found {data.num_steps}"
        )
    if reference.num_steps != len(EXPECTED_TIMES):
        raise RuntimeError(
            f"reference has {reference.num_steps} output states, "
            f"expected {len(EXPECTED_TIMES)}"
        )

    tolerances = {
        "T": (1.0e-3, 2.0e-4),
        "H": (2.0e-3, 2.0e-4),
        "vf_freezing_material_liquid": (3.0e-4, 5.0e-5),
        "vf_freezing_material_solid": (3.0e-4, 5.0e-5),
    }

    for step, expected_time in enumerate(EXPECTED_TIMES):
        if abs(data.time(step) - expected_time) > 1.0e-12:
            raise RuntimeError(
                f"output {step} has time {data.time(step):g}, "
                f"expected {expected_time:g}"
            )
        if abs(reference.time(step) - expected_time) > 1.0e-12:
            raise RuntimeError(
                f"reference output {step} has time {reference.time(step):g}, "
                f"expected {expected_time:g}"
            )

        for name, tolerance in tolerances.items():
            compare_field(
                data.field(step, name),
                reference.field(step, name),
                name,
                step,
                tolerance[0 if step == 1 else 1],
            )

        if step > 0:
            pressure = data.field(step, "pressure")
            velocity = data.field(step, "velocity")
            reference_pressure = reference.field(step, "pressure")
            reference_velocity = reference.field(step, "velocity")
            finite = np.isfinite(pressure) & np.isfinite(velocity).all(axis=1)
            reference_finite = (
                np.isfinite(reference_pressure)
                & np.isfinite(reference_velocity).all(axis=1)
            )
            if not np.array_equal(finite, reference_finite):
                raise RuntimeError(
                    f"active flow-cell mask differs at t={EXPECTED_TIMES[step]:g}"
                )
            if np.max(np.abs(pressure[finite])) > 1.0e-6:
                raise RuntimeError(
                    f"pressure is not analytic zero at t={EXPECTED_TIMES[step]:g}"
                )
            if np.max(np.abs(velocity[finite, 0] - 1.0)) > 1.0e-6:
                raise RuntimeError(
                    f"velocity is not analytic plug flow at t={EXPECTED_TIMES[step]:g}"
                )
            if np.max(np.abs(velocity[finite, 1:])) > 1.0e-6:
                raise RuntimeError(
                    f"transverse velocity is not zero at t={EXPECTED_TIMES[step]:g}"
                )

        liquid = data.field(step, PHASE_FIELDS[0])
        solid = data.field(step, PHASE_FIELDS[1])
        if np.min(liquid) < -1.0e-12 or np.min(solid) < -1.0e-12:
            raise RuntimeError(f"step {step}: negative phase fraction")
        if np.max(np.abs(liquid + solid - 1.0)) > 1.0e-12:
            raise RuntimeError(f"step {step}: phase fractions do not sum to one")

    initial_liquid = data.field(0, "vf_freezing_material_liquid")
    initial_solid = data.field(0, "vf_freezing_material_solid")
    if np.max(np.abs(initial_liquid - 1.0)) > 1.0e-14 or np.max(np.abs(initial_solid)) > 1.0e-14:
        raise RuntimeError("initial material is not entirely liquid")

    final_solid = data.field(len(EXPECTED_TIMES) - 1, "vf_freezing_material_solid")
    if np.max(final_solid) <= 1.0e-2:
        raise RuntimeError("cooling did not produce a solidified region")

    reference_nsteps = [
        int(reference.field(step, "NStep", association="field")[0])
        for step in range(len(EXPECTED_TIMES))
    ]
    actual_nsteps = [
        int(data.field(step, "NStep", association="field")[0])
        for step in range(len(EXPECTED_TIMES))
    ]
    nstep_relative_tolerance = 0.02
    for step, (actual, expected) in enumerate(zip(actual_nsteps, reference_nsteps)):
        if expected == 0:
            error = actual != expected
        else:
            error = abs(actual - expected) / expected > nstep_relative_tolerance
        if error:
            raise RuntimeError(
                f"NStep at t={EXPECTED_TIMES[step]:g} is {actual}; "
                f"reference is {expected} (relative tolerance "
                f"{nstep_relative_tolerance:g})"
            )

    return actual_nsteps


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} TRUCHAS_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable = Path(sys.argv[1]).resolve()
    input_file = Path(sys.argv[2]).resolve()
    reference_file = Path(__file__).with_name("reference") / "out.vtkhdf"
    data, output_dir = run_case(executable, input_file, sys.argv[3])
    reference = TruchasVTKHDFData(reference_file)
    nsteps = check_case(data, reference)

    print("PASS: inviscid freezing plug flow agrees with four-process reference")
    print(f"      accepted steps: {nsteps}")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        sys.exit(1)
