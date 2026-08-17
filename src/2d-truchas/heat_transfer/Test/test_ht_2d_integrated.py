#!/usr/bin/env python3

"""Small end-to-end tests for the JSON-driven ht_2d application."""

import pathlib
import re
import subprocess
import sys
import tempfile

import numpy as np

source_root = pathlib.Path(__file__).resolve().parents[4]
sys.path.insert(0, str(source_root / "src/python/truchas"))
from TruchasVTKHDFData import TruchasVTKHDFData


def half_plane_fraction(vertices, point, normal):
    """Return the fraction of a convex polygon in the selected half-plane."""
    vertices = vertices[:, :2]
    clipped = []
    for start, end in zip(vertices, np.roll(vertices, -1, axis=0)):
        start_value = np.dot(normal, start - point)
        end_value = np.dot(normal, end - point)
        if start_value <= 0.0:
            clipped.append(start)
        if (start_value <= 0.0) != (end_value <= 0.0):
            clipped.append(start + start_value / (start_value - end_value) * (end - start))

    def area(polygon):
        if len(polygon) < 3:
            return 0.0
        return abs(np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1)) -
                   np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1))) / 2.0

    return area(np.asarray(clipped)) / area(vertices)


def main():
    if len(sys.argv) not in (4, 6):
        print(
            f"usage: {sys.argv[0]} HT_2D JSON_INPUT CASE [NPROC MPIEXEC]",
            file=sys.stderr,
        )
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    case = sys.argv[3]
    nproc = int(sys.argv[4]) if len(sys.argv) == 6 else 1
    mpiexec = sys.argv[5] if len(sys.argv) == 6 else None

    output_dir = pathlib.Path(
        tempfile.mkdtemp(prefix=f"ht_2d_{case}_{nproc}p_")
    )
    try:
        command = [str(executable), str(input_file)]
        if nproc > 1:
            command = [mpiexec, "-n", str(nproc)] + command
        result = subprocess.run(
            command,
            cwd=output_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        (output_dir / "stdout.txt").write_text(result.stdout)

        if result.returncode != 0:
            print(result.stdout, end="")
            print(f"FAIL: ht_2d returned {result.returncode}")
            return 1

        log_file = output_dir / "run.log"
        if not log_file.exists():
            print("FAIL: ht_2d did not produce run.log")
            return 1
        log = log_file.read_text()

        if "unrecoverable integration failure" in log:
            print("FAIL: unrecoverable integration failure")
            return 1
        if re.search(r"\b(?:nan|inf)\b", log, re.IGNORECASE):
            print("FAIL: non-finite value reported in run.log")
            return 1

        completed = re.findall(r"Completed integration to T =\s*([0-9.E+-]+)", log)
        if not completed:
            print("FAIL: completion record not found in run.log")
            return 1
        final_time = float(completed[-1])
        if abs(final_time - 1.0e-2) > 1.0e-12:
            print(f"FAIL: final time {final_time:g} differs from 0.01")
            return 1

        output_file = output_dir / "out.vtkhdf"
        if not output_file.exists():
            print("FAIL: ht_2d did not produce out.vtkhdf")
            return 1
        vtkhdf = TruchasVTKHDFData(output_file)
        final_step = vtkhdf.num_steps - 1
        if abs(vtkhdf.time(final_step) - final_time) > 1.0e-12:
            print("FAIL: VTKHDF final time disagrees with run.log")
            return 1
        temperature = vtkhdf.field(final_step, "T")

        expected_temperature = {
            "uniform": 2.0,
            "source": 2.01,
            "nonlinear": (1.0 + 2.0e-2) ** 0.5 - 1.0,
        }
        if case in expected_temperature:
            expected = expected_temperature[case]
            error = max(abs(value - expected) for value in temperature)
            tolerance = 1.0e-8 if case == "nonlinear" else (1.0e-9 if case == "source" else 1.0e-10)
            if error > tolerance:
                print(f"FAIL: {case} temperature error {error:g}")
                return 1
        elif case == "linear":
            expected = 1.0 - vtkhdf.cell_centers(final_step)[:, 0]
            error = max(abs(temperature - expected))
            if error > 1.0e-10:
                print(f"FAIL: steady linear profile error {error:g}")
                return 1
        elif case == "multimaterial_source":
            x = vtkhdf.cell_centers(final_step)[:, 0]
            expected = np.where(x <= 0.5, 2.0 + 1.0e-2, 2.0 + 1.0e-2 / 6.0)
            error = max(abs(temperature - expected))
            if error > 1.0e-9:
                print(f"FAIL: mixed material source temperature error {error:g}")
                return 1
        elif case == "multimaterial_mixed_source":
            x = vtkhdf.cell_centers(final_step)[:, 0]
            low_fraction = np.clip((0.4375 - (x - 0.0625)) / 0.125, 0.0, 1.0)
            heat_capacity = low_fraction + 3.0 * (1.0 - low_fraction)
            expected = 2.0 + 1.0e-2 / heat_capacity
            error = max(abs(temperature - expected))
            if error > 1.0e-9:
                print(f"FAIL: mixed-cell source temperature error {error:g}")
                return 1
        elif case == "multimaterial_mixed_tri_source":
            vertices = vtkhdf.cell_vertices(final_step)
            if not any(len(cell) == 3 for cell in vertices) or \
                    not any(len(cell) == 4 for cell in vertices):
                print("FAIL: mixed-cell source mesh does not contain triangles and quads")
                return 1
            point = np.array([0.53, 0.0])
            normal = np.array([1.0, 0.6])
            low_fraction = np.array(
                [half_plane_fraction(cell, point, normal) for cell in vertices]
            )
            if not np.any((low_fraction > 0.0) & (low_fraction < 1.0)):
                print("FAIL: mixed-cell source mesh has no mixed cells")
                return 1
            heat_capacity = low_fraction + 3.0 * (1.0 - low_fraction)
            expected = 2.0 + 1.0e-2 / heat_capacity
            error = max(abs(temperature - expected))
            if error > 1.0e-5:
                print(f"FAIL: mixed triangle/quad source temperature error {error:g}")
                return 1
        elif case == "multimaterial_conduction":
            flux = 1.0 / (0.5 / 1.0 + 0.5 / 10.0)
            x = vtkhdf.cell_centers(final_step)[:, 0]
            expected = np.where(
                x <= 0.5,
                1.0 - flux * x,
                1.0 - 0.5 * flux - flux / 10.0 * (x - 0.5),
            )
            error = max(abs(temperature - expected))
            if error > 1.0e-9:
                print(f"FAIL: mixed material steady profile error {error:g}")
                return 1

        print(f"PASS: {case} reached final time {final_time:g}")
        print(f"      artifacts: {output_dir}")
        return 0

    finally:
        # Keep the run directory available for post-test diagnosis.
        pass


if __name__ == "__main__":
    sys.exit(main())
