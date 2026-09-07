#!/usr/bin/env python3

"""Analytic four-process regression for a perturbed fluid/VOID plug."""

import numpy as np

from test_pipe import run_case


def clip_polygon_x(polygon, x_cut, keep_right):
    clipped = []
    for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
        start_offset = start[0] - x_cut
        end_offset = end[0] - x_cut
        start_inside = start_offset >= 0.0 if keep_right else start_offset <= 0.0
        end_inside = end_offset >= 0.0 if keep_right else end_offset <= 0.0
        if start_inside:
            clipped.append(start)
        if start_inside != end_inside:
            clipped.append(start + start_offset / (start_offset - end_offset) * (end - start))
    return np.asarray(clipped)


def polygon_area(polygon):
    return 0.5 * abs(
        np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1))
        - np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1))
    )


def expected_water_volume_fraction(data, step, time):
    lower, upper = 0.5 + time, 1.5 + time
    result = []
    areas = []
    for vertices in data.cell_vertices(step):
        polygon = vertices[:, :2]
        cell_area = polygon_area(polygon)
        areas.append(cell_area)
        clipped = clip_polygon_x(polygon, lower, keep_right=True)
        if len(clipped):
            clipped = clip_polygon_x(clipped, upper, keep_right=False)
        result.append(polygon_area(clipped) / cell_area if len(clipped) else 0.0)
    return np.asarray(result), np.asarray(areas)


def check_solution(data):
    expected_times = (0.0, 1.0, 3.0)
    if data.num_steps != len(expected_times):
        raise RuntimeError(
            f"expected {len(expected_times)} output states, got {data.num_steps}"
        )

    for step, time in enumerate(expected_times):
        observed_time = data.time(step)
        if abs(observed_time - time) > 1.0e-14:
            raise RuntimeError(f"step {step}: time={observed_time:g}, expected {time:g}")

        water = data.field(step, "vf_water")
        expected, cell_areas = expected_water_volume_fraction(data, step, time)
        difference = water - expected
        expected_volume = np.sum(cell_areas * expected)
        water_error = np.sum(cell_areas * np.abs(difference)) / expected_volume
        conservation_error = abs(np.sum(cell_areas * difference)) / expected_volume
        # The perturbed coarse mesh does not reproduce pure translation to
        # roundoff: this checks local interface-reconstruction error, unlike
        # the separate signed conservation check below.
        if water_error > 1.0e-2:
            raise RuntimeError(
                f"step {step}: water-volume-fraction L1 error={water_error:g}"
            )
        if conservation_error > 1.0e-4:
            raise RuntimeError(
                f"step {step}: water-volume conservation error={conservation_error:g}"
            )

        velocity = data.field(step, "velocity")
        active = np.isfinite(velocity[:, 0])
        velocity_error = np.max(np.abs(velocity[active, :2] - [1.0, 0.0]))
        if velocity_error > 1.0e-4:
            raise RuntimeError(f"step {step}: velocity error={velocity_error:g}")

        pressure = data.field(step, "pressure")
        active = np.isfinite(pressure)
        pressure_error = np.max(np.abs(pressure[active]))
        if pressure_error > 1.0e-3:
            raise RuntimeError(f"step {step}: pressure error={pressure_error:g}")


def main():
    import pathlib
    import sys

    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} TRUCHAS_2D JSON_INPUT MPIEXEC", file=sys.stderr)
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    try:
        data, output_dir = run_case(executable, input_file, sys.argv[3])
        check_solution(data)
    except (RuntimeError, ValueError) as error:
        print(f"FAIL: {error}")
        return 1

    print("PASS: noisy free-surface plug matches the analytic solution")
    print(f"      artifacts: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
