#!/usr/bin/env python3

"""Analytic checks shared by the axial gravity fill and drain tests."""

import numpy as np


def clip_polygon_x(polygon, x_cut, keep_left):
    clipped = []
    for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
        start_inside = start[0] <= x_cut if keep_left else start[0] >= x_cut
        end_inside = end[0] <= x_cut if keep_left else end[0] >= x_cut
        if start_inside:
            clipped.append(start)
        if start_inside != end_inside:
            fraction = (x_cut - start[0]) / (end[0] - start[0])
            clipped.append(start + fraction * (end - start))
    return np.asarray(clipped)


def polygon_area(polygon):
    if len(polygon) < 3:
        return 0.0
    return 0.5 * abs(
        np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1))
        - np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1))
    )


def expected_water(data, step, front):
    values = []
    areas = []
    for vertices in data.cell_vertices(step):
        polygon = vertices[:, :2]
        area = polygon_area(polygon)
        clipped = clip_polygon_x(polygon, front, keep_left=True)
        values.append(polygon_area(clipped) / area if len(clipped) else 0.0)
        areas.append(area)
    return np.asarray(values), np.asarray(areas)


def check_solution(data, initial_front, acceleration, step_dt):
    expected_times = (0.0, 0.25, 0.5, 0.75, 1.0)
    if data.num_steps != len(expected_times):
        raise RuntimeError(f"expected {len(expected_times)} output states, got {data.num_steps}")

    for step, time in enumerate(expected_times):
        observed_time = data.time(step)
        if abs(observed_time - time) > 1.0e-14:
            raise RuntimeError(f"step {step}: time={observed_time:g}, expected {time:g}")

        # Material transport precedes the momentum update, so each fixed
        # step uses u(t_n).  This is the exact trajectory of that split
        # scheme, rather than the continuum trajectory.
        front = initial_front + 0.5 * acceleration * time * (time - step_dt)
        water = data.field(step, "vf_water")
        expected, areas = expected_water(data, step, front)
        difference = water - expected
        expected_volume = np.sum(areas * expected)
        water_error = np.sum(areas * np.abs(difference)) / expected_volume
        conservation_error = abs(np.sum(areas * difference)) / expected_volume
        if water_error > 2.0e-2:
            raise RuntimeError(f"step {step}: water-fraction L1 error={water_error:g}")
        if conservation_error > 5.0e-4:
            raise RuntimeError(f"step {step}: water-volume conservation error={conservation_error:g}")

        velocity = data.field(step, "velocity")
        active = np.isfinite(velocity[:, 0])
        velocity_error = np.max(np.abs(velocity[active, :2] - [acceleration * time, 0.0]))
        if velocity_error > 1.0e-6:
            raise RuntimeError(f"step {step}: velocity error={velocity_error:g}")

        pressure = data.field(step, "pressure")
        active = np.isfinite(pressure)
        pressure_error = np.max(np.abs(pressure[active]))
        if pressure_error > 1.0e-10:
            raise RuntimeError(f"step {step}: pressure error={pressure_error:g}")


def run_test(executable, input_file, mpiexec, initial_front, acceleration, step_dt, label):
    import pathlib
    import subprocess
    import tempfile
    from TruchasVTKHDFData import TruchasVTKHDFData

    output_dir = pathlib.Path(tempfile.mkdtemp(prefix=f"flow_free_surface_gravity_{label}_4p_"))
    result = subprocess.run(
        [mpiexec, "-n", "4", str(executable), "--simulation", "flow", "--output-dir", ".", "--force", str(input_file)],
        cwd=output_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        print(result.stdout, end="")
        raise RuntimeError(f"flow returned {result.returncode}")
    data = TruchasVTKHDFData(output_dir / "out.vtkhdf")
    check_solution(data, initial_front, acceleration, step_dt)
    print(f"PASS: gravity {label} matches the analytic solution")
    print(f"      artifacts: {output_dir}")
