#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import check_close, execute, finish, run_from_argv


def half_plane_fraction(vertices, point, normal):
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


def test():
    run = run_from_argv(Path(__file__).with_name("input.json"))
    vertices = run.data.cell_vertices(run.final_step)
    if not any(len(cell) == 3 for cell in vertices) or not any(len(cell) == 4 for cell in vertices):
        raise AssertionError("mesh does not contain both triangles and quads")
    point = np.array([0.53, 0.0])
    normal = np.array([1.0, 0.6])
    low_fraction = np.array([half_plane_fraction(cell, point, normal) for cell in vertices])
    if not np.any((low_fraction > 0.0) & (low_fraction < 1.0)):
        raise AssertionError("mesh has no mixed cells")
    heat_capacity = low_fraction + 3.0 * (1.0 - low_fraction)
    expected = 2.0 + 1.0e-2 / heat_capacity
    check_close(run.data.field(run.final_step, "T"), expected, 1.0e-5, "mixed triangle/quad source")
    finish("multimaterial mixed triangle/quad source", run)


if __name__ == "__main__":
    sys.exit(execute(test))
