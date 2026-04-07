#!/usr/bin/env python3

import os

import h5py
import numpy
import truchas


def block_by_name(vtkhdf, name):
    for key in vtkhdf["VTKHDF"]:
        if not key.startswith("Group_"):
            continue
        block = vtkhdf["VTKHDF"][key]
        block_name = block.attrs["Name"]
        if isinstance(block_name, bytes):
            block_name = block_name.decode()
        if block_name == name:
            return block
    raise KeyError(name)


def check(condition, label):
    print(("PASS" if condition else "FAIL") + ": " + label)
    return 0 if condition else 1


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(1, "moving-block.inp")
    filename = os.path.join(output.directory, "moving-block.vtkhdf")

    with h5py.File(filename, "r") as vtkhdf:
        fixed = block_by_name(vtkhdf, "fixed_block")
        moving = block_by_name(vtkhdf, "moving_block")

        fixed_static_mesh = fixed["Steps"].attrs["StaticMesh"]
        moving_static_mesh = moving["Steps"].attrs["StaticMesh"]
        nfail += check(fixed_static_mesh == 1, "fixed block has a static mesh")
        nfail += check(moving_static_mesh == 0, "moving block has a moving mesh")

        times = moving["Steps/Values"][:]
        point_offsets = moving["Steps/PointOffsets"][:]
        point_count = moving["NumberOfPoints"][0]
        point_rows = moving["Points"].shape[0]

        stationary = times <= 0.5 + 1.0e-12
        moved = times > 0.5 + 1.0e-12
        nfail += check(numpy.all(point_offsets[stationary] == 0),
                       "stationary moving PointOffsets are reused")
        nfail += check(numpy.all(numpy.diff(point_offsets) >= 0),
                       "moving PointOffsets are monotone")
        nfail += check(numpy.any(numpy.diff(point_offsets[moved]) > 0),
                       "moving PointOffsets advance during motion")
        nfail += check(point_rows == max(point_offsets) + point_count,
                       "moving Points length matches written geometries")

    truchas.report_summary(nfail)
    return nfail


if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
