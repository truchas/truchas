#!/usr/bin/env python3

import os
import subprocess

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


def vtkhdf_filename(output, input_name):
    return os.path.join(output.directory, os.path.splitext(input_name)[0] + ".vtkhdf")


def run_truchas_expect_failure(tenv, nprocs, input_file, output_name):
    input_path = os.path.join(tenv.input_directory(), input_file)
    output_dir = os.path.join(tenv.working_directory(), output_name)
    command = [
        tenv._mpiexec, "-n", str(nprocs), tenv._truchas_executable,
        "-f", "-o:" + output_dir, input_path,
    ]
    process = subprocess.run(command, encoding="utf-8", stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    return process.returncode, process.stdout + process.stderr


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(1, "moving-block.inp")
    filename = vtkhdf_filename(output, "moving-block.inp")

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

    stdout, output = tenv.truchas(1, "moving-block-cycle-max.inp")
    filename = vtkhdf_filename(output, "moving-block-cycle-max.inp")

    with h5py.File(filename, "r") as vtkhdf:
        moving = block_by_name(vtkhdf, "moving_block")
        times = moving["Steps/Values"][:]
        nfail += check(len(times) == 2 and numpy.allclose(times, [0.0, 0.25]),
                       "cycle_max rescue write does not duplicate a scheduled VTKHDF step")

    stdout, output = tenv.truchas(1, "moving-block-static-thermal.inp")
    filename = vtkhdf_filename(output, "moving-block-static-thermal.inp")

    with h5py.File(filename, "r") as vtkhdf:
        fixed = block_by_name(vtkhdf, "fixed_block")
        moving = block_by_name(vtkhdf, "moving_block")

        nfail += check("CellData/temperature" in fixed,
                       "fixed static temperature is registered as static data")
        nfail += check("Steps/CellDataOffsets/temperature" not in fixed,
                       "fixed static temperature has no temporal offsets")
        nfail += check("CellData/temperature" in moving,
                       "moving static temperature is registered")
        nfail += check("Steps/CellDataOffsets/temperature" in moving,
                       "moving static temperature is registered as temporal data")

    returncode, output = run_truchas_expect_failure(
        tenv, 1, "moving-block-invalid-id.inp", "moving-block-invalid-id_output")
    nfail += check(returncode != 0,
                   "invalid moving block id fails")
    nfail += check("move-block-ids contains unknown mesh block id 999" in output,
                   "invalid moving block id reports the bad id")

    truchas.report_summary(nfail)
    return nfail


if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
