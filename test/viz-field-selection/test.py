#!/usr/bin/env python3

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import h5py
import truchas


TEST_ROOT = Path(__file__).resolve().parents[1]


def check(condition, label):
    print(("PASS" if condition else "FAIL") + ": " + label)
    return 0 if condition else 1


def write_input_with_fields(src, dst, fields):
    text = src.read_text()
    line = "  fields = " + ", ".join(f"'{field}'" for field in fields) + "\n"
    dst.write_text(text.replace("&OUTPUTS\n", "&OUTPUTS\n" + line, 1))


def write_input_with_replacements_and_fields(src, dst, replacements, fields):
    text = src.read_text()
    for old, new in replacements:
        text = text.replace(old, new)
    line = "  fields = " + ", ".join(f"'{field}'" for field in fields) + "\n"
    dst.write_text(text.replace("&OUTPUTS\n", "&OUTPUTS\n" + line, 1))


def prepare_case(work_dir, name, source_dir, input_name, mesh_files, fields):
    case_dir = work_dir / name
    case_dir.mkdir()

    source_dir = TEST_ROOT / source_dir
    for mesh_file in mesh_files:
        shutil.copy2(source_dir / mesh_file, case_dir / mesh_file)

    input_path = case_dir / input_name
    if fields is None:
        shutil.copy2(source_dir / input_name, input_path)
    else:
        write_input_with_fields(source_dir / input_name, input_path, fields)
    return case_dir, input_path


def prepare_case_with_replacements(work_dir, name, source_dir, input_name,
                                   mesh_files, replacements, fields):
    case_dir = work_dir / name
    case_dir.mkdir()

    source_dir = TEST_ROOT / source_dir
    for mesh_file in mesh_files:
        shutil.copy2(source_dir / mesh_file, case_dir / mesh_file)

    input_path = case_dir / input_name
    write_input_with_replacements_and_fields(source_dir / input_name, input_path,
                                             replacements, fields)
    return case_dir, input_path


def vtkhdf_blocks(vtkhdf):
    for key, block in vtkhdf["VTKHDF"].items():
        if key.startswith("Group_"):
            yield block


def block_name(block):
    name = block.attrs["Name"]
    return name.decode() if isinstance(name, bytes) else name


def block_by_name(vtkhdf, name):
    for block in vtkhdf_blocks(vtkhdf):
        if block_name(block) == name:
            return block
    raise KeyError(name)


def has_dataset(vtkhdf, path):
    return any(path in block for block in vtkhdf_blocks(vtkhdf))


def vtkhdf_filename(output, input_path):
    return Path(output.directory) / (input_path.stem + ".vtkhdf")


def run_truchas_raw(tenv, nprocs, input_path, output_dir):
    command = "{:s} -n {:d} {:s} -o:{:s} -f {:s}" \
        .format(tenv._mpiexec, nprocs, tenv._truchas_executable,
                str(output_dir), str(input_path))
    print(command)
    return subprocess.run(command, shell=True, universal_newlines=True,
                          encoding='utf-8',
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def check_flow_pressure_only(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "flow-pressure-only", "hydrostatic", "hydrostatic-1a.inp",
        ["mesh1.gen"], ["pressure"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(has_dataset(vtkhdf, "CellData/pressure"),
                       "selected pressure field is registered")
        nfail += check(has_dataset(vtkhdf, "Steps/CellDataOffsets/pressure"),
                       "selected pressure offsets are registered")
        nfail += check(not has_dataset(vtkhdf, "CellData/velocity"),
                       "unselected velocity field is omitted")
        nfail += check(not has_dataset(vtkhdf, "Steps/CellDataOffsets/velocity"),
                       "unselected velocity offsets are omitted")
        nfail += check(not has_dataset(vtkhdf, "CellData/ProcessIds"),
                       "unselected process id field is omitted")

    return nfail


def check_static_temperature_only(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "static-temperature-only", "hydrostatic", "hydrostatic-1a.inp",
        ["mesh1.gen"], ["temperature"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(has_dataset(vtkhdf, "CellData/temperature"),
                       "selected static temperature field is registered")
        nfail += check(not has_dataset(vtkhdf, "Steps/CellDataOffsets/temperature"),
                       "selected static temperature has no temporal offsets")
        nfail += check(not has_dataset(vtkhdf, "CellData/enthalpy"),
                       "unselected enthalpy field is omitted")
        nfail += check(not has_dataset(vtkhdf, "CellData/pressure"),
                       "unselected pressure field is omitted")
        nfail += check(not has_dataset(vtkhdf, "CellData/velocity"),
                       "unselected velocity field is omitted")

    return nfail


def check_volume_fraction_only(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "volume-fraction-only", "hydrostatic", "hydrostatic-3a.inp",
        ["mesh1.gen"], ["vfrac_water"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(has_dataset(vtkhdf, "CellData/vfrac_water"),
                       "selected water volume fraction field is registered")
        nfail += check(has_dataset(vtkhdf, "Steps/CellDataOffsets/vfrac_water"),
                       "selected water volume fraction offsets are registered")
        nfail += check(not has_dataset(vtkhdf, "CellData/vfrac_oil"),
                       "unselected oil volume fraction field is omitted")
        nfail += check(not has_dataset(vtkhdf, "CellData/pressure"),
                       "unselected pressure field is omitted")
        nfail += check(not has_dataset(vtkhdf, "CellData/velocity"),
                       "unselected velocity field is omitted")

    return nfail


def check_duplicate_volume_fraction_selection(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "duplicate-volume-fraction-selection", "hydrostatic",
        "hydrostatic-3a.inp", ["mesh1.gen"], ["vfrac_water", "vfrac_water"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(has_dataset(vtkhdf, "CellData/vfrac_water"),
                       "duplicate selected volume fraction field is registered")
        nfail += check(has_dataset(vtkhdf, "Steps/CellDataOffsets/vfrac_water"),
                       "duplicate selected volume fraction offsets are registered")
        nfail += check(not has_dataset(vtkhdf, "CellData/vfrac_oil"),
                       "duplicate selection omits unselected oil volume fraction")

    return nfail


def check_volume_fraction_preserves_phase_name_case(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case_with_replacements(
        work_dir, "volume-fraction-preserves-phase-name-case", "hydrostatic",
        "hydrostatic-3a.inp", ["mesh1.gen"],
        [("'water'", "'Water'"), ("'oil'", "'Oil'")], ["vfrac_Water"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(has_dataset(vtkhdf, "CellData/vfrac_Water"),
                       "volume fraction field preserves phase name case")
        nfail += check(has_dataset(vtkhdf, "Steps/CellDataOffsets/vfrac_Water"),
                       "volume fraction offsets preserve phase name case")
        nfail += check(not has_dataset(vtkhdf, "CellData/vfrac_water"),
                       "lower-case volume fraction field is not registered")
        nfail += check(not has_dataset(vtkhdf, "CellData/vfrac_Oil"),
                       "unselected mixed-case oil volume fraction is omitted")

    return nfail


def check_default_volume_fractions(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "default-volume-fractions", "hydrostatic", "hydrostatic-3a.inp",
        ["mesh1.gen"], None)
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(has_dataset(vtkhdf, "CellData/vfrac_water"),
                       "default output includes water volume fraction")
        nfail += check(has_dataset(vtkhdf, "CellData/vfrac_oil"),
                       "default output includes oil volume fraction")
        nfail += check(has_dataset(vtkhdf, "CellData/pressure"),
                       "default output includes pressure")
        nfail += check(has_dataset(vtkhdf, "CellData/velocity"),
                       "default output includes velocity")
        nfail += check(has_dataset(vtkhdf, "CellData/temperature"),
                       "default output includes temperature")
        nfail += check(has_dataset(vtkhdf, "CellData/enthalpy"),
                       "default output includes enthalpy")
        nfail += check(has_dataset(vtkhdf, "CellData/ProcessIds"),
                       "default output includes process ids")

    return nfail


def check_all_selector(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "all-selector", "hydrostatic", "hydrostatic-3a.inp",
        ["mesh1.gen"], ["all"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(has_dataset(vtkhdf, "CellData/vfrac_water"),
                       "all selector includes water volume fraction")
        nfail += check(has_dataset(vtkhdf, "CellData/vfrac_oil"),
                       "all selector includes oil volume fraction")
        nfail += check(has_dataset(vtkhdf, "CellData/pressure"),
                       "all selector includes pressure")
        nfail += check(has_dataset(vtkhdf, "CellData/velocity"),
                       "all selector includes velocity")
        nfail += check(has_dataset(vtkhdf, "CellData/temperature"),
                       "all selector includes temperature")
        nfail += check(has_dataset(vtkhdf, "CellData/enthalpy"),
                       "all selector includes enthalpy")
        nfail += check(has_dataset(vtkhdf, "CellData/ProcessIds"),
                       "all selector includes process ids")

    return nfail


def check_process_id_selection(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "process-id-selection", "hydrostatic", "hydrostatic-1a.inp",
        ["mesh1.gen"], ["process_id"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(has_dataset(vtkhdf, "CellData/ProcessIds"),
                       "process_id field selects process ids")
        nfail += check(not has_dataset(vtkhdf, "CellData/pressure"),
                       "process_id field omits unselected pressure field")
        nfail += check(not has_dataset(vtkhdf, "CellData/velocity"),
                       "process_id field omits unselected velocity field")

    return nfail


def check_process_id_selection_serial(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "process-id-selection-serial", "hydrostatic", "hydrostatic-1a.inp",
        ["mesh1.gen"], ["process_id"])
    stdout, output = tenv.truchas(1, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(not has_dataset(vtkhdf, "CellData/ProcessIds"),
                       "process_id field is omitted in serial output")
        nfail += check(not has_dataset(vtkhdf, "CellData/pressure"),
                       "serial process_id field omits unselected pressure field")
        nfail += check(not has_dataset(vtkhdf, "CellData/velocity"),
                       "serial process_id field omits unselected velocity field")

    return nfail


def check_all_selector_rejects_mixed_fields(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "all-selector-mixed", "hydrostatic", "hydrostatic-1a.inp",
        ["mesh1.gen"], ["all", "temperature"])
    process = run_truchas_raw(tenv, 4, input_path, case_dir / "output")
    output = process.stdout + process.stderr

    nfail += check("may not be combined with other field names" in output,
                   "all selector mixed with field names is rejected")

    return nfail


def check_volume_fraction_shorthand(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "volume-fraction-shorthand", "hydrostatic", "hydrostatic-3a.inp",
        ["mesh1.gen"], ["vfrac"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(has_dataset(vtkhdf, "CellData/vfrac_water"),
                       "vfrac shorthand includes water volume fraction")
        nfail += check(has_dataset(vtkhdf, "CellData/vfrac_oil"),
                       "vfrac shorthand includes oil volume fraction")
        nfail += check(not has_dataset(vtkhdf, "CellData/pressure"),
                       "vfrac shorthand omits unselected pressure field")
        nfail += check(not has_dataset(vtkhdf, "CellData/velocity"),
                       "vfrac shorthand omits unselected velocity field")

    return nfail


def check_gap_displacement_only(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "gap-displacement-only", "contact-box-open",
        "contact-box-open.inp", ["rotated-cube.exo"], ["contact_normal_displ"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        gap = block_by_name(vtkhdf, "gap-interface-10")
        nfail += check("PointData/contact_normal_displ" in gap,
                       "selected gap displacement field is registered")
        nfail += check("Steps/PointDataOffsets/contact_normal_displ" in gap,
                       "selected gap displacement offsets are registered")
        nfail += check("PointData/contact_pressure" not in gap,
                       "unselected gap normal traction field is omitted")
        nfail += check(not has_dataset(vtkhdf, "PointData/displacement"),
                       "unselected solid displacement field is omitted")
        nfail += check(not has_dataset(vtkhdf, "CellData/strain"),
                       "unselected solid strain field is omitted")
        nfail += check(not has_dataset(vtkhdf, "CellData/stress"),
                       "unselected solid stress field is omitted")

    return nfail


def check_inactive_provider_warning(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "inactive-provider", "hydrostatic", "hydrostatic-1a.inp",
        ["mesh1.gen"], ["displacement"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    warning = 'requested visualization field "displacement" is unavailable'
    nfail += check(warning in stdout,
                   "inactive solid mechanics field emits warning")

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(not has_dataset(vtkhdf, "PointData/displacement"),
                       "inactive solid mechanics field is omitted")
        nfail += check(not has_dataset(vtkhdf, "CellData/pressure"),
                       "unrequested active-provider pressure field is omitted")
        nfail += check(not has_dataset(vtkhdf, "CellData/velocity"),
                       "unrequested active-provider velocity field is omitted")

    return nfail


def check_ustruc_in_situ_output(tenv, work_dir):
    nfail = 0
    case_dir, input_path = prepare_case(
        work_dir, "ustruc-in-situ-output", "ustruc", "ustruc-gl-temp.inp",
        ["wedge.exo"], ["temperature", "ustruc"])
    stdout, output = tenv.truchas(4, str(input_path),
                                  output_dir=str(case_dir / "output"))

    with h5py.File(vtkhdf_filename(output, input_path), "r") as vtkhdf:
        nfail += check(has_dataset(vtkhdf, "CellData/temperature"),
                       "selected temperature field is registered with ustruc")
        nfail += check(has_dataset(vtkhdf, "CellData/ustruc1_gl_g"),
                       "active ustruc vector field is registered")
        nfail += check(has_dataset(vtkhdf, "CellData/ustruc1_gl_l"),
                       "active ustruc length field is registered")
        nfail += check(has_dataset(vtkhdf, "CellData/ustruc1_gl_t_sol"),
                       "active ustruc solidification-time field is registered")
        nfail += check(not has_dataset(vtkhdf, "CellData/enthalpy"),
                       "unselected diffusion field is omitted with ustruc")

    return nfail


def run_test(tenv):
    nfail = 0
    with tempfile.TemporaryDirectory(prefix="vtkhdf-field-selection-",
                                     dir=os.getcwd()) as tmp:
        work_dir = Path(tmp)
        nfail += check_flow_pressure_only(tenv, work_dir)
        nfail += check_static_temperature_only(tenv, work_dir)
        nfail += check_volume_fraction_only(tenv, work_dir)
        nfail += check_duplicate_volume_fraction_selection(tenv, work_dir)
        nfail += check_volume_fraction_preserves_phase_name_case(tenv, work_dir)
        nfail += check_default_volume_fractions(tenv, work_dir)
        nfail += check_all_selector(tenv, work_dir)
        nfail += check_process_id_selection(tenv, work_dir)
        nfail += check_process_id_selection_serial(tenv, work_dir)
        nfail += check_all_selector_rejects_mixed_fields(tenv, work_dir)
        nfail += check_volume_fraction_shorthand(tenv, work_dir)
        nfail += check_gap_displacement_only(tenv, work_dir)
        nfail += check_inactive_provider_warning(tenv, work_dir)
        nfail += check_ustruc_in_situ_output(tenv, work_dir)

    truchas.report_summary(nfail)
    return nfail


if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
