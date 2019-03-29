#!/usr/bin/env python3

import truchas

def test_read_restart(tenv):
    """Test creating and reading a restart file"""
    # run truchas & write restart
    stdout, output = tenv.truchas(4, "ds11.inp", output_dir="read_test")
    restart_filename = tenv._working_dir + "/read_test_restart/restart.75"
    tenv.write_restart(output.filename, 75, restart_filename)

    # run from restart
    stdout, output = tenv.truchas(4, "ds11.inp", restart_file=restart_filename, \
                                  output_dir="read_test_restart")

    print("PASS: creating and reading a restart file")
    return 0


def test_fields(tenv):
    """Test comparing final field datasets when restarting"""
    out1 = tenv.open_data("read_test/ds11.h5")
    out2 = tenv.open_data("read_test_restart/ds11.h5")

    # check that the data in the last series is identical
    series_id1 = out1.num_series()
    series_id2 = out2.num_series()

    nfail = 0
    for field_name in ("VOF", "Z_ENTHALPY", "Z_TEMP"):
        f1 = out1.field(series_id1, field_name)
        f2 = out2.field(series_id2, field_name)
        err = abs(f1 - f2).max()
        if err > 1e-3:
            print(field_name, " err: ", err)
            nfail = 1

    status = "PASS" if nfail == 0 else "FAIL"
    print("{:s}: comparing final field datasets when restarting".format(status))
    return nfail


def test_repartition(tenv):
    """Test restarting on a different number of processors"""
    # run with 2 processors & write restart
    stdout, output = tenv.truchas(2, "ds11.inp", output_dir="repartition_np2")
    restart_filename = tenv._working_dir + "/repartition_np2/restart.75"
    tenv.write_restart(output.filename, 75, restart_filename)

    # run from restart with 4 processors
    stdout, output = tenv.truchas(4, "ds11.inp", restart_file=restart_filename, \
                                  output_dir="repartition_restart_np4")
    print("PASS: restarting on a different number of processors")
    return 0


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = 0
    nfail += test_read_restart(tenv)
    nfail += test_fields(tenv)
    nfail += test_repartition(tenv)
    assert nfail == 0
