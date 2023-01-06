#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ustruc-ldrd-restart.inp", restart_file="ustruc-ldrd.restart.141")
    golden = tenv.output("ustruc-ldrd_golden/ustruc-ldrd.h5")

    time = output.time(2) # final time

    # temperature
    test = output.field(2, "Z_TEMP")
    gold = golden.field(3, "Z_TEMP")
    nfail += truchas.compare_max_rel(test, gold, 5e-5, "temp", time)

    # G
    atol = 200
    rtol = 1e-4
    test = np.ma.masked_values(output.field(2, "ustruc-G"), 0)
    gold = np.ma.masked_values(golden.field(3, "ustruc-G"), 0)
    error = abs((test-gold) / (atol + rtol*gold))
    nfail += truchas.compare_max(error, 0, 1.0, "G", time)

    # V
    atol = 1e-4
    rtol = 1e-4
    test = np.ma.masked_values(output.field(2, "ustruc-V"), 0)
    gold = np.ma.masked_values(golden.field(3, "ustruc-V"), 0)
    error = abs((test-gold) / (atol + rtol* gold))
    nfail += truchas.compare_max(error, 0, 1.0, "V", time)

    # lambda1
    atol = 1e-5
    rtol = 1e-4
    test = np.ma.masked_values(output.field(2, "ustruc-lambda1"), 0)
    gold = np.ma.masked_values(golden.field(3, "ustruc-lambda1"), 0)
    error = abs((test-gold) / (atol + rtol*gold))
    nfail += truchas.compare_max(error, 0, 1.0, "lambda1", time)

    # lambda2
    atol = 2e-6
    rtol = 1e-4
    test = np.ma.masked_values(output.field(2, "ustruc-lambda2"), 0)
    gold = np.ma.masked_values(golden.field(3, "ustruc-lambda2"), 0)
    error = abs((test-gold) / (atol + rtol*gold))
    nfail += truchas.compare_max(error, 0, 1.0, "lambda2", time)

    # ustruc
    test = output.field(2, "ustruc-type")
    gold = golden.field(3, "ustruc-type")
    success = (test == gold).all()
    print("{:s}: {:s} at t={:8.2e}".format("PASS" if success else "FAIL", "ustruc", time))
    if not success: nfail += 1

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
