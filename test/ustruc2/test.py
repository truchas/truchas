#!/usr/bin/env python3

import scipy as sp

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ustruc2.inp", restart_file="ustruc1.restart.177")
    golden = tenv.output("ustruc2_golden/ustruc2.h5")

    time = output.time(2) # final time

    # temperature
    test = output.field(2, "Z_TEMP")
    gold = golden.field(3, "Z_TEMP")
    nfail += truchas.compare_max_rel(test, gold, 5e-5, "temp", time)

    # G
    atol = 200
    rtol = 1e-4
    test = sp.ma.masked_values(output.field(2, "uStruc-G"), 0)
    gold = sp.ma.masked_values(golden.field(3, "uStruc-G"), 0)
    error = abs((test-gold) / (atol + rtol*gold))
    nfail += truchas.compare_max(error, 0, 1.0, "G", time)

    # V
    atol = 1e-4
    rtol = 1e-4
    test = sp.ma.masked_values(output.field(2, "uStruc-V"), 0)
    gold = sp.ma.masked_values(golden.field(3, "uStruc-V"), 0)
    error = abs((test-gold) / (atol + rtol* gold))
    nfail += truchas.compare_max(error, 0, 1.0, "V", time)

    # lambda1
    atol = 1e-5
    rtol = 1e-4
    test = sp.ma.masked_values(output.field(2, "uStruc-gv1-lambda1"), 0)
    gold = sp.ma.masked_values(golden.field(3, "uStruc-gv1-lambda1"), 0)
    error = abs((test-gold) / (atol + rtol*gold))
    nfail += truchas.compare_max(error, 0, 1.0, "lambda1", time)

    # lambda2
    atol = 2e-6
    rtol = 1e-4
    test = sp.ma.masked_values(output.field(2, "uStruc-gv1-lambda2"), 0)
    gold = sp.ma.masked_values(golden.field(3, "uStruc-gv1-lambda2"), 0)
    error = abs((test-gold) / (atol + rtol*gold))
    nfail += truchas.compare_max(error, 0, 1.0, "lambda2", time)

    # ustruc
    test = output.field(2, "uStruc-gv1-ustruc")
    gold = golden.field(3, "uStruc-gv1-ustruc")
    success = (test == gold).all()
    print("{:s}: {:s} at t={:8.2e}".format("PASS" if success else "FAIL", "ustruc", time))
    if not success: nfail += 1

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
