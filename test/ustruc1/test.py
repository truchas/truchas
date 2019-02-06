#!/usr/bin/env python3

import scipy as sp

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ustruc1.inp")
    golden = tenv.output("ustruc1_pgolden/ustruc1.h5")

    # final time
    sid = 3
    time = output.time(sid)

    # temperature
    test = output.field(sid, "Z_TEMP")
    gold = golden.field(sid, "Z_TEMP")
    nfail += truchas.compare_max_rel(test, gold, 5e-5, "temp", time)

    # G
    tol = 1e-4
    test = sp.ma.masked_values(output.field(sid, "uStruc-G"), 0)
    gold = sp.ma.masked_values(golden.field(sid, "uStruc-G"), 0)
    error = abs((test-gold) / (200/tol + gold))
    nfail += truchas.compare_max(error, 0, tol, "G", time)

    # V
    tol = 1e-4
    test = sp.ma.masked_values(output.field(sid, "uStruc-V"), 0)
    gold = sp.ma.masked_values(golden.field(sid, "uStruc-V"), 0)
    error = abs((test-gold) / (1e-4/tol + gold))
    nfail += truchas.compare_max(error, 0, tol, "V", time)

    # lambda1
    tol = 1e-4
    test = sp.ma.masked_values(output.field(sid, "uStruc-gv1-lambda1"), 0)
    gold = sp.ma.masked_values(golden.field(sid, "uStruc-gv1-lambda1"), 0)
    error = abs((test-gold) / (2e-6/tol + gold))
    nfail += truchas.compare_max(error, 0, tol, "lambda1", time)

    # lambda2
    tol = 1e-4
    test = sp.ma.masked_values(output.field(sid, "uStruc-gv1-lambda2"), 0)
    gold = sp.ma.masked_values(golden.field(sid, "uStruc-gv1-lambda2"), 0)
    error = abs((test-gold) / (2e-6/tol + gold))
    nfail += truchas.compare_max(error, 0, tol, "lambda2", time)

    # ustruc
    test = output.field(sid, "uStruc-gv1-ustruc")
    gold = golden.field(sid, "uStruc-gv1-ustruc")
    success = (test == gold).all()
    print("{:s}: {:s} at t={:8.2e}".format("PASS" if success else "FAIL", "ustruc", time))
    if not success: nfail += 1

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
