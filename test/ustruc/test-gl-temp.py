#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ustruc-gl-temp.inp")
    golden = tenv.output("ustruc-gl-temp_golden/ustruc-gl-temp.h5")

    sid = 2
    time = output.time(sid)

    # Temperature
    test = output.field(sid, "Z_TEMP")
    gold = golden.field(sid, "Z_TEMP")
    nfail += truchas.compare_max_rel(test, gold, 1e-4, "temp", time)

    # G
    atol = 20
    rtol = 1e-4
    test = np.ma.masked_values(output.field(sid, "ustruc-G"), 0)
    gold = np.ma.masked_values(golden.field(sid, "ustruc-G"), 0)
    error = abs((test-gold) / (atol + rtol*abs(gold)))
    nfail += truchas.compare_max(error, 0, 1.0, "G", time)

    # L (> 0)
    atol = 0
    rtol = 3e-2
    test = np.ma.masked_values(output.field(sid, "ustruc-L"), 0)
    gold = np.ma.masked_values(golden.field(sid, "ustruc-L"), 0)
    error = abs((test-gold) / (atol + rtol* gold))
    nfail += truchas.compare_max(error, 0, 1.0, "L", time)

    # Solidification time (> 0)
    atol = 0
    rtol = 1e-4
    test = np.ma.masked_values(output.field(sid, "ustruc-t_sol"), 0)
    gold = np.ma.masked_values(golden.field(sid, "ustruc-t_sol"), 0)
    error = abs((test-gold) / (atol + rtol* gold))
    nfail += truchas.compare_max(error, 0, 1.0, "t_sol", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
