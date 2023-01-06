#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ustruc-gl-restart.inp", restart_file="ustruc-gl-temp.restart.143")
    golden = tenv.output("ustruc-gl-temp_golden/ustruc-gl-temp.h5")

    time = output.time(2) # final time

    # temperature
    test = output.field(2, "Z_TEMP")
    gold = golden.field(3, "Z_TEMP")
    nfail += truchas.compare_max_rel(test, gold, 1e-5, "temp", time)

    # G
    atol = 20
    rtol = 1e-4
    test = np.ma.masked_values(output.field(2, "ustruc-G"), 0)
    gold = np.ma.masked_values(golden.field(3, "ustruc-G"), 0)
    error = abs((test-gold) / (atol + rtol*abs(gold)))
    nfail += truchas.compare_max(error, 0, 1.0, "G", time)

    # L
    atol = 0
    rtol = 3e-2
    test = np.ma.masked_values(output.field(2, "ustruc-L"), 0)
    gold = np.ma.masked_values(golden.field(3, "ustruc-L"), 0)
    error = abs((test-gold) / (atol + rtol* gold))
    nfail += truchas.compare_max(error, 0, 1.0, "L", time)

    # Solidification time
    atol = 0
    rtol = 2e-4
    test = np.ma.masked_values(output.field(2, "ustruc-t_sol"), 0)
    gold = np.ma.masked_values(golden.field(3, "ustruc-t_sol"), 0)
    print("%d" % test.size)
    print("%d" % gold.size)
    error = abs((test-gold) / (atol + rtol*gold))
    nfail += truchas.compare_max(error, 0, 1.0, "t_sol", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
