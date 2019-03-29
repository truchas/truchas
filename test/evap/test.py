#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "evap.inp")
    golden = tenv.output("evap_golden/evap.h5")

    for sid in (1, 2, 3):
        test = output.field(sid, "Z_TEMP")
        gold = golden.field(sid, "Z_TEMP")
        nfail += truchas.compare_max_rel(test, gold, 1e-7, "temp", output.time(sid))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
