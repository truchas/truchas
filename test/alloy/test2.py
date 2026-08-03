#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "alloy2.inp")
    golden = tenv.output("alloy2_golden/alloy2.h5")

    for j in range(4):
        test = output.field(j+1, "Z_TEMP")
        gold = golden.field(j+1, "Z_TEMP")
        nfail += truchas.compare_max_rel(test, gold, 1e-6, "temp", output.time(j+1))

        test = output.field(j+1, "VOF")[:,1]
        gold = golden.field(j+1, "VOF")[:,1]
        nfail += truchas.compare_max(test, gold, 2e-6, "vof", output.time(j+1))

    return nfail

if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
