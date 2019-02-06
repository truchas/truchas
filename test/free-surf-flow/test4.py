#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "free-surf-flow-4.inp")
    golden = tenv.output("free-surf-flow-4_golden/free-surf-flow-4.h5")

    # initial
    time = output.time(1)

    test = output.field(1, "VOF")[:,0]
    gold = golden.field(1, "VOF")[:,0]
    nfail += truchas.compare_max(test, gold, 0, "vof", time)

    test = output.field(1, "Z_P")
    gold = golden.field(1, "Z_P")
    nfail += truchas.compare_max(test, gold, 1e-13, "pressure", time)

    # final time
    time = output.time(2)

    test = output.field(2, "VOF")[:,0]
    gold = golden.field(2, "VOF")[:,0]
    nfail += truchas.compare_max(test, gold, 1e-13, "vof", time)

    test = output.field(2, "Z_P")
    gold = golden.field(2, "Z_P")
    nfail += truchas.compare_max(test, gold, 1e-13, "pressure", time)

    test = output.field(2, "Z_VC")
    gold = golden.field(2, "Z_VC")
    nfail += truchas.compare_max(test, gold, 1e-13, "velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
