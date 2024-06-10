#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "adv1.inp")
    golden = tenv.output("adv1_golden/adv1.h5")

    for i in range(2,6):
        test = output.field(i, "phi1")
        gold = golden.field(i, "phi1")
        nfail += truchas.compare_max_rel(test, gold, 2e-3, "conc", output.time(i))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
