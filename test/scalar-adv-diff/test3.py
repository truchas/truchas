#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "adv3.inp")
    golden = tenv.output("adv3_golden/adv3.h5")

    test = output.field(2, "phi1")
    gold = golden.field(2, "phi1")
    nfail += truchas.compare_max_rel(test, gold, 1e-3, "conc1", output.time(2))

    test = output.field(2, "phi2")
    gold = golden.field(2, "phi2")
    nfail += truchas.compare_max_rel(test, gold, 1e-3, "conc2", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
