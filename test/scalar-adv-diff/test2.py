#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "adv2.inp")
    golden = tenv.output("adv2_golden/adv2.h5")

    test = output.field(2, "phi1")
    gold = golden.field(2, "phi1")
    nfail += truchas.compare_max_rel(test, gold, 1e-5, "conc", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
