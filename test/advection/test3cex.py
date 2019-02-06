#!/usr/bin/env python3

import scipy as sp

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "advection-3c-exact.inp")

    # From a special run to compute vof for final exact configuration
    gold = tenv.output("advection-3c-exact_pgolden/advection-3c-exact.h5")
    vof_gold = gold.field(1, "VOF")[:,0]

    # verify volume fractions at final time
    time = output.time(2)
    vof = output.field(2, "VOF")[:,0]

    # max <= l2 <= sqrt(ncell)*max = sqrt(3)*31*max
    err = abs(vof - vof_gold)
    nfail += truchas.compare_max(err, 0, 0.3, "vof", time)
    nfail += truchas.compare_l2(err / sp.sqrt(len(err)), 0, 0.02, "vof", time)

    truchas.report_summary(nfail)
    return nfail

if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
