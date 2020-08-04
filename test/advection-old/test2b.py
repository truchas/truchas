#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "advection-old-2b.inp")

    # From a special run to compute vof for final exact configuration
    gold = tenv.output("advection-old-2b_golden/advection-old-2b.h5")
    vof_gold = gold.field(1, "VOF")[:,0]

# verify volume fractions at final time
    sid = output.num_series()
    time = output.time(sid)
    vof = output.field(sid, "VOF")[:,0]

    # max <= l2 <= sqrt(ncell)*max = 31*max
    nfail += truchas.compare_max(vof, vof_gold, 0.09, "vof", time)
    nfail += truchas.compare_l2(vof, vof_gold, 4*0.09, "vof", time)

    truchas.report_summary(nfail)
    return nfail

if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
