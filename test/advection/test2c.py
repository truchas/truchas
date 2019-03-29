#!/usr/bin/env python3

import scipy as sp

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "advection-2c.inp")

    # From a special run to compute vof for final exact configuration
    gold = tenv.output("advection-2c_golden/advection-2c.h5")
    vof_gold = gold.field(1, "VOF")[:,0]

    # verify volume fractions at final time
    sid = output.num_series()
    time = output.time(sid)
    vof = output.field(sid, "VOF")[:,0]

    # max <= l2 <= sqrt(ncell)*max = sqrt(3)*31*max
    nfail += truchas.compare_max(vof, vof_gold, 0.3, "vof", time)
    nfail += truchas.compare_l2(vof, vof_gold, 4*0.3, "vof", time)

    truchas.report_summary(nfail)
    return nfail

if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
