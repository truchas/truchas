#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "advection-3c.inp")
    gold = tenv.output("advection-3c_pgolden/advection-3c.h5")
    vof_gold = gold.field(2, "VOF")

    # verify volume fractions at final time
    time = output.time(2)
    vof = output.field(2, "VOF")
    nfail += truchas.compare_max(vof, vof_gold, 1e-9, "vof", time)

    truchas.report_summary(nfail)
    return nfail

if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
