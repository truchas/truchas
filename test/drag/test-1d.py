#!/usr/bin/env python3

import truchas
import math

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(1, "drag-1d.inp")

    time = output.time(2)
    uex  = 1.5*(1-math.exp(-2*time))

    test = output.field(2, "Z_VC")
    nfail += truchas.compare_max(test[:,0], uex, 4e-5, "x-velocity", time)
    nfail += truchas.compare_max(test[:,1], 0, 1e-12, "y-velocity", time)
    nfail += truchas.compare_max(test[:,2], 0, 1e-12, "z-velocity", time)

    time = output.time(3)
    uex  = 1.5*(1-math.exp(-2*time))

    test = output.field(3, "Z_VC")
    nfail += truchas.compare_max(test[:,0], uex, 1e-5, "x-velocity", time)
    nfail += truchas.compare_max(test[:,1], 0, 1e-12, "y-velocity", time)
    nfail += truchas.compare_max(test[:,2], 0, 1e-12, "z-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
