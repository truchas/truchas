#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "pipe-flow-1b.inp")

    xc = output.centroids()
    uex = (2 - (xc[:,1] + xc[:,2])**2) / 4

    time = output.time(2)
    test = output.field(2, "Z_VC")
    nfail += truchas.compare_max(test[:,0], uex, 1.3e-3, "x-velocity", time)
    nfail += truchas.compare_max(test[:,1], 0, 5e-11, "y-velocity", time)
    nfail += truchas.compare_max(test[:,2], 0, 5e-11, "z-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
