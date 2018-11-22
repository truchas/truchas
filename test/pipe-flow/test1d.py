#!/usr/bin/env python3

import scipy as sp

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "pipe-flow-1d.inp")

    xc = output.centroids()
    uex = (2 - (xc[:,0]-xc[:,1])**2) / (4*sp.sqrt(2))

    time = output.time(2)
    test = output.field(2, "Z_VC")
    nfail += truchas.compare_max(test[:,0], uex, 1.3e-3, "x-velocity", time)
    nfail += truchas.compare_max(test[:,1], uex, 1.3e-3, "y-velocity", time)
    nfail += truchas.compare_max(test[:,2], 0, 5e-12, "z-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
