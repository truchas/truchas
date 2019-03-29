#!/usr/bin/env python3

import scipy as sp

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "advection-old-2a.inp")

    xc = output.centroids()

    for sid, tol in ((1, 1e-10), (2, 0.043)):
        vof = output.field(sid, "VOF")[:,0]
        time = output.time(sid)

        p = -6 + 6*time
        vofex = sp.array([1 if x < p-0.75 \
                          else 0 if x > p+0.75 \
                          else 15/16 if x < p-0.25 \
                          else 1/16 if x > p+0.25 \
                          else 0.5
                          for x in xc[:,0] + xc[:,1]/2])
        nfail += truchas.compare_max(vof, vofex, tol, "vof", time)

    truchas.report_summary(nfail)
    return nfail

if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
