#!/usr/bin/env python3

import scipy as sp

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ds9.inp")

    xc = output.centroids()

    # test final temp
    nseries = output.num_series()
    T = output.field(nseries, "Z_TEMP")
    Tref = 9 + 6*xc[:,0]*xc[:,1] - xc[:,0]**2 - xc[:,1]**2
    nfail += truchas.compare_max_rel(T, Tref, 2e-3, "temperature", output.time(nseries))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
