#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "inflow-bc-3.inp")

    xc = output.centroids()
    time = output.time(2)

    nfail += truchas.compare_max(output.field(2, "Z_P"), 0, 1e-15, "pressure", time)
    nfail += truchas.compare_max(output.field(2, "Z_TEMP"), 2, 1e-6, "temperature", time)

    vel = output.field(2, "Z_VC")
    nfail += truchas.compare_max(vel[:,0], 0.25, 1e-15, "x-velocity", time)
    nfail += truchas.compare_max(vel[:,1], 0, 1e-15, "y-velocity", time)
    nfail += truchas.compare_max(vel[:,2], 0, 1e-15, "z-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
