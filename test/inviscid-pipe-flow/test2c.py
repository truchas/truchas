#!/usr/bin/env python3

import numpy as np

import truchas


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "inviscid-pipe-flow-2c.inp")

    # mask everything for flow region
    flow_region = output.region(1)

    # analytic expression for exact pressure
    xc = output.centroids()[flow_region]
    pressure_ex = 6*(0.5 - (xc[:,0]+xc[:,2])/np.sqrt(2))

    # pressure
    for i in range(1,4):
        pressure = output.field(i, "Z_P")[flow_region]
        nfail += truchas.compare_max(pressure, pressure_ex, 1e-13, "pressure", output.time(i))

    # velocity
    for i in range(2,4):
        velocity = output.field(i, "Z_VC")[flow_region]
        nfail += velocity_test(velocity, 1e-14, output.time(i))

    truchas.report_summary(nfail)
    return nfail


def velocity_test(velocity, tol, time):
    nfail = 0
    velex = np.sqrt(2)*time
    nfail += truchas.compare_max_rel(velocity[:,0], velex, tol, "x-velocity", time)
    nfail += truchas.compare_max(velocity[:,1], 0, tol, "y-velocity", time)
    nfail += truchas.compare_max_rel(velocity[:,2], velex, tol, "z-velocity", time)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
