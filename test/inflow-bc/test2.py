#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "inflow-bc-2.inp")

    xc = output.centroids()
    time = output.time(2)

    # pressure
    pressure = output.field(2, "Z_P")
    nfail += truchas.compare_max(pressure, 0, 1e-15, "pressure", time)

    # vof
    vof = output.field(2, "VOF")[:,0]
    p = 0.25 * time - 0.375
    vofex = np.array([1 if x + 0.125 < p
                      else 0 if x - 0.125 > p
                      else 4*(p-(x-0.125))
                      for x in xc[:,0]])
    nfail += truchas.compare_max(vof, vofex, 1e-15, "vof", time)

    # the x-velocity is 0.25 in cells containing fluid
    vel = output.field(2, "Z_VC")
    velex = np.array([[0.25 if vf > 0 else 0, 0, 0] for vf in vof])
    nfail += truchas.compare_max(vel, velex, 1e-15, "velocity", time)

    # temperature is 2 in cells containing fluid
    temp = output.field(2, "Z_TEMP")
    tempex = np.array([2 if vf > 0 else 0 for vf in vof])
    nfail += truchas.compare_max(temp, tempex, 1e-15, "temperature", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
