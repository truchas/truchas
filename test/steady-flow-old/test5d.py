#!/usr/bin/env python3

import numpy as np

import truchas


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "steady-flow-old-5d.inp")

    # Test VOFs against exact
    xc = output.centroids()

    t = output.time(1)
    vof = output.field(1, "VOF")[:,0]
    nfail += truchas.compare_max(vof, vof_ex(xc, t), 1e-12, "VOF", t)

    t = output.time(2)
    vof = output.field(2, "VOF")[:,0]
    nfail += truchas.compare_max(vof, vof_ex(xc, t), 1e-9, "VOF", t)

    t = output.time(3)
    vof = output.field(3, "VOF")[:,0]
    nfail += truchas.compare_max(vof, vof_ex(xc, t), 6e-3, "VOF", t)

    # Test final velocity against exact
    velocity = output.field(3, "Z_VC")
    velocity[:,0] -= 2.82842712474619
    velocity[:,1] -= 2.82842712474619
    nfail += truchas.compare_max(velocity, 0, 1e-13, "velocity", output.time(3))

    # Test pressure against exact
    nfail += truchas.compare_max(output.field(1, "Z_P"), 0, 1e-10, "pressure", output.time(1))
    nfail += truchas.compare_max(output.field(3, "Z_P"), 0, 1e-10, "pressure", output.time(3))

    truchas.report_summary(nfail)
    return nfail


def vof_ex(xc, t):
    p = -4 + 4*t
    return np.array([0 if p-x < -0.75
                     else (p - x + 0.75)**2 if p-x < -0.25
                     else 0.25 + (p - x + 0.25) if p-x < 0.25
                     else 1 - (p - x - 0.75)**2 if p-x < 0.75
                     else 1
                     for x in (xc[:,0]+xc[:,1])/np.sqrt(2) + (xc[:,1] - xc[:,0])/(2*np.sqrt(2))])


if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
