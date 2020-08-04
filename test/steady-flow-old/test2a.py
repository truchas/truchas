#!/usr/bin/env python3

import numpy as np

import truchas


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "steady-flow-old-2a.inp")

    # Test VOFs against exact
    xc = output.centroids()

    t = output.time(2)
    vof = output.field(2, "VOF")[:,0]
    nfail += truchas.compare_max(vof, vof_ex(xc, t), 2e-13, "VOF", t)

    t = output.time(3)
    vof = output.field(3, "VOF")[:,0]
    nfail += truchas.compare_max(vof, vof_ex(xc, t), 2e-13, "VOF", t)

    # Test final velocity against exact
    velocity = output.field(3, "Z_VC")
    velocity[:,0] -= 4
    nfail += truchas.compare_max(velocity, 0, 1e-13, "velocity", output.time(3))

    # Test pressure against exact
    nfail += truchas.compare_max(output.field(1, "Z_P"), 0, 1e-10, "pressure", output.time(1))
    nfail += truchas.compare_max(output.field(3, "Z_P"), 0, 1e-10, "pressure", output.time(3))

    truchas.report_summary(nfail)
    return nfail


def vof_ex(xc, t):
    p = -4 + 4*t
    return np.array([1 if x < p-0.5 else 0 if x > p+0.5 else p-(x-0.5) for x in xc[:,0]])


if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
