#!/usr/bin/env python3

import scipy as sp

import truchas


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "free-surf-flow-old-3.inp")

    xc = output.centroids()

    # initial conditions
    vof = output.field(1, "VOF")
    nfail += vof_test(vof, xc, 1e-8, output.time(1))
    nfail += truchas.compare_max(output.field(1, "Z_P"), 0, 1e-11, "pressure", output.time(1))
    nfail += velocity_test(output.field(1, "Z_VC"), vof[:,0], 1e-13, output.time(1))

    # intermediate time
    vof = output.field(2, "VOF")
    nfail += vof_test(vof, xc, 1e-8, output.time(2))
    nfail += truchas.compare_max(output.field(2, "Z_P"), 0, 1e-14, "pressure", output.time(2))
    nfail += velocity_test(output.field(2, "Z_VC"), vof[:,0], 1e-13, output.time(2))

    # final time
    vof = output.field(3, "VOF")
    nfail += vof_test(vof, xc, 1e-8, output.time(3))
    nfail += truchas.compare_max(output.field(3, "Z_P"), 0, 1e-14, "pressure", output.time(3))
    nfail += velocity_test(output.field(3, "Z_VC"), vof[:,0], 1e-13, output.time(3))

    truchas.report_summary(nfail)
    return nfail


def vof_test(vof, xc, tol, time):
    # analytic vof solution at cell centroids
    pl = -2 + time
    pr = -1 + time
    vof_ex = sp.array([0 if x <= pl-0.1 or x >= pr+0.1
                       else 1 if pl+0.1 <= x <= pr-0.1
                       else 5*( x - pl + 0.1) if x < pl+0.1
                       else 5*(-x + pr + 0.1)
                       for x in xc[:,0]])
    return truchas.compare_max(vof[:,0], vof_ex, tol, "VOF", time)


def velocity_test(velocity, vof, tol, time):
    # the x-velocity is 1 in cells containing fluid
    nfail = 0
    uerror = max(abs(u - 1) if vf > 0 else abs(u) for u,vf in zip(velocity[:,0],vof))
    nfail += truchas.compare_max(uerror, 0, tol, "x-velocity", time)
    nfail += truchas.compare_max(velocity[:,1], 0, tol, "y-velocity", time)
    nfail += truchas.compare_max(velocity[:,2], 0, tol, "z-velocity", time)
    return nfail


if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
