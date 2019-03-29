#!/usr/bin/env python3

import scipy as sp

import truchas


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "inviscid-pipe-flow-old-1d.inp")

    # early velocity
    time = output.time(2)
    velocity = output.field(2, "Z_VC")
    velex = sp.sqrt(2)*time
    nfail += truchas.compare_max_rel(velocity[:,0], velex, 1e-11, "x-velocity", time)
    nfail += truchas.compare_max_rel(velocity[:,1], velex, 1e-11, "y-velocity", time)
    nfail += truchas.compare_max(velocity[:,2], 0, 1e-12, "z-velocity", time)

    # final velocity
    time = output.time(3)
    velocity = output.field(3, "Z_VC")
    velex = sp.sqrt(2)*time
    nfail += truchas.compare_max(velocity[:,0], velex, 1e-12, "x-velocity", time)
    nfail += truchas.compare_max(velocity[:,1], velex, 1e-12, "y-velocity", time)
    nfail += truchas.compare_max(velocity[:,2], 0, 1e-13, "z-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
