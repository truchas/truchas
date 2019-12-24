#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "freezing-flow-3.inp")
    golden = tenv.output("freezing-flow-3_pgolden/freezing-flow-3.h5")

    # early time
    time = output.time(2)
    vof = output.field(2, "VOF")[:,0]
    nfail += truchas.compare_max(output.field(2, "Z_TEMP"), golden.field(2, "Z_TEMP"),
                                 1e-7, "temp", time)
    nfail += truchas.compare_max(output.field(2, "Z_P"), 0, 1e-10, "pressure", time)
    nfail += compare_velocity(output.field(2, "Z_VC"), vof, 1e-10, time)
    nfail += truchas.compare_max(vof, golden.field(2, "VOF")[:,0], 5e-9, "vof", time)

    # final time
    time = output.time(3)
    vof = output.field(3, "VOF")[:,0]
    nfail += truchas.compare_max(output.field(3, "Z_TEMP"), golden.field(3, "Z_TEMP"),
                                 1e-5, "temp", time)
    nfail += truchas.compare_max(output.field(3, "Z_P"), 0, 1e-10, "pressure", time)
    nfail += compare_velocity(output.field(3, "Z_VC"), vof, 1e-10, time)
    nfail += truchas.compare_max(vof, golden.field(3, "VOF")[:,0], 5e-6, "vof", time)

    truchas.report_summary(nfail)
    return nfail


def compare_velocity(vel, vof, tol, time):
    # the x-direction velocity is 1 in cells which aren't solidified
    uerror = max(abs(u - 1.) if vf < 1. else abs(u) for u,vf in zip(vel[:,0],vof))
    verror = max(abs(vel[:,1]))
    werror = max(abs(vel[:,2]))

    nfail = 0
    nfail += truchas.compare_max(uerror, 0, tol, "x-velocity", time)
    nfail += truchas.compare_max(verror, 0, tol, "y-velocity", time)
    nfail += truchas.compare_max(werror, 0, tol, "z-velocity", time)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
