#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "tangential-surface-tension-couette.inp")

    xc = output.centroids()

    # verify final fields
    sid = 2
    time = output.time(sid)

    # temperature
    test = output.field(sid, "Z_TEMP")
    minx = min(xc[:,0]) - 0.5
    Tref = 1 + (xc[:,0] - minx)
    nfail += truchas.compare_max(test, Tref, 1e-8, "temp", time)

    # pressure
    test = output.field(sid, "Z_P")
    nfail += truchas.compare_max(test, 0, 1e-8, "pressure", time)

    # velocity
    test = output.field(sid, "Z_VC")
    dsig_dx = -1
    viscosity = 20
    minz = min(xc[:,2]) - 0.5
    nfail += truchas.compare_max(test[:,0], dsig_dx / viscosity * (xc[:,2] - minz),
                                 1e-7, "x-velocity", time)
    nfail += truchas.compare_max(test[:,1], 0, 1e-10, "y-velocity", time)
    nfail += truchas.compare_max(test[:,2], 0, 1e-10, "z-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
