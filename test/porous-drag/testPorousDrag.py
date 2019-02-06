#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "porous-drag.inp")

    sid = output.series_id(114)
    time = output.time(sid)

    # velocity
    velocity = output.field(sid, "Z_VC")
    uex = 0.5
    nfail += truchas.compare_max_rel(velocity[:,0], uex, 5e-8, "x-velocity", time)
    nfail += truchas.compare_max(velocity[:,1], 0, 5e-9, "y-velocity", time)
    nfail += truchas.compare_max(velocity[:,2], 0, 5e-9, "z-velocity", time)

    # pressure
    pressure = output.field(sid, "Z_P")
    xc = output.centroids()
    pex = 20000 * (1 - xc[:,0] / 20)
    nfail += truchas.compare_max(pressure, pex, 5e-6, "pressure", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
