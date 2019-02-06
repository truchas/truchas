#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "diverging-duct-1b.inp")
    xc = output.centroids()
    time = output.time(2)

    # pressure
    pressure = output.field(2, "Z_P")
    pex = 2.5 - 2 / (xc[:,0]**2 + 0.5*(xc[:,1]+xc[:,2])**2)
    nfail += truchas.compare_max_rel(pressure, pex, 0.17, "pressure", time)

    # velocity
    velocity = output.field(2, "Z_VC")
    uex = xc[:,0] / (xc[:,0]**2 + 0.5*(xc[:,1]+xc[:,2])**2)
    vex = 0.5*(xc[:,1]+xc[:,2]) / (xc[:,0]**2 + 0.5*(xc[:,1]+xc[:,2])**2)
    uerr = (velocity[:,0] - uex) / max(abs(uex))
    verr = (velocity[:,1] - vex) / max(abs(vex))
    werr = (velocity[:,2] - vex) / max(abs(vex))
    nfail += truchas.compare_max(uerr, 0, 0.02, "x-velocity", time)
    nfail += truchas.compare_max(verr, 0, 0.1, "y-velocity", time)
    nfail += truchas.compare_max(werr, 0, 0.1, "z-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
