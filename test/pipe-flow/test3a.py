#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "pipe-flow-3a.inp")

    flow_region = output.region(1)
    xc = output.centroids()[flow_region]

    # pressure
    pex = 3*(0.5 - xc[:,0])
    nfail += truchas.compare_max(output.field(1, "Z_P")[flow_region], pex, 1e-13, "pressure", output.time(1))
    nfail += truchas.compare_max(output.field(2, "Z_P")[flow_region], pex, 1e-11, "pressure", output.time(2))

    # velocity
    test = output.field(2, "Z_VC")[flow_region]
    uex = (1 - xc[:,1]**2) / 2

    time = output.time(2)
    nfail += truchas.compare_max(test[:,0], uex, 1.3e-3, "x-velocity", time)
    nfail += truchas.compare_max(test[:,1], 0, 1e-13, "y-velocity", time)
    nfail += truchas.compare_max(test[:,2], 0, 1e-13, "z-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
