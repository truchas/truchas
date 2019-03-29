#!/usr/bin/env python3

import scipy as sp

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "hydrostatic-old-2c.inp")

    xc = output.centroids()

    # pressure
    pex = -2*xc[:,1] + sp.sqrt(2) * (xc[:,0] + xc[:,2])
    pex -= sp.mean(pex)
    for sid in (1, 2):
        pressure = output.field(sid, "Z_P")
        pressure -= sp.mean(pressure)

        nfail += truchas.compare_max(pressure, pex, 4e-9, "pressure", output.time(sid))

    # velocity zero everywhere
    nfail += truchas.compare_max(output.field(2, "Z_VC"), 0, 1e-13, "velocity", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
