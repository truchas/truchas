#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "hydrostatic-9c.inp")

    fluid_region = output.region(2)
    xc = output.centroids()[fluid_region]

    # pressure
    pex = 2*(xc[:,0]+xc[:,2])/np.sqrt(2)
    pex -= np.mean(pex)
    for sid in (1, 2):
        pressure = output.field(sid, "Z_P")[fluid_region]
        pressure -= np.mean(pressure)

        nfail += truchas.compare_max(pressure, pex, 1e-12, "pressure", output.time(sid))

    # velocity zero everywhere
    nfail += truchas.compare_max(output.field(2, "Z_VC"), 0, 1e-13, "velocity", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
