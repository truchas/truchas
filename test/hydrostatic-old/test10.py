#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "hydrostatic-old-10.inp")

    xc = output.centroids()

    # pressure
    pex = (-2 + 0.05 * xc[:,1]) * xc[:,1]
    pex -= np.mean(pex)
    Tex = 1.5 + xc[:,1] / 10
    for sid in (1, 2):
        time = output.time(sid)

        pressure = output.field(sid, "Z_P")
        pressure -= np.mean(pressure)
        nfail += truchas.compare_max(pressure, pex, 1e-10, "pressure", time)

        nfail += truchas.compare_max(output.field(sid, "Z_TEMP"), Tex, 1e-14, "temperature", time)

    # final velocity zero everywhere
    nfail += truchas.compare_max(output.field(2, "Z_VC"), 0, 1e-13, "velocity", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
