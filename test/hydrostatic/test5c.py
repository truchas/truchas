#!/usr/bin/env python3

import numpy as np
import numpy.linalg as npla

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "hydrostatic-5c.inp")

    xc = output.centroids()

    # pressure
    pex = np.array([-2*y if y < 0 else -y-0.5 for y in xc[:,1]])
    pex -= np.mean(pex)
    for sid in (1, 2):
        pressure = output.field(sid, "Z_P")
        pressure -= np.mean(pressure)

        error = npla.norm(pressure - pex) / pex.size
        nfail += truchas.compare_l2(error, 0, 8e-3, "pressure", output.time(sid))

    # velocity zero everywhere
    nfail += truchas.compare_max(output.field(2, "Z_VC"), 0, 1e-13, "velocity", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
