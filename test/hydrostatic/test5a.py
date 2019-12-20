#!/usr/bin/env python3

import scipy as sp
import scipy.linalg as spla

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "hydrostatic-5a.inp")

    xc = output.centroids()

    # initial vof
    vof = output.field(1, "VOF")[:,0]
    vofex = sp.array([0 if y > 1 else 1 if y < 0 else 0.5 for y in xc[:,1]])
    nfail += truchas.compare_max(vof - vofex, 0, 1e-13, "vof", output.time(1))

    # pressure
    pex = sp.array([-2*y if y < 0 else -y-0.5 for y in xc[:,1]])
    pex -= sp.mean(pex)
    for sid in (1, 2):
        pressure = output.field(sid, "Z_P")
        pressure -= sp.mean(pressure)

        error = spla.norm(pressure - pex) / pex.size
        nfail += truchas.compare_l2(error, 0, 8e-3, "pressure", output.time(sid))

    # velocity zero everywhere
    nfail += truchas.compare_max(output.field(2, "Z_VC"), 0, 1e-13, "velocity", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
