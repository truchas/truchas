#!/usr/bin/env python3

import scipy as sp
import scipy.linalg as spla

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "hydrostatic-old.inp")

    xc = output.centroids()
    cycle = 20
    sid = output.series_id(cycle)
    time = output.time(sid)

    # pressure
    vof = output.field(sid, "VOF")[:,0]
    pressure = sp.array([p for p, vf in zip(output.field(sid, "Z_P"), vof) if vf > 0.99])
    pressure -= sp.mean(pressure)

    dpdz = -9.81e3
    pex = [dpdz*z for z, vf in zip(xc[:,2], vof) if vf > 0.99]
    pex -= sp.mean(pex)

    error = spla.norm(pressure - pex) / sp.sqrt(pex.size)
    nfail += truchas.compare_l2(error, 0, 1e-5, "pressure", time)

    # velocity zero everywhere
    nfail += truchas.compare_max(output.field(sid, "Z_VC"), 0, 1e-9, "velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
