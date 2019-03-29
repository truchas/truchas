#!/usr/bin/env python3

import truchas


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "1d-ss-hc-xlinear-wedge.inp")

    xc = output.centroids()
    sid = output.num_series()
    time = output.time(sid)
    Tref = 1 + xc[:,0]

    T = output.field(sid, "Z_TEMP")
    err = abs((T - Tref) / Tref)
    nfail += truchas.compare_max(err, 0, 1e-6, "temperature", time)

    H = output.field(sid, "Z_ENTHALPY")
    err = abs((H/2 - Tref) / Tref)
    nfail += truchas.compare_max(err, 0, 1e-6, "enthalpy", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
