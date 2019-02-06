#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "diverging-duct-old-0.inp")
    xc = output.centroids()
    sid = output.series_id(274)

    # pressure
    pressure = output.field(sid, "Z_P")
    exact = 1 - 0.5 / (1 + 0.025 * xc[:,0])**2
    print(max(exact), min(exact))
    nfail += truchas.compare_max_rel(pressure, exact, 5e-3, "pressure", output.time(2))

    # velocity
    velx = output.field(sid, "Z_VC")[:,0]
    exact = 1 / (1 + 0.025 * xc[:,0])
    print(max(exact), min(exact))
    nfail += truchas.compare_max_rel(velx, exact, 5e-3, "x-velocity", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
