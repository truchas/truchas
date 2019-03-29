#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "static-wall.inp")

    # verify zero velocity
    sid = output.series_id(1)
    test = output.field(sid, "Z_VC")
    nfail += truchas.compare_max(test, 0, 1e-11, "velocity", output.time(sid))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
