#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout1, output1 = tenv.truchas(4, "flux.inp")
    stdout2, output2 = tenv.truchas(4, "vflux.inp")

    sid = output1.num_series()
    time = output1.time(sid)
    nfail += truchas.compare_max_rel(output1.field(sid, "Z_TEMP"),
                                     output2.field(sid, "Z_TEMP"),
                                     1e-10, "temp", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
