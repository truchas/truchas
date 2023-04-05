#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "test.inp", restart_file="restart.bin")

    sid = output.num_series()
    time = output.time(sid)

    # Ensure the void did not significantly grow
    void0 = output.field(1, "VOF")[:,2]
    void1 = output.field(sid, "VOF")[:,2]
    nfail += truchas.compare_max(void0, void1, 1e-4, "void-vof", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
