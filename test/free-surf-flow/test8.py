#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "free-surf-flow-8.inp")
    golden = tenv.output("free-surf-flow-8_golden/free-surf-flow-8.h5")

    for sid in (1, 2):
        time = output.time(sid)

        vof = output.field(sid, "VOF")[:,0]
        gold = golden.field(sid, "VOF")[:,0]
        nfail += truchas.compare_max(vof, gold, 1e-13, "vof", time)

        test = output.field(sid, "Z_P")
        gold = golden.field(sid, "Z_P")
        nfail += truchas.compare_max(test, gold, 2e-11, "pressure", time)

        # ensure pressure is 0 in void
        void_error = max(abs(p) for p,vf in zip(test,vof) if vf == 0)
        nfail += truchas.compare_max(void_error, 0, 1e-12, "void-pressure", time)

        test = output.field(sid, "Z_VC")
        gold = golden.field(sid, "Z_VC")
        nfail += truchas.compare_max(test, gold, 1e-11, "velocity", time)

        # the velocity is 0 in purely void cells
        uerror = max(max(abs(u)) for u,vf in zip(test,vof) if vf == 0)
        nfail += truchas.compare_max(uerror, 0, 1e-11, "void-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
