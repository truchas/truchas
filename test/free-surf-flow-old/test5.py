#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "free-surf-flow-old-5.inp")
    golden = tenv.output("free-surf-flow-old-5_golden/free-surf-flow-old-5.h5")

    # final time
    time = output.time(2)

    vof = output.field(2, "VOF")[:,0]
    gold = golden.field(2, "VOF")[:,0]
    nfail += truchas.compare_max(vof, gold, 4e-2, "vof", time)

    test = output.field(2, "Z_P")
    gold = golden.field(2, "Z_P")
    nfail += truchas.compare_max(test, gold, 1e-10, "pressure", time)

    # the x-velocity is 1 in cells containing fluid
    test = output.field(2, "Z_VC")
    uerror = max(abs(u - 1.) if vf > 0.0 else abs(u) for u,vf in zip(test[:,0],vof))
    nfail += truchas.compare_max(uerror, 0, 1e-11, "x-velocity", time)
    nfail += truchas.compare_max(test[:,1], 0, 1e-11, "y-velocity", time)
    nfail += truchas.compare_max(test[:,2], 0, 1e-11, "z-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
