#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "stretch-mixed.inp")
    golden = tenv.output("stretch-mixed_golden/stretch-mixed.h5")

    sid = 1
    nfail += truchas.compare_max(output.field(sid, "Displacement"),
                                 golden.field(sid, "Displacement"),
                                 1e-17, "displacement", 0)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
