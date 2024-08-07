#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "stretch-hex.inp")
    golden = tenv.output("stretch-hex_golden/stretch-hex.h5")

    sid = 1
    nfail += truchas.compare_max(output.field(sid, "Displacement"),
                                 golden.field(sid, "Displacement"),
                                 2e-10, "displacement", 0)

    nfail += truchas.compare_max(output.field(sid, "sigma"),
                                 golden.field(sid, "sigma"),
                                 1.0, "stress", 0)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
