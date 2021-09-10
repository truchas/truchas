#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "contact-0.inp")
    golden = tenv.output("contact-0_golden/contact-0.h5")

    sid = 1
    nfail += truchas.compare_max(output.field(sid, "Displacement"),
                                 golden.field(sid, "Displacement"),
                                 1e-10, "displacement", 0)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
