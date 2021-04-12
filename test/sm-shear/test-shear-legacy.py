#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "shear-legacy.inp")
    golden = tenv.output("shear-legacy_golden/shear-legacy.h5")

    # test final displacement
    sid = 2
    displ = output.field(sid, "Displacement")
    gold = golden.field(sid, "Displacement")
    nfail += truchas.compare_max(displ, gold, 1e-10, "displacement", output.time(sid))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
