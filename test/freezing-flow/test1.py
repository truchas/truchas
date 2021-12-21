#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "freezing-flow-1.inp")
    golden = tenv.output("freezing-flow-1_golden/freezing-flow-1.h5")

    # early time
    time = output.time(2)
    nfail += truchas.compare_max(output.field(2, "Z_TEMP"), golden.field(2, "Z_TEMP"),
                                 5.2e-7, "temp", time)
    nfail += truchas.compare_max(output.field(2, "Z_P"), 0, 1e-10, "temp", time)
    nfail += truchas.compare_max(output.field(2, "Z_VC"), 0, 1e-10, "velocity", time)
    nfail += truchas.compare_max(output.field(2, "VOF")[:,0], golden.field(2, "VOF")[:,0],
                                 5e-7, "vof", time)

    # final time
    time = output.time(3)
    nfail += truchas.compare_max(output.field(3, "Z_TEMP"), golden.field(3, "Z_TEMP"),
                                 5e-6, "temp", time)
    nfail += truchas.compare_max(output.field(3, "Z_P"), 0, 1e-10, "temp", time)
    nfail += truchas.compare_max(output.field(3, "Z_VC"), 0, 1e-10, "velocity", time)
    nfail += truchas.compare_max(output.field(3, "VOF")[:,0], golden.field(3, "VOF")[:,0],
                                 1e-10, "vof", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
