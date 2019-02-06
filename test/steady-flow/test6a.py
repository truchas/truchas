#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "steady-flow-6a.inp")

    # From a special run to compute vof for final exact configuration
    gold = tenv.output("steady-flow-6a_golden/steady-flow-6a.h5")

    # Test final VOFs against golden output
    vof_output = output.field(2, "VOF")[:,0]
    vof_gold = gold.field(1, "VOF")[:,0]
    nfail += truchas.compare_max(vof_output, vof_gold, 0.09, "VOF", output.time(2))

    # max <= l2 <= sqrt(ncell)*max = 31*max
    nfail += truchas.compare_l2(vof_output, vof_gold, 3*0.09, "VOF", output.time(2))

    # Test final velocity against exact
    velocity = output.field(2, "Z_VC")
    velocity[:,0] -= 4
    velocity[:,1] -= 3
    nfail += truchas.compare_max(velocity, 0, 1e-13, "velocity", output.time(2))

    # Test pressure against exact
    nfail += truchas.compare_max(output.field(1, "Z_P"), 0, 1e-10, "pressure", output.time(1))
    nfail += truchas.compare_max(output.field(2, "Z_P"), 0, 1e-10, "pressure", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
