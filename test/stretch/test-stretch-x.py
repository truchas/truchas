#!/usr/bin/env python3

# This test ensures the normal-direction displacement BC produces the same
# results as a hardwired x/y/z BC in the same direction. It runs both versions
# and checks that they produce the same result.

import truchas

def run_test(tenv):
    nfail = 0
    stdout1, output1 = tenv.truchas(4, "stretch-x.inp")
    stdout2, output2 = tenv.truchas(4, "stretch-nx.inp")

    # test final displacement
    sid = 2
    time = output1.time(sid)
    displ1 = output1.field(sid, "Displacement")
    displ2 = output2.field(sid, "Displacement")
    nfail += truchas.compare_max(displ1, displ2, 0.0, "displacement", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
