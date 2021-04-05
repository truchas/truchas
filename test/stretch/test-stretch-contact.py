#!/usr/bin/env python3

# This test sets up a box stretched in the x direction, with a gap down the
# middle (in the YZ plane). It tests compression and stretching. In the case of
# stretching, one half of the box should be displaced with a constant value
# equal to the BC, and the other half should have zero displacement equal to the
# other BC. In the case of compression, the gap is touching, and the
# displacement field should be equivalent to the case where there is no gap at
# all.

import truchas

def run_test(tenv):
    nfail = 0
    stdout1, output1 = tenv.truchas(4, "stretch-x.inp")
    stdout2, output2 = tenv.truchas(4, "stretch-contact-compress.inp")
    stdout3, output3 = tenv.truchas(4, "stretch-contact-split.inp")

    # test final displacement
    sid = 2
    displ1 = output1.field(sid, "Displacement")
    displ2 = output2.field(sid, "Displacement")
    nfail += truchas.compare_max(displ1, displ2, 1e-10, "displacement", output.time(sid))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
