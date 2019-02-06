#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "nat-conv-tet-old.inp", restart_file="restart-tet.bin")
    golden = tenv.output("nat-conv-tet-old_pgolden/nat-conv-tet-old.h5")

    test = output.field(2, "Z_VC")
    gold = golden.field(3, "Z_VC")

    uerror = (test[:,0] - gold[:,0]) / max(abs(gold[:,0]))
    werror = (test[:,2] - gold[:,2]) / max(abs(gold[:,2]))

    nfail += truchas.compare_max(uerror, 0, 8e-5, "x-velocity", output.time(2))
    nfail += truchas.compare_max(werror, 0, 8e-5, "x-velocity", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
