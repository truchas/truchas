#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "diag-advect.inp")
    golden = tenv.output("diag-advect_pgolden/diag-advect.h5")

    # Test final cycle number
    nseries = output.num_series()
    cycle = output.cycle(nseries)
    cycleg = golden.cycle(nseries)
    status = "PASS" if cycle == cycleg else "FAIL"
    print("{:s} : matching final cycle numbers".format(status))
    if cycle != cycleg: nfail += 1

    # Test final fluid volume fraction
    vof = output.field(2, "VOF")[:,0] # comp 0 is circle
    vofex = golden.field(2, "VOF")[:,0]
    nfail += truchas.compare_max(vof, vofex, 1e-6, "vof", output.time(nseries))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
