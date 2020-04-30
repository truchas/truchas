#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ds1.inp")
    golden = tenv.output("ds1_golden/ds1.h5")

    # test final cycle number
    cycle = output.cycle(2)
    cycleg = golden.cycle(2)
    status = "PASS" if cycle == cycleg else "FAIL"
    print("{:s}: matching final cycle numbers".format(status))
    if cycle != cycleg: nfail += 1

    # test final concentration
    cnctr = output.field(2, "phi1")
    exact = golden.field(2, "phi1")
    nfail += truchas.compare_max_rel(cnctr, exact, 2e-7, "conc", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
