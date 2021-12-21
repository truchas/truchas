#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "vfrad1-tet.inp")
    golden = tenv.output("vfrad1-tet_golden/vfrad1-tet.h5")

    # cycle number
    cycle = output.cycle(2)
    cycleg = golden.cycle(2)
    status = "PASS" if cycle == cycleg else "FAIL"
    print("{:s}: matching cycle numbers {:d}".format(status, cycle))
    if cycle != cycleg: nfail += 1

    # temperature
    test = output.field(2, "Z_TEMP")
    gold = golden.field(2, "Z_TEMP")
    nfail += truchas.compare_max_rel(test, gold, 1e-7, "temp", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
