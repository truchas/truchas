#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ds4.inp")
    golden = tenv.output("ds4_pgolden/ds4.h5")

    # test final cycle number
    cycle = output.cycle(2)
    cycleg = golden.cycle(2)
    status = "PASS" if cycle == cycleg else "FAIL"
    print("{:s}: matching final cycle numbers".format(status))
    if cycle != cycleg: nfail += 1

    # test final temperature
    test = output.field(2, "Z_TEMP")
    gold = golden.field(2, "Z_TEMP")
    nfail += truchas.compare_max_rel(test, gold, 1e-6, "temp", output.time(2))

    # test final temperature
    test = output.field(2, "VOF")[:,2]
    gold = golden.field(2, "VOF")[:,2]
    nfail += truchas.compare_max(test, gold, 1e-6, "vof", output.time(2))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
