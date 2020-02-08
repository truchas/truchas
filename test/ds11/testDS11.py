#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ds11.inp")
    golden = tenv.output("ds11_golden/ds11.h5")

    test_region = output.region(1, 2)
    gold_region = golden.region(1, 2)

    # cycle number
    for sid in (2, 4):
        cycle = output.cycle(sid)
        cycleg = golden.cycle(sid)
        status = "PASS" if cycle == cycleg else "FAIL"
        print("{:s}: matching cycle numbers {:d}".format(status, cycle))
        if cycle != cycleg: nfail += 1

    # fields
    for sid in (3, 4):
        time = output.time(sid)

        test = output.field(sid, "Z_TEMP")[test_region]
        gold = golden.field(sid, "Z_TEMP")[gold_region]
        nfail += truchas.compare_max_rel(test, gold, 1e-6, "temp", time)

        test = output.field(sid, "VOF")[:,2] # comp 2 is fluid
        gold = golden.field(sid, "VOF")[:,2]
        nfail += truchas.compare_max(test, gold, 1e-6, "vof", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
