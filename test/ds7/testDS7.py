#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ds7.inp")
    golden = tenv.output("ds7_golden/ds7.h5")

    for sid in (2, 3):
        time = output.time(sid)

        # cycle number
        cycle = output.cycle(sid)
        cycleg = golden.cycle(sid)
        status = "PASS" if cycle == cycleg else "FAIL"
        print("{:s}: matching cycle numbers {:d}".format(status, cycle))
        if cycle != cycleg: nfail += 1

        # temperature
        test = output.field(sid, "Z_TEMP")
        gold = golden.field(sid, "Z_TEMP")
        nfail += truchas.compare_max_rel(test, gold, 1e-8, "temp", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
