#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "sloshing-flow-old.inp")
    golden = tenv.output("sloshing-flow-old_pgolden/sloshing-flow-old.h5")

    sid = 2
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
    nfail += truchas.compare_max(test, gold, 1e-6, "temp", time)

    # fluid volume fraction
    test = output.field(sid, "VOF")[:,0]
    gold = golden.field(sid, "VOF")[:,0]
    nfail += truchas.compare_max(test, gold, 1e-6, "vof", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
