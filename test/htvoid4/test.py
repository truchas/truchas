#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "htvoid4.inp")
    golden = tenv.output("htvoid4_golden/htvoid4.h5")

    # checking final values
    sid = output.num_series()
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
    nfail += truchas.compare_max(test, gold, 2e-9, "temp", time)

    # solid fraction
    test = output.field(sid, "VOF")[:,1]
    gold = golden.field(sid, "VOF")[:,1]
    nfail += truchas.compare_max(test, gold, 1e-12, "vof", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
