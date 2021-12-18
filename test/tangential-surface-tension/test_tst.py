#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "tangential-surface-tension.inp")
    golden = tenv.output("tangential-surface-tension_golden/tangential-surface-tension.h5")

    # cycle numbers
    for sid in (1, 2, 3):
        cycle = output.cycle(sid)
        cycleg = golden.cycle(sid)
        status = "PASS" if cycle == cycleg else "FAIL"
        print("{:s}: matching cycle numbers {:d}".format(status, cycle))
        if cycle != cycleg: nfail += 1

    # verify early fields
    sid = 1
    time = output.time(sid)

    # temperature
    test = output.field(sid, "Z_TEMP")
    gold = golden.field(sid, "Z_TEMP")
    nfail += truchas.compare_max(test, gold, 1e-10, "temp", time)

    # velocity
    test = output.field(sid, "Z_VC")
    gold = golden.field(sid, "Z_VC")
    nfail += truchas.compare_max(test, gold, 1e-10, "velocity", time)

    # verify final fields
    sid = 3
    time = output.time(sid)

    # temperature
    test = output.field(sid, "Z_TEMP")
    gold = golden.field(sid, "Z_TEMP")
    nfail += truchas.compare_max(test, gold, 1e-6, "temp", time)

    # velocity
    test = output.field(sid, "Z_VC")
    gold = golden.field(sid, "Z_VC")
    nfail += truchas.compare_max(test, gold, 1e-6, "velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
