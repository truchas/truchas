#!/usr/bin/env python3

import numpy as np

import truchas


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "broken-dam.inp")
    gold = tenv.output("broken-dam_golden/broken-dam.h5")
    time = output.time(2)

    # test final cycle number
    cycle = output.cycle(2)
    cycleg = gold.cycle(2)
    status = "PASS" if cycle == cycleg else "FAIL"
    print("{:s} : matching final cycle numbers".format(status))
    if cycle != cycleg: nfail += 1

    # test final fluid volume fraction
    vof = output.field(2, "VOF")[:,0]
    vofg = gold.field(2, "VOF")[:,0]
    nfail += truchas.compare_max(vof, vofg, 5e-9, "vof", time)

    # test final velocity
    vel = output.field(2, "Z_VC")
    velg = gold.field(2, "Z_VC")
    nfail += truchas.compare_max(vel[:,0], velg[:,0], 2e-9, "vel-x", time)
    nfail += truchas.compare_max(vel[:,1], velg[:,1], 3e-9, "vel-y", time)

    truchas.report_summary(nfail)
    return nfail


if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
