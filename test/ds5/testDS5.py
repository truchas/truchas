#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ds5.inp")
    golden = tenv.output("ds5_pgolden/ds5.h5")

    for sid in (2, 3):
        # cycle number
        cycle = output.cycle(sid)
        cycleg = golden.cycle(sid)
        status = "PASS" if cycle == cycleg else "FAIL"
        print("{:s}: matching cycle numbers {:d}".format(status, cycle))
        if cycle != cycleg: nfail += 1

        # phis
        time = output.time(sid)

        phi1ex = golden.field(sid, "phi1")
        phi1 = output.field(sid, "phi1")
        phi2 = output.field(sid, "phi2") / 1.5
        phi3 = output.field(sid, "phi3") / 2

        nfail += truchas.compare_max_rel(phi1, phi1ex, 1e-6, "phi1", time)
        nfail += truchas.compare_max_rel(phi2, phi1, 1e-14, "phi2", time)
        nfail += truchas.compare_max_rel(phi3, phi1, 1e-14, "phi3", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
