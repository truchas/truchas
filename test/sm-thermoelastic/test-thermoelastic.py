#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "thermoelastic-cooling.inp")
    golden = tenv.output("thermoelastic-cooling_golden/thermoelastic-cooling.h5")

    # test initial and final displacement
    for sid in (1, output.num_series()):
        time = output.time(sid)
        displt = output.field(sid, "Displacement")
        displg = golden.field(sid, "Displacement")
        nfail += truchas.compare_max(displt, displg, 1e-6, "displacement", output.time(sid))

        sigmat = output.field(sid, "sigma")
        sigmag = golden.field(sid, "sigma")
        nfail += truchas.compare_max(sigmat, sigmag, 5e3, "stress", output.time(sid))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
