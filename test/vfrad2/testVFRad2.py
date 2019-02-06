#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0

    # run with centered and shifted and compare
    stdout, output_c = tenv.truchas(4, "run-centered.inp")
    stdout, output_s = tenv.truchas(4, "run-shifted.inp")

    # temperature
    Tcent = output_c.field(6, "Z_TEMP")
    Tshft = output_s.field(6, "Z_TEMP")
    nfail += truchas.compare_max(Tcent, Tshft, 1e-13, "temp", output_c.time(6))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
