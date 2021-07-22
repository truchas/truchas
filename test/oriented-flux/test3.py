#!/usr/bin/env python3

# In this test, a laser heats a plate straight-on, and another at an
# angle. The total enthalpy added should be roughly the same in each
# case.

import numpy as np

import truchas


def run_test(tenv):
    nfail = 0
    stdout1, output1 = tenv.truchas(4, "irrad-direct.inp")
    stdout2, output2 = tenv.truchas(4, "irrad-rot.inp")

    total_added_enthalpy1 = sum(output1.field(2, "Z_ENTHALPY")) - sum(output1.field(1, "Z_ENTHALPY"))
    total_added_enthalpy2 = sum(output2.field(2, "Z_ENTHALPY")) - sum(output2.field(1, "Z_ENTHALPY"))
    total_added_enthalpy1 *= np.sin(np.pi/4) # volume scaling

    rel_err = abs(total_added_enthalpy1 - total_added_enthalpy2) / total_added_enthalpy1
    tol = 1e-6
    truchas.report_test("total enthalpy", output1.time(2), rel_err, tol, "rel")
    if rel_err > tol: nfail += 1

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
