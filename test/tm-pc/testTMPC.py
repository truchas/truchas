#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "tm-pc.inp")
    golden = tenv.output("tm-pc_pgolden/tm-pc.h5")

    # test final fields
    sid = 2
    time = output.time(sid)

    # temperature
    temp = output.field(sid, "Z_TEMP")
    gold = golden.field(sid, "Z_TEMP")
    nfail += truchas.compare_max_rel(temp, gold, 2e-4, "temperature", time)

    # thermal strain against value calculated from temperature
    #   Here we compare with a calculated solution where the material has not seen
    #   any non-isothermal phase change, and a golden solution where it has.
    #   The non-isothermal thermal strain is a simple explicit integration of
    #   a nonlinear equation.
    ncell = 18
    cte = [2.2e-5, 2.1e-5, 2.0e-5]
    T0 = 500
    Ti = 400
    Th = 375
    Tl = 350
    pcstrain = [1.15013007E-03, -1.57897042E-03]

    epstherm = output.field(sid, "epstherm")[:,0]
    gold = golden.field(sid, "epstherm")[:,0]

    epsref = np.array([     (T-T0)*cte[0]               if T >= Ti
                       else pcstrain[0] + (T-Ti)*cte[1] if T >= Th
                       else epsgold                     if T >= Tl
                       else pcstrain[1] + (T-Tl)*cte[2]
                       for T,epsgold in zip(temp, gold)])
    nfail += truchas.compare_max_rel(epstherm, epsref, 1e-8, "thermal strain", time)

    # stress
    test = output.field(sid, "sigma")
    l1 = 5.20e+10
    l2 = 2.60e+10
    XXref = -epstherm * (2*l1 + 2*l2 - 2*l1**2 / (l1 + 2*l2))

    nfail += truchas.compare_max_rel(test[:,0], XXref, 1e-8, "sigma_xx", time)
    nfail += truchas.compare_max_rel(test[:,1], XXref, 1e-8, "sigma_yy", time)
    nfail += truchas.compare_max(test[:,2], 0, 5e-1, "sigma_zz", time)
    nfail += truchas.compare_max(test[:,3], 0, 5e-2, "sigma_xy", time)
    nfail += truchas.compare_max(test[:,4], 0, 5e-2, "sigma_xz", time)
    nfail += truchas.compare_max(test[:,5], 0, 5e-1, "sigma_xx", time)

    # strain
    test = output.field(sid, "epsilon")
    ZZref = epstherm * (1 + 2 * l1/(l1 + 2*l2))

    nfail += truchas.compare_max(test[:,0], 0, 1e-10, "epsilon_xx", time)
    nfail += truchas.compare_max(test[:,1], 0, 1e-10, "epsilon_yy", time)
    nfail += truchas.compare_max_rel(test[:,2], ZZref, 1e-8, "epsilon_zz", time)
    nfail += truchas.compare_max(test[:,3], 0, 1e-10, "epsilon_xy", time)
    nfail += truchas.compare_max(test[:,4], 0, 1e-10, "epsilon_xz", time)
    nfail += truchas.compare_max(test[:,5], 0, 1e-10, "epsilon_xx", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
