#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "shear-constr-heat.inp")

    # verify fields at initial time
    sid = 1
    time = output.time(sid)

    # stress
    test = output.field(sid, "sigma")
    nfail += truchas.compare_max(test[:,0], 0, 1e-8, "sigma_xx", time)
    nfail += truchas.compare_max(test[:,1], 0, 1e-8, "sigma_yy", time)
    nfail += truchas.compare_max_rel(test[:,2], -152.533333333, 1e-8, "sigma_zz", time)
    nfail += truchas.compare_max(test[:,3], 0, 1e-8, "sigma_xy", time)
    nfail += truchas.compare_max_rel(test[:,4], 500, 1e-8, "sigma_xz", time)
    nfail += truchas.compare_max(test[:,5], 0, 1e-8, "sigma_xx", time)

    # strain
    test = output.field(sid, "epsilon")
    nfail += truchas.compare_max_rel(test[:,0], 2.93333333333e-3, 1e-8, "epsilon_xx", time)
    nfail += truchas.compare_max_rel(test[:,1], 2.93333333333e-3, 1e-8, "epsilon_yy", time)
    nfail += truchas.compare_max(test[:,2], 0, 1e-8, "epsilon_zz", time)
    nfail += truchas.compare_max(test[:,3], 0, 1e-8, "epsilon_xy", time)
    nfail += truchas.compare_max_rel(test[:,4], 500 / 5.2e4, 1e-8, "epsilon_xz", time)
    nfail += truchas.compare_max(test[:,5], 0, 1e-8, "epsilon_xx", time)

    # displacement
    xn = output.node_coordinates()
    xgold = 2.93333333333e-3 * xn[:,0] + (2*500.0 / 5.2e4) * xn[:,2]
    ygold = 2.93333333333e-3 * xn[:,1]

    test = output.field(sid, "Displacement")
    nfail += truchas.compare_max(test[:,0], xgold, 1e-10, "x-displacement", time)
    nfail += truchas.compare_max(test[:,1], ygold, 1e-10, "y-displacement", time)
    nfail += truchas.compare_max(test[:,2], 0, 1e-10, "z-displacement", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
