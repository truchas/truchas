#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "viscoplastic-ring.inp")
    golden = tenv.output("viscoplastic-ring_golden/viscoplastic-ring.h5")

    # initial time
    sid = 1
    time = output.time(sid)

    # stress
    test = output.field(sid, "sigma")
    gold = golden.field(sid, "sigma")
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], 5e-3, f"sigma{j+1:1d}", time)

    # strain
    test = output.field(sid, "epsilon")
    gold = golden.field(sid, "epsilon")
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], 1e-10, f"epsilon{j+1:1d}", time)

    # normal traction
    test = output.field(sid, "Gap Normal Traction")
    gold = golden.field(sid, "Gap Normal Traction")
    test = test[~np.isnan(test)] # remove NaNs in non-gap areas
    gold = gold[~np.isnan(gold)] # remove NaNs in non-gap areas
    nfail += truchas.compare_max(test, gold, 5e-3, "normal traction", time)

    # final time
    sid = output.series_id(5)
    time = output.time(sid)

    # stress
    test = output.field(sid, "sigma")
    gold = golden.field(sid, "sigma")
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], 2e2, f"sigma{j+1:1d}", time)

    # strain
    test = output.field(sid, "epsilon")
    gold = golden.field(sid, "epsilon")
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], 5e-9, f"epsilon{j+1:1d}", time)

    # plastic strain
    test = output.field(sid, "e_plastic")
    gold = golden.field(sid, "e_plastic")
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], 5e-9, f"eps_plas{j+1:1d}", time)

    truchas.report_summary(nfail)
    return nfail

# def eps_eff(eps):
#     return np.sqrt(2/9 * ((eps[:,0]-eps[:,1])**2 + (eps[:,1]-eps[:,2])**2 + (eps[:,2]-eps[:,0])**2)
#                    + 4/3 * (eps[:,3]**2 + eps[:,4]**2 + eps[:,5]**2))

if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
