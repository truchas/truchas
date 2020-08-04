#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "viscoplastic-ring.inp")
    golden = tenv.output("viscoplastic-ring_pgolden/viscoplastic-ring.h5")
    true_region = output.region(1, 2)
    gap_region_node = output.region_node(3)

    # initial time
    sid = 1
    time = output.time(sid)

    # stress
    test = output.field(sid, "sigma")[true_region]
    gold = golden.field(sid, "sigma")[true_region]
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], 0.1, "sigma{:1d}".format(j+1), time)

    # strain
    test = output.field(sid, "epsilon")[true_region]
    gold = golden.field(sid, "epsilon")[true_region]
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], 1e-10, "epsilon{:1d}".format(j+1), time)

    # normal traction
    test = output.field(sid, "NTRAC_04")[gap_region_node]
    gold = golden.field(sid, "NTRAC_04")[gap_region_node]
    nfail += truchas.compare_max(test, gold, 0.1, "normal traction", time)

    # final time
    sid = output.series_id(5)
    time = output.time(sid)

    # stress (loose tolerance)
    test = output.field(sid, "sigma")[true_region]
    gold = golden.field(sid, "sigma")[true_region]
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], 5e3, "sigma{:1d}".format(j+1), time)

    # strain
    test = output.field(sid, "epsilon")[true_region]
    gold = golden.field(sid, "epsilon")[true_region]
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], 5e-9, "epsilon{:1d}".format(j+1), time)

    # plastic strain
    test = eps_eff(output.field(sid, "e_plastic")[true_region])
    gold = eps_eff(golden.field(sid, "e_plastic")[true_region])
    nfail += truchas.compare_max(test, gold, 5e-9, "eps_plas", time)

    truchas.report_summary(nfail)
    return nfail

def eps_eff(eps):
    return np.sqrt(2/9 * ((eps[:,0]-eps[:,1])**2 + (eps[:,1]-eps[:,2])**2 +
                           (eps[:,2]-eps[:,0])**2) +
                    4/3 * (eps[:,3]**2 + eps[:,4]**2 + eps[:,5]**2))

if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
