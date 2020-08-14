#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "simple-gap.inp")
    golden = tenv.output("simple-gap_golden/simple-gap.h5")
    true_region = output.region(1, 2)
    #gap_region_node = output.region_node(3, 4)

    # verify initial
    sid = 1
    time = output.time(sid)

    # stress
    test = output.field(sid, "sigma")[true_region]
    gold = golden.field(sid, "sigma")[true_region]
    for j in range(6):
        nfail += truchas.compare_max_rel(test[:,j], gold[:,j], 1e-6, "sigma{:1d}".format(j+1), time)

    # strain
    test = output.field(sid, "epsilon")[true_region]
    gold = golden.field(sid, "epsilon")[true_region]
    for j in range(6):
        nfail += truchas.compare_max_rel(test[:,j], gold[:,j], 1e-7, "epsilon{:1d}".format(j+1), time)

    # traction sideset 7
    test = output.field(sid, "NTRAC_07")
    gold = np.ma.masked_values(golden.field(sid, "NTRAC_07"), 0)
    error = abs((test-gold)/gold).max()
    nfail += truchas.compare_max(error, 0, 1e-6, "normal traction 7", time)

    # traction sideset 8
    test = output.field(sid, "NTRAC_08")
    gold = np.ma.masked_values(golden.field(sid, "NTRAC_08"), 0)
    error = abs((test-gold)/gold).max()
    nfail += truchas.compare_max(error, 0, 1e-6, "normal traction 8", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
