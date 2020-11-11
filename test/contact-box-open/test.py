#!/usr/bin/env python3

import numpy as np
import numpy.linalg as npla

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "contact-box-open.inp")
    golden = tenv.output("contact-box-open_golden/contact-box-open.h5")
    true_region = output.region(1, 2)
    gap_region_node = output.region_node(3)

    # VERIFY INITIAL
    sid = 1
    time = output.time(sid)

    # stress
    tol = 1e-3
    test = output.field(sid, "sigma")[true_region]
    gold = [-2.288e7, -5.720e6, -1.716e7, 0.0, 0.0, -9.90733e6]
    name = "sigma{:1d}"
    for j in range(6):
        if j == 5: tol = 2e0
        nfail += truchas.compare_max(test[:,j], gold[j], tol, name.format(j+1), time)

    # strain
    tol = 1e-8
    test = output.field(sid, "epsilon")[true_region]
    gold = [0.0, 3.3e-4, 1.1e-4, 0.0, 0.0, -1.90526e-4]
    name = "epsilon{:1d}"
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[j], tol, name.format(j+1), time)

    # normal traction
    tol = 1e-6
    gold = -2.288e7
    test = output.field(sid, "NTRAC_02")[gap_region_node]
    nfail += truchas.compare_max_rel(test, gold, tol, "normal-traction", time)

    # VERIFY JUST BEFORE GAP OPENING
    sid = 2
    time = output.time(sid)

    # stress
    tol = 1e1
    test = output.field(sid, "sigma")[true_region]
    gold = golden.field(sid, "sigma")[true_region]
    name = "sigma{:1d}"
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], tol, name.format(j+1), time)

    # strain
    tol = 2e-8
    test = output.field(sid, "epsilon")[true_region]
    gold = golden.field(sid, "epsilon")[true_region]
    name = "epsilon{:1d}"
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], tol, name.format(j+1), time)

    # displacement magnitude
    abs_tol = 1e-9
    rel_tol = 1e-4
    test = output.field(sid, "Displacement")
    gold = golden.field(sid, "Displacement")
    dtest = npla.norm(test, axis=1)
    dgold = npla.norm(gold, axis=1)
    err = abs(dtest - dgold) / (abs_tol + rel_tol*abs(dgold))
    nfail += truchas.compare_max(err, 0, 1, "displacement", time)

    # temperature
    tol = 1e-5
    test = output.field(sid, "Z_TEMP")[true_region]
    gold = golden.field(sid, "Z_TEMP")[true_region]
    nfail += truchas.compare_max_rel(test, gold, tol, "temperature", time)

    # VERIFY JUST AFTER GAP OPENING
    sid = 3
    time = output.time(sid)

    # stress
    tol = 1e1
    test = output.field(sid, "sigma")[true_region]
    gold = golden.field(sid, "sigma")[true_region]
    name = "sigma{:1d}"
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], tol, name.format(j+1), time)

    # strain
    tol = 2e-8
    test = output.field(sid, "epsilon")[true_region]
    gold = golden.field(sid, "epsilon")[true_region]
    name = "epsilon{:1d}"
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[:,j], tol, name.format(j+1), time)

    # normal traction
    tol = 1e-4
    test = output.field(sid, "NTRAC_02")[gap_region_node]
    nfail += truchas.compare_max(test, 0, tol, "normal-traction", time)

    # VERIFY FINAL
    sid = 4
    time = output.time(sid)

    # stress
    tol = 5e1
    test = output.field(sid, "sigma")[true_region]
    gold = [7.62665e6, 1.90666e6, 5.71999e6, 3.81333e6, 6.60487e6, 3.30244e6]
    name = "sigma{:1d}"
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[j], tol, name.format(j+1), time)

    # strain
    tol = 1e-9
    test = output.field(sid, "epsilon")[true_region]
    gold = [-1.46667e-4, -2.56666e-4, -1.83333e-4, 7.33332e-5, 1.27017e-4, 6.35085e-5]
    name = "epsilon{:1d}"
    for j in range(6):
        nfail += truchas.compare_max(test[:,j], gold[j], tol, name.format(j+1), time)

    # normal traction
    tol = 1e-4
    test = output.field(sid, "NTRAC_02")[gap_region_node]
    nfail += truchas.compare_max(test, 0, tol, "normal-traction", time)

    # temperature
    tol = 1e-8
    test = output.field(sid, "Z_TEMP")[true_region]
    gold = 288
    nfail += truchas.compare_max_rel(test, gold, tol, "temperature", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
