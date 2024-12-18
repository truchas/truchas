#!/usr/bin/env python3

import numpy as np

import truchas


def run_test(tenv):
    nfail = 0
    stdout, output_base = tenv.truchas(4, "wg1.inp")
    golden_base = tenv.output("wg1_golden/wg1.h5")
    output = output_base.em_data()
    golden = golden_base.em_data()

    test = output.field("|E|")
    gold = golden.field("|E|")

    nfail += truchas.compare_max(test, gold, 5e-7, "E", 0.0)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
