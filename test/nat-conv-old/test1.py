#!/usr/bin/env python3

import truchas
import os
import numpy

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "nat-conv-old.inp", restart_file="restart.bin")

    tol = 0.015
    report = "{:s}: max {:s} velocity rel error = {:8.2e} (tol={:8.2e})"

    # horizontal velocity
    filename = os.path.join(output.directory,"VhMax.dat")
    vel = numpy.loadtxt(filename)[-1,1]
    gold = 7.585e-5
    error = abs((vel - gold) / gold)
    print(report.format("FAIL" if error > tol else "PASS", "horizontal", error, tol))
    if error > tol: nfail += 1

    # vertical velocity
    filename = os.path.join(output.directory,"VvMax.dat")
    vel = numpy.loadtxt(filename)[-1,3]
    gold = 7.685e-5
    error = abs((vel - gold) / gold)
    print(report.format("FAIL" if error > tol else "PASS", "vertical", error, tol))
    if error > tol: nfail += 1

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
