#!/usr/bin/env python3

import truchas
import os
import numpy

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "test.inp")
    filename = os.path.join(output.directory, "probe.dat")
    dataA = numpy.loadtxt(filename)
    filename = os.path.join(tenv._input_dir, "golden.dat")
    dataB = numpy.loadtxt(filename)
    nfail += truchas.compare_max_rel(dataA[:,1], dataB[:,1], 5e-6, 'temperature', -1.0)
    return nfail

if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
