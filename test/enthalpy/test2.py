#!/usr/bin/env python3

import truchas
import os
import numpy

def run_test(outA, outB):
    nfail = 0

    # test final cycle number
    if outA.cycle(2) == outB.cycle(2):
      status = "PASS"
    else:
      nfail += 1
      status = "FAIL"
    print("{:s}: final cycle numbers match".format(status))

    # test final temperature
    tempA = outA.field(2,"Z_TEMP")
    tempB = outB.field(2,"Z_TEMP")
    nfail += truchas.compare_max_rel(tempA, tempB, 2e-8, "temperature", outA.time(2))

    # test the enthalpy probe data
    filename = os.path.join(outA.directory, "probe1.dat")
    dataA = numpy.loadtxt(filename)
    filename = os.path.join(outB.directory, "probe1.dat")
    dataB = numpy.loadtxt(filename)
    nfail += truchas.compare_max_rel(dataA[:,1], dataB[:,1], 2e-8, 'enthalpy', -1.0)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    stdout, out0 = tenv.truchas(1, "1-phase-1.inp", output_dir="test2_out0")
    stdout, out1 = tenv.truchas(1, "2-phase-nolh-1.inp", output_dir="test2_out1")
    stdout, out2 = tenv.truchas(1, "2-phase-nolh-2.inp", output_dir="test2_out2")
    stdout, out3 = tenv.truchas(1, "2-phase-nolh-3.inp", output_dir="test2_out3")
    stdout, out4 = tenv.truchas(1, "2-phase-nolh-4.inp", output_dir="test2_out4")
    print("Comparing 1-phase-1 and 2-phase-nolh-1")
    nfail = run_test(out0, out1)
    print("Comparing 2-phase-nolh-1 and 2-phase-nolh-2")
    nfail = nfail + run_test(out1, out2)
    print("Comparing 2-phase-nolh-1 and 2-phase-nolh-3")
    nfail = nfail + run_test(out1, out3)
    print("Comparing 2-phase-nolh-1 and 2-phase-nolh-4")
    nfail = nfail + run_test(out1, out4)
    assert nfail == 0
