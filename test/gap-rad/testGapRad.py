#!/usr/bin/env python3

import truchas
import os
import numpy


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "gap-rad.inp")

    probe_names = ["left_end", "right_end", "gap_left", "gap_right"]
    probe_golden = {"left_end":  [1.499998474,  1.371828885],
                    "right_end": [0.500001526,  0.628171115],
                    "gap_left":  [1.498533630,  1.370837833],
                    "gap_right": [0.5014663700, 0.629162167]}

    # verify initial probe data
    for pname in probe_names:
        filename = os.path.join(output.directory,pname + ".dat")
        data = numpy.loadtxt(filename)
        time = data[0,0]
        probe_data = data[0,1]
        nfail += truchas.compare_max(probe_data, probe_golden[pname][0], 1e-9, pname, time)

    # verify final probe data
    # get the last cycle
    for pname in probe_names:
        filename = os.path.join(output.directory,pname + ".dat")
        data = numpy.loadtxt(filename)
        time = data[-1,0]
        probe_data = data[-1,1]
        nfail += truchas.compare_max(probe_data, probe_golden[pname][1], 5e-5, pname, time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
