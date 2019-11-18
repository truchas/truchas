#!/usr/bin/env python3

import truchas
import os
import numpy

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "turbulence.inp")

    # Empirically determined velocities :( MAC
    # The pressure values are from the linear pressure profile
    # The min/max velocity values are essentially "golden" values
    umin_gold =  632.852
    umax_gold = 2062.587
    pmin_gold = 25000.0
    pmax_gold = 75000.0

    # umin velocity
    # final x-velocity: order is (cycle#,time,vx,vy,vz)
    filename = os.path.join(output.directory,"umin.dat")
    data = numpy.loadtxt(filename)
    time = data[-1, 0]
    umin = data[-1, 1]
    nfail += truchas.compare_max_rel(umin, umin_gold, 1.2e-6, "umin x-velocity", time)

    # umax velocity
    # final x-velocity: order is (cycle#,time,vx,vy,vz)
    filename = os.path.join(output.directory,"umax.dat")
    data = numpy.loadtxt(filename)
    time = data[-1, 0]
    umax = data[-1, 1]
    nfail += truchas.compare_max_rel(umax, umax_gold, 1.2e-6, "umax x-velocity", time)

    # umax pressure
    # final pressure: order is (cycle#,time,p)
    filename = os.path.join(output.directory,"umaxp.dat")
    data = numpy.loadtxt(filename)
    time = data[-1, 0]
    pmax = data[-1, 1]
    nfail += truchas.compare_max_rel(pmax, pmax_gold, 1e-9, "umax pressure", time)

    # pdown pressure
    # final pressure: order is (cycle#,time,p)
    filename = os.path.join(output.directory,"pdown.dat")
    data = numpy.loadtxt(filename)
    time = data[-1, 0]
    pmin = data[-1, 1]
    nfail += truchas.compare_max_rel(pmin, pmin_gold, 1e-9, "pdown pressure", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
