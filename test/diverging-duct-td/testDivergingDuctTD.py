#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "diverging-duct-td.inp")
    xc = output.centroids()

    # time 0
    time = output.time(3)

    # pressure
    pressure = output.field(3, "Z_P")
    pex = exact_pressure(xc[:,0], time)
    nfail += truchas.compare_max_rel(pressure, pex, 3e-3, "pressure", time)

    # velocity
    velx = output.field(3, "Z_VC")[:,0]
    velxex = exact_velocity(xc[:,0], time)
    nfail += truchas.compare_max_rel(velx, velxex, 1e-3, "x-velocity", time)

    # time 1
    time = output.time(6)

    # velocity
    velx = output.field(6, "Z_VC")[:,0]
    velxex = exact_velocity(xc[:,0], time)
    nfail += truchas.compare_max_rel(velx, velxex, 3e-3, "x-velocity", time)

    # time 2
    time = output.time(11)

    # pressure
    pressure = output.field(11, "Z_P")
    pex = exact_pressure(xc[:,0], time)
    nfail += truchas.compare_max_rel(pressure, pex, 1e-3, "pressure", time)

    # velocity
    velx = output.field(11, "Z_VC")[:,0]
    velxex = exact_velocity(xc[:,0], time)
    nfail += truchas.compare_max_rel(velx, velxex, 1e-3, "x-velocity", time)

    truchas.report_summary(nfail)
    return nfail

def vin(t):
    return 1 if t <= 5 else 0.5 if t >= 15 else 1 - 0.05*(t-5)

def exact_velocity(x, t):
    return vin(t) / (1 + x/40)

def exact_pressure(x, t):
    return 0.68 + 0.5*vin(t)**2 * (0.64 - 1 / (1 + x/40)**2)

if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
