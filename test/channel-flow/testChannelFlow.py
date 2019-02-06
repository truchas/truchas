#!/usr/bin/env python3

import scipy as sp

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "channel-flow.inp", restart_file="restart.bin")

    # verify x-velocity
    sid = output.series_id(51031)
    ncelly = 7 # ncellx = 2
    time = output.time(sid)

    # regular ncelly by ncellx grid -- pulling x-vel from first column
    velx = output.field(sid, "Z_VC")[:ncelly,0]
    velx_ex = analytic_solution(ncelly)

    err = (velx - velx_ex) / velx_ex
    nfail += truchas.compare_max(err, 0, 1e-12, "x-velocity", time)

    truchas.report_summary(nfail)
    return nfail


def analytic_solution(ncelly):
    """Analytic solution to balance between viscous stress and pressure
    gradient. Assumes equal size cells along each axis."""
    grad_pressure = (0 - 1000) / 1
    mu = 1
    dy = 1 / ncelly
    rhs = dy**2 * grad_pressure / mu

    velx = sp.empty(ncelly)
    gamma = sp.empty(ncelly)

    beta = -3
    BCoeff = -2
    velx[0] = rhs / beta
    gamma[0] = 1

    # forward elimination loop
    for j in range(1,ncelly):
        if j==ncelly-1: BCoeff = -3
        gamma[j] = 1 / beta
        beta = BCoeff - gamma[j]
        velx[j] = (rhs - velx[j-1]) / beta

    # back substitution loop
    for j in reversed(range(ncelly-1)):
        velx[j] -= gamma[j+1] * velx[j+1]

    return velx


if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
