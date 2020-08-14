#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "ds10.inp")

    xc = output.centroids()

    # test final temp
    nseries = output.num_series()
    T = output.field(nseries, "Z_TEMP")

    error = np.empty(T.shape[0])
    for j in range(error.shape[0]):
        if min(abs(xc[j,:2])) < 1e-3:
            # exclude gap cells from error calc
            error[j] = 0
        else:
            Tref = sum(xc[j,:2])
            if xc[j,0] < 0:
                if xc[j,1] < 0:
                    Tref += 2.5
                else:
                    Tref += 2.7
            else:
                if xc[j,1] < 0:
                    Tref += 3.3
                else:
                    Tref += 3.5
            error[j] = abs(T[j] - Tref) / Tref

    nfail += truchas.compare_max(error, 0, 1e-7, "temperature", output.time(nseries))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
