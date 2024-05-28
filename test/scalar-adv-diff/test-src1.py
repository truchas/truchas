#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "src1.inp")

    sid = 2
    t = output.time(sid)
    xc = output.centroids()
    exact1 = 1 + 4*t*xc[:,0]*(1-xc[:,0])
    exact2 = 2 + (t**3)*(xc[:,0]+xc[:,1])
    
    nfail += truchas.compare_max_rel(output.field(sid, "phi1"), exact1, 1e-2, "conc1", t)
    nfail += truchas.compare_max_rel(output.field(sid, "phi2"), exact2, 1e-4, "conc2", t)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
