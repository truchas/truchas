#!/usr/bin/env python3

# ZJJ, March 2019. I've converted to new test structure and kept old
# comments. I don't understand this test either.

# NNC, Sept 2013.  I've kept comments by DAK about the test for reference.
# I've migrated the test as it was.  I don't understant it at all.

import scipy as sp

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "void-collapse.inp")

    # ZJJ: This velocity test looks for cycle 113, but the test only runs to cycle 30.
    # It was also written with a function name that was overwritten by the later
    # velocity test, so it was never run. As of March 2019 it would fail if renamed,
    # so I've commented it out here to skip as was done indirectly before.
    # # velocity
    # sid = output.series_id(113)
    # time = output.time(sid)

    # test = output.field(sid, "Z_VC")
    # nfail += truchas.compare_max_rel(test[:,0], 0.5, 5e-8, "x-velocity", time)
    # nfail += truchas.compare_max(test[:,1], 0, 5e-10, "y-velocity", time)
    # nfail += truchas.compare_max(test[:,2], 0, 5e-10, "z-velocity", time)

    # pressure at cycle 1
    #   This test seems to check some numerical discretization
    #   transient that I do not understand. DAK 8-24-10
    sid = output.series_id(1)
    time = output.time(sid)

    vof = output.field(sid, "VOF")[:,0]
    test = sp.ma.masked_where(vof < 0.01, output.field(sid, "Z_P")).compressed()
    P_anal = sp.array([65.9934,55.9944,45.9954,35.9964,25.9974,15.9984,5.9994,0.49995])
    nfail += truchas.compare_max(test, P_anal, 1e-3, "pressure", time)

    # pressure at cycle 2
    sid = output.series_id(2)
    time = output.time(sid)
    test = sp.ma.masked_where(vof < 0.01, output.field(sid, "Z_P")).compressed()
    nfail += truchas.compare_max(test, 0, 1e-3, "pressure", time)

    # void volume fraction at cycle 30
    sid = output.series_id(30)
    time = output.time(sid)
    test = output.field(sid, "VOF")[9,1] # void volume fraction in cell 10
    nfail += truchas.compare_max(test, 3e-4, 1e-7, "void fraction", time)

    # x-velocity at cycle 30
    test = output.field(sid, "Z_VC")[:,0]
    nfail += truchas.compare_max(test, 0.9999, 1e-4, "x-velocity", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
