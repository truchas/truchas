#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "viscoplastic-1d.inp")

    strain_golden = np.array([1e-3, 0, 0, 0, 0, 0])
    vmstrain_golden = 1e-3
    stress_golden = np.array([4e6, 2e6, 2e6, 0, 0, 0])
    vmstress_golden = 2e6
    strain_rate_golden = np.sqrt(2) * 1e-9
    strain_plastic_rate_golden = np.array([strain_rate_golden,
                                           -strain_rate_golden/2,
                                           -strain_rate_golden/2,
                                           0, 0, 0])
    vmpstrain_golden = np.sqrt(4.5) * 1e-9

    # The loose-ish tolerances on stress/strain are so the tiny influence
    # of increasing plastic strain on the strain and stress fields is neglected.
    # Using the Von Mises transformation makes this problem rotationally
    # invariant.
    for sid in range(1,output.num_series()+1):
        time = output.time(sid)

        strain_plastic_golden = time * strain_plastic_rate_golden

        vmstress = von_mises(output.field(sid, "sigma"))
        nfail += truchas.compare_max_rel(vmstress, vmstress_golden, 1e-5, "stress", time)

        vmstrain = von_mises(output.field(sid, "epsilon"))
        nfail += truchas.compare_max_rel(vmstrain, vmstrain_golden, 1e-5, "strain", time)

        epsdot = output.field(sid, "epsdot")
        nfail += truchas.compare_max(epsdot, strain_rate_golden, 1e-12, "epsdot", time)

        vmpstrain = output.field(sid, "e_plastic")
        nfail += truchas.compare_max(vmpstrain, vmpstrain_golden, 1e-8, "plastic strain", time)
        print()

    truchas.report_summary(nfail)
    return nfail

def von_mises(eps):
    return np.sqrt(((eps[:,0]-eps[:,1])**2 + (eps[:,1]-eps[:,2])**2 + (eps[:,2]-eps[:,0])**2
                    + 6 * (eps[:,3]**2 + eps[:,4]**2 + eps[:,5]**2))/2)

if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
