#!/usr/bin/env python3

import truchas


def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "gap-rad.inp")

    probe_names = ["left end", "right end", "gap left", "gap right"]
    probe_golden = {"left end":  [1.499998474,  1.371828885],
                    "right end": [0.500001526,  0.628171115],
                    "gap left":  [1.498533630,  1.370837833],
                    "gap right": [0.5014663700, 0.629162167]}

    # verify initial probe data
    time = output.probe(probe_names[0], "TEMP")[0, 1] # 0 is first cycle, 1 is time index
    for pname in probe_names:
        probe_data = output.probe(pname, "TEMP")[0,-1] # 0 for first cycle, -1 for temp data
        nfail += truchas.compare_max(probe_data, probe_golden[pname][0], 1e-9, pname, time)

    # verify final probe data
    # get the last cycle
    time = output.probe(probe_names[0], "TEMP")[-1, 1] # -1 is last cycle, 1 is time index
    for pname in probe_names:
        probe_data = output.probe(pname, "TEMP")[-1,-1] # -1 for last cycle, -1 for temp data
        nfail += truchas.compare_max(probe_data, probe_golden[pname][1], 5e-5, pname, time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
