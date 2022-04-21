#!/usr/bin/env python3

import numpy as np

import truchas

# Run the 1D problem, then run again from a restart, and compare the final-time output.
def run_test(tenv):
    nfail = 0
    stdout1, output1 = tenv.truchas(4, "viscoplastic-1d.inp", output_dir="complete")
    restart_filename = tenv.working_directory() + "vp1d.restart"
    output1.write_restart(restart_filename, 4)
    stdout2, output2 = tenv.truchas(4, "viscoplastic-1d.inp", restart_file=restart_filename,
                                    output_dir="restarted")

    sid1 = output1.num_series()
    sid2 = output2.num_series()
    time = output1.time(sid1)
    for fieldname in ("sigma", "epsilon", "epsdot", "e_plastic"):
        nfail += truchas.compare_max_rel(output1.field(sid1, fieldname),
                                         output2.field(sid2, fieldname),
                                         1e-9, fieldname, time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
