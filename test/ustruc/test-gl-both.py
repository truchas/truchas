#!/usr/bin/env python3

import numpy as np

import truchas

def run_test(tenv):
    nfail = 0
    stdout1, output1 = tenv.truchas(4, "ustruc-gl-both.inp", output_dir="complete")
    
    sid = 2
    time = output1.time(sid)
    
    # Compare some of the results with golden results from those generated with
    # individual simulations. This is is just a quick sanity check. These are
    # the same tests as in test-gl-temp.py and test-gl-frac.py.
    golden = tenv.output("ustruc-gl-temp_golden/ustruc-gl-temp.h5")
    gold = np.ma.masked_values(golden.field(sid, "ustruc1-L"), 0)

    # Fraction-based L (> 0); see test-gl-frac.py
    atol = 0
    rtol = 6e-2
    test = np.ma.masked_values(output1.field(sid, "ustruc2-L"), 0)
    error = abs((test-gold) / (atol + rtol* gold))
    nfail += truchas.compare_max(error, 0, 1.0, "frac L", time)

    # Temperature-based L (> 0); see test-gl-temp.py
    atol = 0
    rtol = 3e-2
    test = np.ma.masked_values(output1.field(sid, "ustruc1-L"), 0)
    error = abs((test-gold) / (atol + rtol* gold))
    nfail += truchas.compare_max(error, 0, 1.0, "temp L", time)

    # Now restart the problem from an intermediate time and compare with the
    # original output.
    restart_filename = tenv.working_directory() + "both.restart"
    output1.write_restart(restart_filename, 1)
    stdout2, output2 = tenv.truchas(4, "ustruc-gl-both.inp", restart_file=restart_filename,
                                    output_dir="restarted")

    sid1 = output1.num_series()
    sid2 = output2.num_series()
    time = output1.time(sid1)
    for fieldname in ("ustruc1-G", "ustruc1-L", "ustruc1-t_sol",
                      "ustruc2-G", "ustruc2-L", "ustruc2-t_sol"):
        nfail += truchas.compare_max_rel(np.ma.masked_values(output1.field(sid1, fieldname), 0),
                                         np.ma.masked_values(output2.field(sid2, fieldname), 0),
                                         1e-9, fieldname, time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
