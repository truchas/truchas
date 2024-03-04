#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "em3.inp")
    golden = tenv.output("em3_golden/em3.h5")

    # joule heat
    test = output.field(1, "Joule_P")
    gold = golden.field(1, "Joule_P")
    nfail += truchas.compare_max_rel(test, gold, 1e-4, "joule heat", output.time(1))

    # final temperature
    sid = output.series_id(30)
    test = output.field(sid, "Z_TEMP")
    gold = golden.field(sid, "Z_TEMP")
    nfail += truchas.compare_max_rel(test, gold, 1e-5, "temperature", output.time(sid))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
