#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "em2.inp")
    golden = tenv.output("em2_pgolden/em2.h5")

    # Three Joule heat calculations, used for the following cycle groups.
    group = [range(0,5), range(5,16), range(16,21)]

    for g in group:
        cycle = g[0]
        sid = output.series_id(cycle)

        test = output.field(sid, "Joule_P")
        gold = golden.field(sid, "Joule_P")
        nfail += truchas.compare_max_rel(test, gold, 1e-8, "joule heat", output.time(sid))

        for cycle in g[1:]:
            sid = output.series_id(cycle)
            test1 = output.field(sid, "Joule_P")
            nfail += truchas.compare_max(test1, test, 0, "joule heat", output.time(sid))

    # final temperature
    sid = output.series_id(20)
    test = output.field(sid, "Z_TEMP")
    gold = golden.field(sid, "Z_TEMP")
    nfail += truchas.compare_max_rel(test, gold, 1e-8, "temperature", output.time(sid))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
