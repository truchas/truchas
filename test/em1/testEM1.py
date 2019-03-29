#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "em1.inp")
    golden = tenv.output("em1_pgolden/em1.h5")

    conducting_region = output.region(11, 12)

    # initial joule heat
    sid = 1

    test = output.field(sid, "Joule_P")[conducting_region]
    gold = golden.field(sid, "Joule_P")[conducting_region]
    nfail += truchas.compare_max_rel(test, gold, 1e-8, "joule heat", output.time(sid))

    # cycles 1 - 5 should be identical to the initial joule heat
    gold = test
    for cycle in range(1,6):
        sid = output.series_id(cycle)
        time = output.time(sid)
        test = output.field(sid, "Joule_P")[conducting_region]
        nfail += truchas.compare_max(test, gold, 0, "joule heat", output.time(sid))

    # # scaled joule heat
    # sid = output.series_id(6)
    # gold = 4 * test
    # test = output.field(sid, "Joule_P")[conducting_region]
    # nfail += truchas.compare_max_rel(test, gold, 1e-15, "joule heat", output.time(sid))

    # # cycles 7 - 11 should be identical to the cycle 6 joule heat
    # gold = test
    # for cycle in range(7,12):
    #     sid = output.series_id(cycle)
    #     test = output.field(sid, "Joule_P")[conducting_region]
    #     nfail += truchas.compare_max(test, gold, 0, "joule heat", output.time(sid))

    # # zero joule heat
    # for cycle in range(11,16):
    #     sid = output.series_id(cycle)
    #     test = output.field(sid, "Joule_P")[conducting_region]
    #     nfail += truchas.compare_max(test, 0, 0, "joule heat", output.time(sid))

    # # final joule heat
    # sid = output.series_id(16)
    # test = output.field(sid, "Joule_P")[conducting_region]
    # gold = golden.field(sid, "Joule_P")[conducting_region]
    # nfail += truchas.compare_max_rel(test, gold, 1e-8, "joule heat", output.time(sid))

    # # cycles 17-20 should be identical to the cycle 16 joule heat
    # gold = test
    # for cycle in range(17,21):
    #     sid = output.series_id(cycle)
    #     test = output.field(sid, "Joule_P")[conducting_region]
    #     nfail += truchas.compare_max(test, gold, 0, "joule heat", output.time(sid))

    # # free space joule heat
    # free_space_region = output.region(10)
    # for cycle in range(0,21):
    #     sid = output.series_id(cycle)
    #     test = output.field(sid, "Joule_P")[free_space_region]
    #     nfail += truchas.compare_max(test, 0, 0, "free-space joule heat", output.time(sid))

    # # final temperature
    # test = output.field(sid, "Z_TEMP")[conducting_region]
    # goldf = golden.field(sid, "Z_TEMP")[conducting_region]
    # nfail += truchas.compare_max_rel(test, goldf, 1e-8, "temp", output.time(sid))

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
