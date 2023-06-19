#!/usr/bin/env python3

import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "test.inp")
    golden = tenv.output("test_golden/test.h5")

    for sid in range(1, output.num_series()+1):
      time = output.time(sid)
      
      test = output.field(sid, "Z_TEMP")
      gold = golden.field(sid, "Z_TEMP")
      nfail += truchas.compare_max(test, gold, 1e-6, "temp", time)
      
      test = output.field(sid, "VOF")
      gold = golden.field(sid, "VOF")
      nfail += truchas.compare_max(test, gold, 1e-6, "vof", time)

    truchas.report_summary(nfail)
    return nfail
  

if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
