import truchas

def run_test(tenv):
    nfail = 0
    stdout, output = tenv.truchas(4, "time-power-function.inp")
    golden = tenv.output("power-function_golden/power-function.h5")

    sid = output.num_series()
    time = output.time(sid)
    nfail += truchas.compare_max_rel(output.field(sid, "Z_TEMP"),
                                     golden.field(sid, "Z_TEMP"),
                                     1e-7, "temp", time)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0

