#!/usr/bin/env python3

# This runs a series of comparisons between the new solid mechanics and old, in
# thermoelastic (linear-only) configurations. It ensures both produce the same
# results, within some tolerance. No golden data is used.

import truchas


def run_test(tenv):
    nfail = 0

    # thermal gradients
    for i in range(3):
        name = f"thermoelastic{i}"
        nfail += test_consistency(1, name, f"{name}.inp", f"{name}-legacy.inp", 1e-10)

    for i in range(3,5):
        name = f"thermoelastic{i}"
        nfail += test_consistency(4, name, f"{name}.inp", f"{name}-legacy.inp", 1e-10)

    # TODO: off-axis rotations

    truchas.report_summary(nfail)
    return nfail


def test_consistency(nproc, name, infile, infile_legacy, tol):
    nfail = 0
    stdout1, output1 = tenv.truchas(nproc, infile)
    stdout2, output2 = tenv.truchas(nproc, infile_legacy)
    sid = 1
    time = output1.time(sid)
    out1 = output1.field(sid, "Displacement")
    out2 = output2.field(sid, "Displacement")
    nfail += truchas.compare_max(out1, out2, tol, f"{name}-displacement", time)
    # TODO: strain
    # TODO: stress
    print()
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
