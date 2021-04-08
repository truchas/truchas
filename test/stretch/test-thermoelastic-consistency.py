#!/usr/bin/env python3

# This runs a series of comparisons between the new solid mechanics and old, in
# thermoelastic (linear-only) configurations. It ensures both produce the same
# results, within some tolerance. No golden data is used.

import truchas


def run_test(tenv):
    nfail = 0

    # stretch in the x direction, compare displacements, strains
    for d in ("x", "y", "z"):
        name = f"stretch-n{d}"
        nfail += test_consistency(name, f"{name}.inp", f"{name}-legacy.inp", 1e-12)

    # TODO: thermal gradients
    # TODO: shear BCs
    # TODO: normal stress BCs

    # off-axis rotations
    nfail += test_consistency("stretch-nx-0", "stretch-nx-0.inp", "stretch-nx-0-legacy.inp", 1e-9)
    nfail += test_consistency("stretch-nx-2", "stretch-nx-2.inp", "stretch-nx-2-legacy.inp", 1e-9)

    truchas.report_summary(nfail)
    return nfail


def test_consistency(name, infile, infile_legacy, tol):
    nfail = 0
    stdout1, output1 = tenv.truchas(4, infile)
    stdout2, output2 = tenv.truchas(4, infile_legacy)
    sid = 2
    time = output1.time(sid)
    out1 = output1.field(sid, "Displacement")
    out2 = output2.field(sid, "Displacement")
    nfail += truchas.compare_max(out1, out2, tol, f"{name}-displacement", time)
    print()
    # TODO: strain
    # TODO: stress
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
