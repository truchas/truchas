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
        stdout1, output1 = tenv.truchas(4, f"{name}.inp")
        stdout2, output2 = tenv.truchas(4, f"{name}-legacy.inp")
        sid = 2
        time = output1.time(sid)
        out1 = output1.field(sid, "Displacement")
        out2 = output2.field(sid, "Displacement")
        nfail += truchas.compare_max(out1, out2, 1e-12, f"{name}-displacement", time)
        # TODO: strain
        # TODO: stress
        print()

    # TODO: thermal gradients
    # TODO: shear BCs
    # TODO: normal stress BCs
    # TODO: off-axis rotations

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
