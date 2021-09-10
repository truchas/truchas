#!/usr/bin/env python3

# This test ensures the normal-direction displacement BC produces rotationally
# invariant results.

import numpy as np

import truchas

def rotation_matrix(rotation_angles):
    """This performs the rotations in the opposite order as Truchas."""
    t = [np.pi / 180 * a for a in rotation_angles]
    rotx = np.array([[1, 0, 0],
                     [0, np.cos(t[0]), -np.sin(t[0])],
                     [0, np.sin(t[0]),  np.cos(t[0])]])
    roty = np.array([[np.cos(t[1]), 0, np.sin(t[1])],
                     [0, 1, 0],
                     [-np.sin(t[1]), 0, np.cos(t[1])]])
    rotz = np.array([[np.cos(t[2]), -np.sin(t[2]), 0],
                     [np.sin(t[2]),  np.cos(t[2]), 0],
                     [0, 0, 1]])
    rot = rotx @ roty @ rotz
    return rot

def rotated(displ, rotation_angles):
    rot = rotation_matrix(rotation_angles)
    displr = np.array([ rot @ d for d in displ])
    return displr

def run_test(tenv):
    nfail = 0
    sid = 2 # final output
    stdout0, output0 = tenv.truchas(4, "stretch-nx-0.inp")
    displ0 = output0.field(sid, "Displacement")
    time = output0.time(sid)

    # opposite the rotation angles in the input files.
    rotation_angles = [[0,0,-45],
                       [0,0,-60],
                       [0,0,-90],
                       [-45,-45,-45],
                       ]

    for i, ra in zip(range(1,5), rotation_angles):
        case = f"stretch-nx-{i}"
        stdout, output = tenv.truchas(4, f"{case}.inp")

        displ = output.field(sid, "Displacement")
        displ = rotated(displ, ra)
        nfail += truchas.compare_max(displ0, displ, 1e-9, f"displacement-{i}", time)

        # TODO: compare other fields

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
