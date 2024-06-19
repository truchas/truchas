#!/usr/bin/env python3

# This test ensures the normal-direction displacement BC produces rotationally
# invariant results.

import os

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
    displr = np.array([rot @ d for d in displ])
    return displr


def run_truchas(tenv, rotation_angle):
    replacements = {"rotation_angles": "{:f}, {:f}, {:f}".format(*rotation_angle)}
    identifier = truchas.TruchasDatabase.identifier(replacements)
    input_file = os.path.join(tenv._working_dir, f"{identifier}.inp")
    tenv.generate_input_deck(replacements, "template-stretch-nx.inp", input_file)
    stdout, output = tenv.truchas(4, input_file)
    return output


def run_test(tenv):
    rotation_angles = [[0,0,0],
                       [0,0,45],
                       [0,0,60],
                       [0,0,90],
                       [45,45,45],
                       ]
    sid = 2 # final output
    output0 = run_truchas(tenv, rotation_angles[0])
    displ0 = output0.field(sid, "Displacement")
    time = output0.time(sid)

    nfail = 0
    for i, ra in enumerate(rotation_angles[1:]):
        output = run_truchas(tenv, ra)

        displ = output.field(sid, "Displacement")
        displ = rotated(displ, -np.array(ra))
        nfail += truchas.compare_max(displ0, displ, 1e-9, f"displacement-{i+1}", time)

        # TODO: compare other fields

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
