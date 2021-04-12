#!/usr/bin/env python3

# This test sets up a box stretched in the x direction, with a gap down the
# middle (in the YZ plane). It tests compression and stretching. In the case of
# stretching, one half of the box should be displaced with a constant value
# equal to the BC, and the other half should have zero displacement equal to the
# other BC. In the case of compression, the gap is touching, and the
# displacement field should be equivalent to the case where there is no gap at
# all.

import truchas
import numpy as np

def run_test(tenv):
    nfail = 0
    stdout1, output1 = tenv.truchas(4, "compress-x.inp")
    stdout2, output2 = tenv.truchas(4, "compress-contact.inp")
    stdout3, output3 = tenv.truchas(4, "stretch-contact.inp")

    # test initial displacements (there is no change in final displacement)
    t = 0
    sid = 1
    displ1 = output1.field(sid, "Displacement")

    # 1. test compression
    displ2 = output2.field(sid, "Displacement")
    nodex = output2.node_coordinates()
    displ2 = remove_links(displ2, nodex)
    nfail += truchas.compare_max(displ1, displ2, 1e-10, "displacement", t)

    # 2. test stretching
    displ3 = output3.field(sid, "Displacement")
    nfail += 1 # TODO

    truchas.report_summary(nfail)
    return nfail


# def sort_to_mesh(field, xin, xout):
#     """Sort a field on xin coordinates to xout coordinates. This doesn't do any
#     kind of math, it's really for reordering data"""
#     fout = sp.zeros(field.shape)


def remove_links(field, nodex):
    """Collapse link-paired nodes into one node,
    taking the first value of the link."""
    duplicate = np.zeros(field.shape[0], dtype=bool)
    for i in range(field.shape[0]):
        x = nodex[i]
        for j in range(i+1,field.shape[0]):
            if np.all(x == nodex[j]): duplicate[j] = True

    collapsed = np.array([f for f,d in zip(field,duplicate) if not d])
    # collapsed = np.empty((sum(~duplicate),field.shape[1]))
    # i = 0
    # for j in range(field.shape[0]):
    #     if not duplicate[j]:
    #         collapsed[i] = field[j,:]
    #         i += 1

    return collapsed


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
