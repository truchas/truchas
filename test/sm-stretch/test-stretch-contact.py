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

    # test initial displacements (there is no change in final displacement)
    t = 0
    sid = 1

    # 1. test compression
    displ1 = output1.field(sid, "Displacement")
    displ2 = output2.field(sid, "Displacement")
    nodex1 = output1.node_coordinates()
    nodex2 = output2.node_coordinates()

    displ2, nodex2 = collapse_links(displ2, nodex2)
    displ2, nodex2 = serialize(displ2, nodex2)
    displ1, nodex1 = serialize(displ1, nodex1)

    nfail += truchas.compare_max(displ1[:,0], displ2[:,0], 1e-10, "compress-displacement", t)

    eps1 = output1.field(sid, "epsilon")
    eps2 = output2.field(sid, "epsilon")
    nfail += truchas.compare_max_rel(eps1[:,0], eps2[:,0], 1e-8, "compress-strain", t)

    # 2. test stretching
    stdout3, output3 = tenv.truchas(4, "stretch-contact.inp")
    displ = output3.field(sid, "Displacement")
    strain = output3.field(sid, "epsilon")

    nfail += truchas.compare_max(strain, 0, 1e-10, "stretch-strain", t)
    nfail += truchas.compare_max(displ[:,1:], 0, 1e-10, "stretch-displacement-yz", t)

    # half the mesh has a displacement of 0, the other half a displacement of 1.
    c0 = np.count_nonzero(abs(displ[:,0] - 0) < 1e-10)
    c1 = np.count_nonzero(abs(displ[:,0] - 1) < 1e-10)
    if c0 == 275 and c1 == 275:
        status = "PASS"
    else:
        status = "FAIL"
        nfail += 1
    print(f"{status}: stretch: ncells with displ_x == 0: {c0} (expected 275), ncells with displ_x == 1: {c1} (expected 275)")

    truchas.report_summary(nfail)
    return nfail


def serialize(field, x):
    """Sort a field based on x coordinates."""
    ind = np.lexsort((x[:,2],x[:,1],x[:,0]))
    return field[ind], x[ind]


def collapse_links(field, nodex):
    """Collapse link-paired nodes into one node,
    taking the first value of the link."""
    duplicate = np.zeros(field.shape[0], dtype=bool)
    for i in range(field.shape[0]):
        x = nodex[i]
        for j in range(i+1,field.shape[0]):
            if np.all(x == nodex[j]): duplicate[j] = True

    field = np.array([f for f,d in zip(field,duplicate) if not d])
    nodex = np.array([f for f,d in zip(nodex,duplicate) if not d])

    return field, nodex


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0
