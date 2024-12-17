#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import numpy.linalg as npla


def cell_nodes(cn):
    """Map node IDs from degenerate hex to tet, pyramid, wedge, or hex."""
    cn2 = (   cn[1:5] if cn[0] == cn[1] # TET
        else  cn[:5]  if cn[4] == cn[5] # PYR
        else  cn[:6]  if cn[4] == cn[7] # PYR
        else  cn)                       # HEX
    return cn2


def cell_volume(x):
    """Return the volume for a cell, given its node coordinates."""
    vol = (  tet_volume(x) if 4 == x.shape[0]
        else pyr_volume(x) if 5 == x.shape[0]
        else wed_volume(x) if 6 == x.shape[0]
        else hex_volume(x))
    return vol


def tet_volume(x):
    return npla.det([x[1,:]-x[0,:], x[2,:]-x[0,:], x[3,:]-x[0,:]]) / 6


def pyr_volume(x):
    vol = (tet_volume(x[[0,1,3,4],:])
        +  tet_volume(x[[1,2,0,4],:])
        +  tet_volume(x[[2,3,1,4],:])
        +  tet_volume(x[[3,0,2,4],:])) / 2
    return vol


def wed_volume(x):
    vol = (tet_volume(x[[0,1,2,3],:])
        +  tet_volume(x[[4,3,5,1],:])
        +  tet_volume(x[[1,2,3,5],:])
        +  tet_volume(x[[1,2,0,4],:])
        +  tet_volume(x[[3,5,4,0],:])
        +  tet_volume(x[[0,2,5,4],:])) / 2
    return vol


def hex_volume(x):
    vol = (tet_volume(x[[0,1,3,4],:])
        +  tet_volume(x[[1,2,0,5],:])
        +  tet_volume(x[[2,3,1,6],:])
        +  tet_volume(x[[3,0,2,7],:])
        +  tet_volume(x[[4,7,5,0],:])
        +  tet_volume(x[[5,4,6,1],:])
        +  tet_volume(x[[6,5,7,2],:])
        +  tet_volume(x[[7,6,4,3],:])
        +  tet_volume(x[[0,2,7,5],:])
        +  tet_volume(x[[1,3,4,6],:])) / 2
    return vol
