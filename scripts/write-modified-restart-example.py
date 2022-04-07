#!/usr/bin/env python3

#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import numpy as np

import truchas


# map then modify Truchas data
#tdata = truchas.TruchasMappedData("preheat_output/preheat.h5", "mesh2.gen", 0.001)
tdata = truchas.TruchasMappedData("advection-2b_pgolden/advection-2b.h5", "mesh2a.gen")
sid = tdata.num_series()

# assign values in certain blocks
tdata.assign_value_block(sid, "VOF", 1, [0,1])
tdata.assign_value_block(sid, "Z_TEMP", 1, 0)

# in block 1, overwrite rho with x
xc = tdata.centroids()
tdata.assign_value_block(sid, "Z_RHO", 1, xc[:,0])

# do something more complicated?
xc = tdata.centroids()
rho = tdata.field(sid, "Z_RHO")
rho = np.where(tdata.blockid() == 1, xc[:,0], rho)
tdata.assign_field(sid, "Z_RHO", rho)

tdata.write_restart("preheat.restart", sid)
