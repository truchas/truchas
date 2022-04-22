#!/usr/bin/env python3

#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import numpy as np

import truchas

# Load truchas data from the preheat run, map it to a new mesh, and scale that mesh.
tdata = truchas.TruchasMappedData("preheat_output/preheat.h5", "mesh2.gen", 0.001)

# Select the last cycle in the preheat output. Further calls will read data from
# this cycle, overwrite data at this cycle, and generate a restart with our
# modifications.
sid = tdata.num_series()

# Set the VOF to [0, 1] (purely the final material in the input) in block 2
tdata.assign_value_block(sid, "VOF", 2, [0, 1])

# Compute the function f = 1 + 2*x at cell centers across the mesh.
xc = tdata.centroids()
f = 1 + 2 * xc[0]

# Reassign the temperature to the function f, but only in blocks 1 and 3.
tdata.assign_value_block(sid, "Z_TEMP", 1, f)
tdata.assign_value_block(sid, "Z_TEMP", 3, f)

# Alternatively, we could read the whole field and modify it before reassigning.
temp = tdata.field(sid, "Z_TEMP")
graphite_region = tdata.region(1, 3)
temp[graphite_region] = f[graphite_region]
tdata.assign_field(sid, "Z_TEMP", temp)

# Alternatively, we could use Numpy functions
temp = tdata.field(sid, "Z_TEMP")
temp = np.where(tdata.blockid() == 1 or tdata.blockid() == 3, f, temp)
tdata.assign_field(sid, "Z_TEMP", temp)

# Write the custom restart file
tdata.write_restart("preheat.restart", sid)
