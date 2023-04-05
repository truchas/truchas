#!/usr/bin/env python3

import sys
import numpy as np
import truchas

tdata = truchas.TruchasData(sys.argv[1])
sid = tdata.num_series()

vof = tdata.field(sid, "VOF")
cellid = 40
#solid_vof = 1e-6
solid_vof = 0
void_vof = 0.01
vof[cellid-1,:] = [solid_vof, 1-solid_vof-void_vof, void_vof]
tdata.assign_field(sid, "VOF", vof)

tdata.write_restart(sys.argv[2], sid)
