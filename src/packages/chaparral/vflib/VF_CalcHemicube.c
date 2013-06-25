/****************************************************************************\
 *                                                                          *
 *   Copyright (c) 1995, 2000, 2005 Sandia Corporation.                           *
 *   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, *
 *   the U.S. Government retains certain rights in this software.           *
 *   For more info, see the README file in the top-level directory.         * 
 *                                                                          *
\****************************************************************************/

/*
@(#)++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
@(#)
@(#)    $RCSfile: VF_CalcHemicube.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcHemicube.c,v $
@(#)
@(#)    DESCRIPTION:  Top-level routine to calculate viewfactors.
@(#)
@(#)    encl          enclosure number
@(#)
@(#)    sub_divide    maximum number of surface subdivisions
@(#)    hc_resolution hemicube resolution
@(#)    hc_min_sep
@(#)
@(#)    debug_level   output level for viewfactor calculation
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>

#include "vf.h"

void VF_CalcHemicube (int encl, int hc_sub_divide, 
                      int hc_resolution, double hc_min_sep)
{
  double      clock0, clock1;
  VFenclosure *enclosure=VF_GetEnclosure(encl);
    
  enclosure->vf_method = VF_HEMICUBE;
  if (hc_resolution%2 != 0) hc_resolution++;
  enclosure->hemicube.sub_divide     = hc_sub_divide;
  enclosure->hemicube.min_separation = hc_min_sep;
  enclosure->hemicube.resolution     = hc_resolution;
  enclosure->hemicube.min_distance   = MAXDOUBLE;
  enclosure->hemicube.imax           = 1;
  enclosure->hemicube.jmax           = 1;
  VF_OutputCalcStartBanner();
  clock0 = VF_Clock();
  VF_MatrixCalcHemicube();
  clock1 = VF_Clock();
  enclosure->time_vf_calc = clock1-clock0;
  VF_OutputCalcEndBanner();
}
