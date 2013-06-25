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
@(#)    $RCSfile: VF_CalcRead.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcRead.c,v $
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

void VF_CalcRead (int encl, char* filename)
{
  int         tmp;
  double      clock0, clock1;
  VFenclosure *enclosure=VF_GetEnclosure(encl);
    
  VF_InfoRead(encl, filename);
  tmp = enclosure->vf_method;
  enclosure->vf_method = VF_FILE_READ;
  VF_OutputCalcStartBanner();
  clock0 = VF_Clock();
  VF_MatrixRead(encl, filename);
  enclosure->vf_method = VF_FILE_READ;
  clock1 = VF_Clock();
  enclosure->time_vf_calc = clock1-clock0;
  VF_OutputCalcEndBanner();
  enclosure->vf_method = tmp;
}
