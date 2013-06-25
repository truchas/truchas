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
@(#)    $RCSfile: VF_PartialEnclosure.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_PartialEnclosure.c,v $
@(#)
@(#)    DESCRIPTION:  Compute the last row of viewfactors for a
@(#)    partial enclosure.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include "vf.h"

void VF_PartialEnclosure()
{
  int         i, nrows_g, nrows_l;
  float       *rowVF, *colVF, rowsum;
  double      *areas, inva, scale;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();

  VF_Sync();
  /*========================================================*/
  /* COMPLETE THE VIEWFACTOR MATRIX FOR A PARTIAL ENCLOSURE */
  /*========================================================*/
  nrows_g = enclosure->npatches_g;
  nrows_l = enclosure->npatches_l;
  VF_GetSPbuffer0_ptr(&colVF);
  VF_GetSPbuffer1_ptr(&rowVF);
  VF_GetDPbuffer0_ptr(&areas);
  VF_GetMatrixAreas(areas);
  VF_GetMatrixCol(nrows_g-1, VFLIB_Size-1, colVF);
  if (VFLIB_Rank==VFLIB_Size-1) {
    for (rowsum=0.0, i=0; i<nrows_g-1; i++) {
      scale    = areas[i]*(double)topology->nrotations;
      rowVF[i] = (float)(colVF[i]*scale);
      rowsum  += rowVF[i];
    }
    if (enclosure->debug_level>=VF_OUTPUT_SUMMARY) {
      printf("   Minimum Partial Enclosure Area = %g\n\n",rowsum);
    }
#ifdef NO_AREA
    inva = 1.0/(1.01*rowsum);
    for (rowsum=0.0, i=0; i<nrows_g-1; i++) {
      rowVF[i] *= inva;
      rowsum   += rowVF[i];
    }
#else
    if (enclosure->asink<rowsum) {
      printf("   WARNING!!  Partial Enclosure Area must be >= %g\n\n",rowsum);
    }
    inva = (double)topology->nrotations/enclosure->asink;
    for (rowsum=0.0, i=0; i<nrows_g-1; i++) {
      scale    = inva*areas[i];
      rowVF[i] = (float)(colVF[i]*scale);
      rowsum  += rowVF[i];
    }
#endif
    rowVF[nrows_g-1] = (float)MAX(0.0,1.0-rowsum);
    VF_UpdateMatrixRow(nrows_l-1,rowVF);
    VF_SetRawRowsum(nrows_l-1);
  }
}
