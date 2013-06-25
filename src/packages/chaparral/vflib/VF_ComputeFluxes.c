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
@(#)    $RCSfile: VF_ComputeFluxes.c,v $
@(#)    $Revision: 1.4 $  $Date: 2006/03/17 04:59:11 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_ComputeFluxes.c,v $
@(#)
@(#)    DESCRIPTION:  Compute heat flux from radiosity solution.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

void
VF_ComputeFluxes(double radq[], double radj[], double sol[], 
                 double *total_flux, int debug)
{
  int    i, j, n, n_l, n_g, row, *index;
  float  *vf;
  double qsum, *areas;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  if (debug>=VF_OUTPUT_VERBOSE && VFLIB_Rank==0) {
        printf("   Computing final fluxes\n");
  }
  n_g = enclosure->npatches_g;
  n_l = enclosure->npatches_l;
  if (radj==NULL) {
    for (i=0; i<n_l; i++) {
      qsum  = 0.0;
      n     = enclosure->row[i].array0.cnt;
      vf    = enclosure->row[i].array0.data;
      index = enclosure->row[i].array0.index;
      for (j=0; j<n; j++) {
   	qsum += vf[j]*sol[index[j]];
      }
      n     = enclosure->row[i].array1.cnt;
      vf    = enclosure->row[i].array1.data;
      index = enclosure->row[i].array1.index;
      for (j=0; j<n; j++) {
   	qsum += vf[j]*sol[index[j]];
      }
      row	= enclosure->row[i].global_index;
      radq[row] = qsum-sol[row];
    }
    VF_ExchangeDouble(radq);
  } else {
    for (i=0; i<n_l; i++) {
      qsum  = 0.0;
      n     = enclosure->row[i].array0.cnt;
      vf    = enclosure->row[i].array0.data;
      index = enclosure->row[i].array0.index;
      for (j=0; j<n; j++) {
  	qsum += vf[j]*sol[index[j]];
      }
      n     = enclosure->row[i].array1.cnt;
      vf    = enclosure->row[i].array1.data;
      index = enclosure->row[i].array1.index;
      for (j=0; j<n; j++) {
  	qsum += vf[j]*sol[index[j]];
      }
      row	= enclosure->row[i].global_index;
      radq[row] = qsum-sol[row];
      radj[row] = sol[row];
    }
    VF_ExchangeDouble(radq);
    VF_ExchangeDouble(radj);
  }
  if (debug>=VF_OUTPUT_VERBOSE) {
    VF_GetDPbuffer1_ptr(&areas);
    VF_GetMatrixAreas(areas);
    for (qsum=0.0, i=0; i<n_g; i++) {
      qsum += areas[i]*radq[i];
    }
    *total_flux = qsum;
  } else {
    *total_flux = 0.0;
  }
}


void
VF_ComputeFluxes_Old(double radq[], double sol[], double *total_flux, int debug)
{
  int    i, j, n, n_l, row, *index;
  float  *vf;
  double qsum;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  if (debug>=VF_OUTPUT_VERBOSE && VFLIB_Rank==0) {
        printf("   Computing final fluxes\n");
  }
  n_l = enclosure->npatches_l;
  for (i=0; i<n_l; i++) {
    qsum  = 0.0;
    n     = enclosure->row[i].array0.cnt;
    vf    = enclosure->row[i].array0.data;
    index = enclosure->row[i].array0.index;
    for (j=0; j<n; j++) {
      qsum += vf[j]*sol[index[j]];
    }
    n     = enclosure->row[i].array1.cnt;
    vf    = enclosure->row[i].array1.data;
    index = enclosure->row[i].array1.index;
    for (j=0; j<n; j++) {
      qsum += vf[j]*sol[index[j]];
    }
    row       = enclosure->row[i].global_index;
    radq[row] = qsum-sol[row];
  }
  if (debug>=VF_OUTPUT_VERBOSE) {
    for (qsum=0.0, i=0; i<n_l; i++) {
      row   = enclosure->row[i].global_index;
      qsum += enclosure->row[i].area*radq[row];
    }
    VF_GlobalSumDouble(&qsum);
    *total_flux = qsum;
  } else {
    *total_flux = 0.0;
  }
  VF_ExchangeDouble(radq);
}
