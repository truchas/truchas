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
@(#)    $RCSfile: VF_MakeMatrixSymmetric.c,v $
@(#)    $Revision: 1.4.4.1 $  $Date: 2006/05/10 18:15:34 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_MakeMatrixSymmetric.c,v $
@(#)
@(#)    description:  Force symmetry of the viewfactor matrix sparsity
@(#)                  pattern.  Need to do this because the stochastic
@(#)                  behavior of the hemicube method (jittering and
@(#)                  aliasing) and raycast method can lead to cases
@(#)                  where f(i,j) and f(j,i) are not both zero or both
@(#)                  nonzero.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

void VF_MatrixStats(int*, int*, int*);
void VF_MatrixRowStats(VFsparse_array*, int, int*, int*, int*);

void VF_MakeMatrixSymmetric(int encl, int method, int output)
{
  int    i, j, cnt0, cnt1, cnt2, knt0, knt1, knt2;
  int    nrows_g, nrows_l, icnt, master; 
  float  *colVF, *rowVF;
  double fct1, fct2;
  double clock0, clock1, *areas;
  VFrow  *row;
  VFenclosure *enclosure=VF_GetEnclosure(encl);
  
  if (output>=VF_OUTPUT_SUMMARY) {
    VF_Sync();
    clock0 = VF_Clock();
    if (enclosure->vf_method==VF_PAIRWISE) method = VF_SYMMETRIC_NONE;
    if (VFLIB_Rank==0) {
      switch (method) {
      case VF_SYMMETRIC_SUB:
        printf("     Enforcing reciprocity by subtraction...\n");
        break;
      case VF_SYMMETRIC_ADD:
        printf("     Enforcing reciprocity by addition...\n");
        break;
      case VF_SYMMETRIC_AVG:
        printf("     Enforcing reciprocity by averaging...\n");
        break;
      default:
        printf("     Reciprocity already enforced...\n");
        break;
      }
    }
  }
  VF_GetSPbuffer0_ptr(&colVF);
  VF_GetSPbuffer1_ptr(&rowVF);
  VF_GetDPbuffer0_ptr(&areas);
  nrows_g = enclosure->npatches_g;
  nrows_l = enclosure->npatches_l;
  if (method==VF_SYMMETRIC_NONE) {
    row = enclosure->row;
    for (i=0; i<nrows_l; i++, row++) {
      row->sym_rowsum = row->raw_rowsum;
    }
    if (output>=VF_OUTPUT_VERBOSE) {
      VF_MatrixStats(&cnt0, &cnt1, &cnt2);
      if (VFLIB_Rank==0) {
        printf("       Nonzero lower triangular entries = %d\n",cnt1);
        printf("       Nonzero upper triangular entries = %d\n",cnt2);
      }
      fflush(stdout);
    }
    return;
  }
  if (output>=VF_OUTPUT_VERBOSE) {
    VF_MatrixStats(&cnt0, &cnt1, &cnt2);
  }
  VF_GetMatrixAreas(areas);
  
  row = enclosure->row;
  for (icnt=0, i=0; i<nrows_g; i++) {
    master = -1;
    if (nrows_l>0 && icnt<nrows_l) {
      if (row->global_index==i) master = VFLIB_Rank;
    }
    VF_GlobalMaxInt(&master);
    VF_GetMatrixCol(i, master, colVF);
    if (master==VFLIB_Rank) {
      VF_GetMatrixRow(icnt, &j, rowVF);
      for (j=0; j<i; j++) {
        rowVF[j] = (float)(areas[j]*colVF[j]/areas[i]);
      }
      for (j=i+1; j<nrows_g; j++) {
        if (rowVF[j]==0.0 && colVF[j]==0.0) continue;
        if (rowVF[j]==0.0 || colVF[j]==0.0) {
          switch (method) {
          case VF_SYMMETRIC_SUB:
            if (colVF[j]==0.0) {
              rowVF[j] = 0.0;
            }
            break;
          case VF_SYMMETRIC_ADD:
            if (rowVF[j]==0.0) {
              rowVF[j] = (float)(areas[j]*colVF[j]/areas[i]);
            }
            break;
          case VF_SYMMETRIC_AVG:
             fct1     = areas[i]*areas[i];
             fct2     = areas[j]*areas[j];
             rowVF[j] = (float)((areas[i]*areas[j]*colVF[j]+fct2*rowVF[j])/(fct1+fct2));
             break;
          }
        } else {
          fct1     = areas[i]*areas[i];
          fct2     = areas[j]*areas[j];
          rowVF[j] = (float)((areas[i]*areas[j]*colVF[j]+fct2*rowVF[j])/(fct1+fct2));
        }
      }
      VF_UpdateMatrixRow(icnt,rowVF);
      icnt++;
      row++;
    }
  }
  
  /*===================================================*/
  /* CALCULATE THE NEW ROWSUM AFTER ENFORCING SYMMETRY */
  /*===================================================*/
  for (i=0; i<nrows_l; i++) {
    VF_SetSymmetricRowsum(i);
  }
  if (output>=VF_OUTPUT_SUMMARY) {
    if (output>=VF_OUTPUT_VERBOSE) {
      VF_MatrixStats(&knt0, &knt1, &knt2);
    }
    VF_Sync();
    clock1 = VF_Clock();
    enclosure->time_vf_symmetry = clock1-clock0;
    if (VFLIB_Rank==0) {
      if (output==VF_OUTPUT_SUMMARY) {
        printf("       Elapsed time = %.2f\n",enclosure->time_vf_symmetry);
      } else if (output>=VF_OUTPUT_VERBOSE) {
        printf("       Nonzero lower triangular entries = %d changed to %d\n",cnt1,knt1);
        printf("       Nonzero upper triangular entries = %d changed to %d\n",cnt2,knt2);
        printf("       Elapsed time                     = %.2f\n",enclosure->time_vf_symmetry);
      }
      fflush(stdout);
    }
  }
}

void VF_MatrixStats(int *diag, int *lower, int *upper)
{
  int         i, ii, nrows_l;
  VFenclosure *enclosure=VF_CurrentEnclosure();
  
  *diag   = 0;
  *lower  = 0;
  *upper  = 0;
  nrows_l = enclosure->npatches_l;
  for (i=0; i<nrows_l; i++) {
    ii = enclosure->row[i].global_index;
    VF_MatrixRowStats(&(enclosure->row[i].array0), 
                      ii, diag, lower, upper);
    VF_MatrixRowStats(&(enclosure->row[i].array1), 
                      ii, diag, lower, upper);
  }
  VF_GlobalSumInt(diag);
  VF_GlobalSumInt(lower);
  VF_GlobalSumInt(upper);
}

void VF_MatrixRowStats(VFsparse_array *array, int diag_index,
                       int *diag, int *lower, int *upper)
{
  int i, nonzeros, *index;

  index    = array->index;
  nonzeros = array->cnt;
  for (i=0; i<nonzeros; i++) {
    if (index[i]==diag_index) {
      *diag += 1;
    } else if (index[i]<diag_index) {
      *lower += 1;
    } else {
      *upper += 1;
    }
  }
}
