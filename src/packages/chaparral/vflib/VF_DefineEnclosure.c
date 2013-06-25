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
@(#)    $RCSfile: VF_DefineEnclosure.c,v $
@(#)    $Revision: 1.7.2.2 $  $Date: 2006/04/10 22:51:03 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_DefineEnclosure.c,v $
@(#)
@(#)    DESCRIPTION:  Top-level routine to calculate viewfactors.
@(#)
@(#)    enclID        enclosure ID string
@(#)    noblocking    =0, blocking may exist; =1, nonblocking environment
@(#)    partial       =0, full enclosure; =1, partial enclosure
@(#)    npatches      number of patches on local processor
@(#)    global_ids    global ids of the patches on local processor
@(#)    asink         sink area for partial enclosure
@(#)    debug_level   output level for viewfactor calculation
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "vf.h"

static int VF_RandomizeSurfaces=TRUE;

int VF_DefineEnclosure (char* enclID, int nonblocking,
                        int partial, double asink, 
                        int npatches, int global_ids[], 
                        int debug_level)
{
  int          npatches_l, *gids, *gid_map;
  int          i, ii, n, size, offset=0;
  int          blocks, remainder, more;
  int          index = VF_NewEnclosure(enclID,debug_level);
  float        *buffer;
  double       clock0, clock1;
  VFenclosure* encl  = VF_GetEnclosure(index);

  clock0                 = VF_Clock();
  encl->nonblocking      = nonblocking;
  encl->partial          = partial;
  encl->asink            = asink;
  encl->npatches_g       = npatches;
  encl->smoothed         = FALSE;
  encl->sym_method       = VF_SYMMETRIC_NONE;
  encl->seed             = -41731;
  encl->time_vf_init     = 0.0;
  encl->time_vf_calc     = 0.0;
  encl->time_vf_symmetry = 0.0;
  encl->time_vf_smooth   = 0.0;
  ran2(&encl->seed);
  VF_GlobalSumInt(&(encl->npatches_g));
  VF_InitializeBuffers(encl->npatches_g);
  if (debug_level>=VF_OUTPUT_DEBUG_0) {
    VF_PrintSyncStart();
    printf("Processor %d: Input GIDs --------------------------------\n",VFLIB_Rank);
    for (i=0; i<npatches; i++) {
      printf("  global_ids[%d] = %d\n",i,global_ids[i]);
    }
    VF_PrintSyncEnd();
  }
  VF_GetINTbuffer_ptr(&gids);
  VF_GetSPbuffer0_ptr(&buffer);
  gid_map = (int*)buffer;
  if (VFLIB_Size>1) {
    blocks     = (encl->npatches_g-partial)/2;
    remainder  = (encl->npatches_g-partial)%2;
    npatches_l = 2*(int)(blocks/VFLIB_Size);
    more       = blocks%VFLIB_Size;
    if (VFLIB_Rank<more) npatches_l += 2;
    if (VFLIB_Rank==VFLIB_Size-1) {
      if (remainder) npatches_l++;
      if (partial  ) npatches_l++;
    }
    for (size=0, i=0; i<VFLIB_Size; i++) {
      if (i==VFLIB_Rank) {
        n      = npatches;
        offset = size;
        for (ii=0; ii<npatches; ii++) {
          gid_map[offset+ii] = global_ids[ii];
        }
      }
      VF_BroadcastInt(&n, 1, i);
      VF_BroadcastInt(&gid_map[size], n, i);
      size += n;
    }
    for (size=0, i=0; i<VFLIB_Size; i++) {
      if (i==VFLIB_Rank) {
        n      = npatches_l;
        offset = size;
      }
      VF_BroadcastInt(&n, 1, i);
      size += n;
    }
    VF_SortIntegerArray(gid_map, encl->npatches_g);
    if (VF_RandomizeSurfaces) {
      blocks    = (encl->npatches_g-partial)/2;
      remainder = (encl->npatches_g-partial)%2;
      for (i=0; i<blocks; i++) {
        gids[2*i+0] = gid_map[i];
        gids[2*i+1] = gid_map[(encl->npatches_g-partial-1)-i];
      }
      if (remainder) gids[encl->npatches_g-partial-1] = gid_map[blocks];
      if (partial  ) gids[encl->npatches_g-1] = gid_map[encl->npatches_g-1];
      for (i=blocks-1; i>1; i--) {
        ii           = RANI(0,i);
        n            = gids[2*ii];
        gids[2*ii]   = gids[2*i];
        gids[2*i]    = n;
        n            = gids[2*ii+1];
        gids[2*ii+1] = gids[2*i+1];
        gids[2*i+1]  = n;
      }
      VF_SortIntegerArray(&gids[offset], npatches_l);
      encl->seed = -41731;
      ran2(&encl->seed);
    } else {
      memcpy(gids,gid_map,encl->npatches_g*sizeof(int));
    }
  } else {
    offset     = 0;
    npatches_l = npatches;
    memcpy(gids,global_ids,npatches*sizeof(int));
    VF_SortIntegerArray(gids, npatches);
    memcpy(gid_map,gids,npatches*sizeof(int));
  }

  if (!encl->initialized) {
    encl->initialized    = TRUE;
    encl->host2vflib_map = VF_Newi(npatches);
    encl->row            = (VFrow *)VF_Newv(npatches_l*sizeof(VFrow));
    if (VFLIB_Size>1) {
      encl->comm_plan.npartners = VFLIB_Size;
      encl->comm_plan.partners  = (VFcommPartner*)VF_Newv(VFLIB_Size*sizeof(VFcommPartner));
    } else {
      encl->comm_plan.npartners = VFLIB_Size;
      encl->comm_plan.partners  = NULL;
    }
  } else {
    if (VFLIB_Size>1) {
      for (i=0; i<VFLIB_Size; i++) {
        encl->comm_plan.partners[i].proc   = -1;
        encl->comm_plan.partners[i].cnt    = 0;
        encl->comm_plan.partners[i].offset = 0;
        VF_Free(encl->comm_plan.partners[i].global_index);
      }
    }
    for (i=0; i<encl->npatches_l; i++) {
      VF_FreeSparseArray(&(encl->row[i].array0));
      VF_FreeSparseArray(&(encl->row[i].array1));
    }
    if (encl->host_npatches != npatches) {
      encl->host2vflib_map = VF_ReNewi(encl->host2vflib_map, npatches);
    }
    if (encl->npatches_l != npatches_l) {
      encl->row = (VFrow *)VF_ReNewv((void*)encl->row, npatches_l*sizeof(VFrow));
    }
  }
  encl->host_npatches = npatches;
  for (i=0; i<encl->host_npatches; i++) {
    encl->host2vflib_map[i] = VF_LocateIndexInSortedArray(global_ids[i], gid_map, 
                                                          encl->npatches_g);
  }
  encl->npatches_l = npatches_l;
  for (i=0; i<encl->npatches_l; i++) {
    VF_InitializeSparseArray(&(encl->row[i].array0));
    VF_InitializeSparseArray(&(encl->row[i].array1));
  }
  for (i=0; i<encl->npatches_l; i++) {
    encl->row[i].host_gid     = gids[offset+i];
    encl->row[i].global_index = VF_LocateIndexInSortedArray(encl->row[i].host_gid, 
                                                            gid_map, encl->npatches_g);
    encl->row[i].local_index  = i;
  }
  if (VFLIB_Size>1) {
    for (i=0; i<VFLIB_Size; i++) {
      encl->comm_plan.partners[i].proc = i;
      ii = npatches_l;
      VF_BroadcastInt(&ii,1,i);
      encl->comm_plan.partners[i].cnt = ii;
      ii = offset;
      VF_BroadcastInt(&ii,1,i);
      encl->comm_plan.partners[i].offset = ii;
      encl->comm_plan.partners[i].global_index = VF_Newi(encl->comm_plan.partners[i].cnt);
      for (n=0; n<encl->comm_plan.partners[i].cnt; n++) {
        encl->comm_plan.partners[i].global_index[n] = VF_LocateIndexInSortedArray(gids[ii+n], 
                                                            gid_map, encl->npatches_g);
      }
      VF_SortIntegerArray(encl->comm_plan.partners[i].global_index,
                          encl->comm_plan.partners[i].cnt);
    }
    if (debug_level>=VF_OUTPUT_DEBUG_0) {
      VF_PrintSyncStart();
      printf("Processor %d:\n",VFLIB_Rank);
      for (i=0; i<encl->npatches_g; i++) {
        printf("  gid_map[%2d] = %2d\tgids[%2d] = %2d\n",i,gid_map[i],i,gids[i]);
      }
      for (i=0; i<VFLIB_Size; i++) {
        printf("  CommPlan[%d]:\n",i);
        printf("    proc   = %d \n",encl->comm_plan.partners[i].proc);
        printf("    offset = %d\n",encl->comm_plan.partners[i].offset);
        printf("    cnt    = %d\n",encl->comm_plan.partners[i].cnt);
        for (ii=0; ii<encl->comm_plan.partners[i].cnt; ii++) {
          printf("             %d\n",encl->comm_plan.partners[i].global_index[ii]);
        }
      }
      VF_PrintSyncEnd();
    }
  }
  VF_SortMatrixRows();
  for (i=0; i<encl->npatches_l; i++) {
    encl->row[i].local_index = i;
  }
  if (debug_level>=VF_OUTPUT_DEBUG_0) {
    VF_PrintSyncStart();
    printf("Processor %d --------------------------------\n",VFLIB_Rank);
    for (i=0; i<encl->npatches_l; i++) {
      printf("  row %d:\n",i);
      printf("    host_gid     = %d\n",encl->row[i].host_gid);
      printf("    global_index = %d\n",encl->row[i].global_index);
      printf("    local_index  = %d\n",encl->row[i].local_index);
    }
    VF_PrintSyncEnd();
  }
  clock1 = VF_Clock();
  encl->time_vf_init += clock1-clock0;
  return index;
}

void VF_InitializeEnclosure (int index)
{
  int          i;
  double       clock0, clock1;
  VFenclosure* encl = VF_GetEnclosure(index);

  clock0 = VF_Clock();
  ran2(&encl->seed);
  VF_InitializeBuffers(encl->npatches_g);
  if (!encl->initialized) {
    encl->initialized    = TRUE;
    encl->host2vflib_map = VF_Newi(encl->host_npatches);
    encl->row            = (VFrow *)VF_Newv(encl->npatches_l*sizeof(VFrow));
    if (VFLIB_Size>1) {
      encl->comm_plan.npartners = VFLIB_Size;
      encl->comm_plan.partners  = (VFcommPartner*)VF_Newv(VFLIB_Size*sizeof(VFcommPartner));
    } else {
      encl->comm_plan.npartners = VFLIB_Size;
      encl->comm_plan.partners  = NULL;
    }
  }
  for (i=0; i<encl->npatches_l; i++) {
    VF_InitializeSparseArray(&(encl->row[i].array0));
    VF_InitializeSparseArray(&(encl->row[i].array1));
  }
  for (i=0; i<encl->npatches_l; i++) {
    encl->row[i].host_gid     = -1;
    encl->row[i].global_index = -1;
    encl->row[i].local_index  = -1;
  }
  clock1 = VF_Clock();
  encl->time_vf_init += clock1-clock0;
}

void
VF_RandomizeSurfacesOff(void)
{
  VF_RandomizeSurfaces=FALSE;
}

void
VF_RandomizeSurfacesOn(void)
{
  VF_RandomizeSurfaces=TRUE;
}
