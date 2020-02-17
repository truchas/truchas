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
@(#)    $RCSfile: VF_EnclosureUtils.c,v $
@(#)    $Revision: 1.6.4.1 $  $Date: 2006/04/03 17:57:28 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_EnclosureUtils.c,v $
@(#)
@(#)    DESCRIPTION:  Misc. parallel utility functions.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "vf.h"

static int          num_enclosures=0;
static int          current_enclosure=-1;
static VFenclosure* VF_DATA_Enclosure=NULL;

int
VF_NumEnclosures(void)
{
  return num_enclosures;
}

void
VF_InitializeEnclosures(int nenclosures)
{
  int i, nbytes;

  /* Original function was a no-op if nenenclosures <= num_enclosures;
   * function now always initializes the VF_DATA_Enclosure structure. */
  VF_DeleteEnclosures();
  num_enclosures    = nenclosures;
  nbytes            = sizeof(VFenclosure) * nenclosures;
  VF_DATA_Enclosure = (VFenclosure *)VF_Newv(nbytes);
  for (i=0; i<num_enclosures; i++) {
    VF_DATA_Enclosure[i].id                  = NULL;
    VF_DATA_Enclosure[i].initialized         = FALSE;
    VF_DATA_Enclosure[i].npatches_g          = 0;
    VF_DATA_Enclosure[i].npatches_l          = 0;
    VF_DATA_Enclosure[i].smoothed            = 0;
    VF_DATA_Enclosure[i].sym_method          = VF_SYMMETRIC_NONE;
    VF_DATA_Enclosure[i].row                 = NULL;
#ifdef AZTEC
    VF_DATA_Enclosure[i].data_org            = NULL;
    VF_DATA_Enclosure[i].update_index        = NULL;
    VF_DATA_Enclosure[i].external            = NULL;
    VF_DATA_Enclosure[i].external_index      = NULL;
#endif
    VF_DATA_Enclosure[i].topology            = NULL;
    VF_DATA_Enclosure[i].comm_plan.npartners = 0;
    VF_DATA_Enclosure[i].comm_plan.partners  = NULL;
  }
}

void
VF_DeleteEnclosures(void)
{
  int i, j;
  
  for (i=0; i<num_enclosures; i++) {
    VF_Free(VF_DATA_Enclosure[i].id);
    VF_Free(VF_DATA_Enclosure[i].host2vflib_map);
    if (VFLIB_Size>1) {
      for (j=0; j<VFLIB_Size; j++) {
        VF_Free(VF_DATA_Enclosure[i].comm_plan.partners[j].global_index);
      }
    }
    VF_Free(VF_DATA_Enclosure[i].comm_plan.partners);
    for (j=0; j<VF_DATA_Enclosure[i].npatches_l; j++) {
      VF_FreeSparseArray(&(VF_DATA_Enclosure[i].row[j].array0));
      VF_FreeSparseArray(&(VF_DATA_Enclosure[i].row[j].array1));
    }
    VF_Free(VF_DATA_Enclosure[i].row);
#ifdef AZTEC
    VF_Free(VF_DATA_Enclosure[i].data_org);
    VF_Free(VF_DATA_Enclosure[i].update_index);
    VF_Free(VF_DATA_Enclosure[i].external);
    VF_Free(VF_DATA_Enclosure[i].external_index);
#endif
  }
  VF_Free(VF_DATA_Enclosure);
  num_enclosures = 0;
}

int
VF_NewEnclosure(char* id, int debug_level)
{
  int i;

  for (i=0; i<num_enclosures; i++) {
    if (VF_DATA_Enclosure[i].id == NULL) {
      VF_DATA_Enclosure[i].id = VF_Newc(strlen(id)+1);
      strcpy(VF_DATA_Enclosure[i].id,id);
      break;
    }
    if (strcmp(VF_DATA_Enclosure[i].id, id) == 0) break;
  }
  VF_DATA_Enclosure[i].debug_level      = debug_level;
  VF_DATA_Enclosure[i].smoothed         = FALSE;
  VF_DATA_Enclosure[i].sym_method       = VF_SYMMETRIC_NONE;
  VF_DATA_Enclosure[i].seed             = -41731;
  VF_DATA_Enclosure[i].time_vf_init     = 0.0;
  VF_DATA_Enclosure[i].time_vf_calc     = 0.0;
  VF_DATA_Enclosure[i].time_vf_symmetry = 0.0;
  VF_DATA_Enclosure[i].time_vf_smooth   = 0.0;
  current_enclosure = i;
  return current_enclosure;
}

VFenclosure*
VF_FindEnclosure(char* id)
{
  int i;

  for (i=0; i<num_enclosures; i++) {
    if (VF_DATA_Enclosure[i].id == NULL) {
      VF_DATA_Enclosure[i].id = VF_Newc(strlen(id+1));
      strcpy(VF_DATA_Enclosure[i].id,id);
      break;
    }
    if (strcmp(VF_DATA_Enclosure[i].id, id) == 0) break;
  }
  current_enclosure = i;
  return &VF_DATA_Enclosure[current_enclosure];
}

VFenclosure*
VF_GetEnclosure(int enclosure)
{
  current_enclosure = enclosure;
  return &VF_DATA_Enclosure[enclosure];
}

VFenclosure*
VF_CurrentEnclosure()
{
  return &VF_DATA_Enclosure[current_enclosure];
}

VFtopology*
VF_CurrentTopology()
{
  return VF_DATA_Enclosure[current_enclosure].topology;
}

void
VF_SetEnclosure(int n)
{
  current_enclosure = n;
}

void
VF_RemapEnclosure (int index, int npatches, int global_ids[])
{
  float        *buffer;
  int          i, cnt, proc, size, *all_gids;
  VFenclosure* encl = VF_GetEnclosure(index);
  
  if (VFLIB_Size>1) {
    VF_GetSPbuffer0_ptr(&buffer);
    all_gids = (int*)buffer;
    for (size=0, proc=0; proc<VFLIB_Size; proc++) {
      if (proc==VFLIB_Rank) {
        cnt = npatches;
        for (i=0; i<npatches; i++) {
          all_gids[size+i] = global_ids[i];
        }
      }
      VF_BroadcastInt(&cnt, 1, proc);
      VF_BroadcastInt(&all_gids[size], cnt, proc);
      size += cnt;
    }
    VF_SortIntegerArray(all_gids, encl->npatches_g);
    if (encl->host_npatches != npatches) {
      encl->host2vflib_map = VF_ReNewi(encl->host2vflib_map, npatches);
      encl->host_npatches  = npatches;
    }
    for (i=0; i<encl->host_npatches; i++) {
      encl->host2vflib_map[i] = VF_LocateIndexInSortedArray(global_ids[i], 
                                                            all_gids, 
                                                            encl->npatches_g);
    }
  }
}

int
VF_NumPatches_l(int encl)
{
  VFenclosure *enclosure=VF_GetEnclosure(encl);
  return (enclosure->npatches_l);
}

int
VF_NumPatches_g(int encl)
{
  VFenclosure *enclosure=VF_GetEnclosure(encl);
  return (enclosure->npatches_g);
}
