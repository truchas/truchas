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
@(#)    $RCSfile: VF_Initialize.c,v $
@(#)    $Revision: 1.9 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_Initialize.c,v $
@(#)
@(#)    DESCRIPTION:  Misc. parallel utility functions.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define VF_INITIALIZE
#include "vf.h"
#undef  VF_INITIALIZE

void
VF_Setup(int staged_io, int ncntrls, MPI_Comm Comm)
{
  VFLIB_StagedIO = staged_io;
  VFLIB_Ncntrls  = ncntrls;
#ifndef VF_NO_MPI
  MPI_Comm_dup  (Comm,&VFLIB_Comm);
  MPI_Comm_rank (VFLIB_Comm,&VFLIB_Rank);
  MPI_Comm_size (VFLIB_Comm,&VFLIB_Size);
#else
  VFLIB_Rank = 0;
  VFLIB_Size = 1;
#endif
}

void
VF_SetNumEnclosures(int nenclosures)
{
  VF_InitializeEnclosures(nenclosures);
}

void
VF_SetMaxSurfaces(int max_surfaces)
{
  VF_InitializeBuffers(max_surfaces);
}

void
VF_CleanUp()
{
  VF_DeleteEnclosures();
  VF_DeleteTopology();
  VF_FreeBuffers();
  VF_RadSolveCleanUp();
}

void
VF_GetVersionNumber(char *version)
{
  strcpy(version,VF_VERSION);
}

void
VF_GetVersionDate(char *date)
{
  strcpy(date,VF_DATE);
}
