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
@(#)    $RCSfile: vf_api.h,v $
@(#)    $Revision: 1.7 $  $Date: 2006/03/17 04:59:11 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/vf_api.h,v $
@(#)
@(#)--------------------------------------------------------------------------
*/

#ifndef _VF_API_H_

#define _VF_API_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "vf_mpi.h"

void   VF_Setup(int, int, MPI_Comm);
void   VF_SetNumEnclosures(int);
void   VF_SetMaxSurfaces(int);
void   VF_CleanUp(void);
void   VF_GetVersionNumber(char*);
void   VF_GetVersionDate(char*);
int    VF_NewEnclosure(char*, int);
int    VF_DefineEnclosure (char*, int, int, double, int, int*, int);
void   VF_RemapEnclosure (int, int, int*);
void   VF_RandomizeSurfacesOff(void);
void   VF_RandomizeSurfacesOn(void);
void   VF_DefineTopology(int, int, int, int, 
                         double*, double*, double*, int*, 
                         int, int*, int, int, int, int,
                         int, int, double, int);
void   VF_ResetTopology(int);
void   VF_CalcRead(int, char*);
void   VF_CalcHemicube(int, int, int, double);
void   VF_CalcPairwise (int, int, int, int, int, double, double);
void   VF_JitterOff(void);
void   VF_JitterOn(void);
void   VF_SmoothMatrix(int, double, double, int, int, int);
void   VF_RadSolve(int, double*, double*, double*, double, 
                   double, int, int*, int, int, int);
void   VF_RadSolveAux(int, double*, double*, double*, double*,  
                      double, double, int, int*, int, int, int);
void   VF_WriteGenesis(char*);
void   VF_WriteExodus(char*);
int    VF_WriteExodusAux1(char*, int, double);
void   VF_WriteExodusAux2(int, int, int, double, double*, double*);
void   VF_WriteExodusClose(int);
void   VF_WriteINP(char*);
void   VF_MatrixWrite(int, char*, int);

void   VF_OutputInitBanner(void);

void   VF_OutputCalcStartBanner(void);

void   VF_OutputCalcEndBanner(void);

void   VF_OutputSmoothingBanner(double wt, double tol, int max_iter);

void   VF_OutputMatrixSummaryBanner(void);

int    VF_NumPatches_l(int);

int    VF_NumPatches_g(int);

int    VF_MaxPatches(void);

int    VF_NumEnclosures(void);

double VF_GetMemoryUsage(void);

#ifdef __cplusplus
}
#endif

#endif
