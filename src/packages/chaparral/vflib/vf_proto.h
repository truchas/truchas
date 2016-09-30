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
@(#)    $RCSfile: vf_proto.h,v $
@(#)    $Revision: 1.13 $  $Date: 2006/03/17 04:59:11 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/vf_proto.h,v $
@(#)
@(#)--------------------------------------------------------------------------
*/

#ifndef _VF_PROTO_H_

#define _VF_PROTO_H_

#include <inttypes.h>

#include "vf_defines.h"
#include "vf_api.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN32
void sleep(clock_t);
#endif

void*    VF_Malloc_void(int, char*, int);
char*    VF_Malloc_char(int, char*, int);
int*     VF_Malloc_int(int, char*, int);
long*    VF_Malloc_long(int, char*, int);
float*   VF_Malloc_float(int, char*, int);
double*  VF_Malloc_double(int, char*, int);
void**   VF_MallocPtr_void(int, char*, int);
char**   VF_MallocPtr_char(int, char*, int);
int**    VF_MallocPtr_int(int, char*, int);
long**   VF_MallocPtr_long(int, char*, int);
float**  VF_MallocPtr_float(int, char*, int);
double** VF_MallocPtr_double(int, char*, int);
void*    VF_Realloc_void(void*, int, char*, int);
char*    VF_Realloc_char(char*, int, char*, int);
int*     VF_Realloc_int(int*, int, char*, int);
long*    VF_Realloc_long(long*, int, char*, int);
float*   VF_Realloc_float(float*, int, char*, int);
double*  VF_Realloc_double(double*, int, char*, int);
void     VF_MemoryInfo(char*);


void     VF_Exit(int);
void     VF_Sync(void);
void     VF_PrintSyncStart(void);
void     VF_PrintSyncEnd(void);
void     VF_StagedIO_Start(void);
void     VF_StagedIO_End(void);

void     VF_GlobalSum_Int(int*, char*, int);
void     VF_GlobalSum_UInt64(uint64_t*, char*, int);
void     VF_GlobalSum_Float(float*, char*, int);
void     VF_GlobalSum_Double(double*, char*, int);
void     VF_GlobalMin_Int(int*, char*, int);
void     VF_GlobalMin_Float(float*, char*, int);
void     VF_GlobalMin_Double(double*, char*, int);
void     VF_GlobalMax_Int(int*, char*, int);
void     VF_GlobalMax_Float(float*, char*, int);
void     VF_GlobalMax_Double(double*, char*, int);
void     VF_Sum_Int(int*, int, char*, int);
void     VF_Sum_Float(float*, int, char*, int);
void     VF_Sum_Double(double*, int, char*, int);
void     VF_Min_Int(int*, int, char*, int);
void     VF_Min_Float(float*, int, char*, int);
void     VF_Min_Double(double*, int, char*, int);
void     VF_Max_Int(int*, int, char*, int);
void     VF_Max_Float(float*, int, char*, int);
void     VF_Max_Double(double*, int, char*, int);
void     VF_Broadcast_Char(char*, int, int, char*, int);
void     VF_Broadcast_Int(int*, int, int, char*, int);
void     VF_Broadcast_Float(float*, int, int, char*, int);
void     VF_Broadcast_Double(double*, int, int, char*, int);
void     VF_Allgather_Int(int*, int*, int, char*, int);
void     VF_Allgather_Float(float*, float*, int, char*, int);
void     VF_Allgather_Double(double*, double*, int, char*, int);
void     VF_Allgatherv_Int(int*, int*, int, int*, int*, char*, int);
void     VF_Allgatherv_Double(double*, double*, int, int*, int*, char*, int);
void     VF_Exchange_Int(int*, char*, int);
void     VF_Exchange_Float(float*, char*, int);
void     VF_Exchange_Double(double*, char*, int);
void     VF_Gather_VFs(float*, int, char*, int);
void     VF_Gather_ExodusData(double*, double*, int*, int, int, char*, int);
void     VF_Gather_ExodusVals(double*, double*, double*, double*);
void     VF_Gather_RadiosityVectors(double*, double*, double*, 
                                    double*, double*, double*);
void     VF_Gather_RadiosityVectorsAux(double*, double*, double*, double*,  
                                       double*, double*, double*, double*);

void     VF_Error(char*);
double   VF_Clock(void);
int      VF_LocateIndexInArray(int, int*, int);
int      VF_LocateIndexInSortedArray(int, int*, int);
void     VF_SortIntegerArray(int*, int);
void     VF_InitializeSparseArray(VFsparse_array*);
void     VF_AllocateSparseArray(VFsparse_array*, int);
void     VF_ReallocateSparseArray(VFsparse_array*, int);
void     VF_FreeSparseArray(VFsparse_array*);
void     VF_SortIntegerArray(int*, int);
void     VF_SortSparseArray(VFsparse_array*);
void     VF_SortSparseArrayAux(int*, float*, int);
double   VF_SumSparseArray(VFsparse_array*);
void     VF_ExpandSparseArray(VFsparse_array*, float*);


void     VF_GetSPbuffer0_ptr(float**);
void     VF_GetSPbuffer1_ptr(float**);
void     VF_GetSPbuffer2_ptr(float**);
void     VF_GetDPbuffer0_ptr(double**);
void     VF_GetDPbuffer1_ptr(double**);
void     VF_GetDPbuffer2_ptr(double**);
void     VF_GetDPbuffer3_ptr(double**);
void     VF_GetINTbuffer_ptr(int**);
void     VF_FreeBuffers();


void     VF_InfoRead(int, char*);
void     VF_MatrixRead(int, char*);
void     VF_MatrixWrite(int, char*, int);
void     VF_AllInfoRead(char*);
void     VF_AllMatrixRead(char*, int);
void     VF_AllMatrixWrite(char*, int);


void     VF_LoadMatrixRow(int, float*, double);
void     VF_UpdateMatrixRow(int, float*);
void     VF_GetMatrixRow(int, int*, float*);
void     VF_GetMatrixCol(int, int, float*);
float    VF_GetMatrixEntry(int, int);
int      VF_GetMatrixRowGID(int);
float    VF_GetMatrixRowLastEntry(int);
void     VF_SetRawRowsum(int);
void     VF_SetSymmetricRowsum(int);
void     VF_SetSmoothedRowsum(int);
double   VF_GetRawRowsum(int);
double   VF_GetSmoothedRowsum(int);
double   VF_GetMatrixRowArea(int);
void     VF_GetMatrixAreas(double*);
void     VF_MakeMatrixSymmetric(int, int, int);
void     VF_RowsumStats(double*, double*, double*, double*, 
                        double*, double*, double*, double*,
                        double*, double*, double*, double*);
void     VF_SortMatrixRows(void);
VFrow*   VF_LocateLocalRow(int);
void     VF_ReciprocityStats(double*, double*);
void     VF_FillUpperDiagonal(void);
double   Ddot1(VFenclosure *e, double *v1, double *v2);
double   Ddot0(int n, double *v1, double *v2);
double   Ddot(int n, double *v1, double *v2);
void     Dscal(int n, double a, double r[]);
void     Daxpy1(VFenclosure *e, double a, double *x, double *y);
void     Daxpy2(VFenclosure *e, double a, double *x, double *y);
void     Daxpy(int n, double a, double *x, double *y);
void     Dcopy(int n, double *x, double *y);

void     VF_GetMatrix(int, int*, int*, float*, float*, float*);
void     VF_GetRowCounts(int, int, int*);
void     VF_GetRowCounts_Aux(int, int*);

void     VF_InitializeEnclosure(int);
void     VF_InitializeEnclosures(int);
void     VF_DeleteEnclosures(void);
void     VF_DeleteTopology(void);
void     VF_InitializeBuffers(int);

VFenclosure* VF_GetEnclosure(int);
VFenclosure* VF_FindEnclosure(char*);
VFenclosure* VF_CurrentEnclosure(void);
VFtopology*  VF_CurrentTopology(void);
VFtopology*  VF_GetTopology(void);
void         VF_SortFacets(void);

void   VF_TopologyConcat(int, int, int, double*, double*, double*, int*, int*,
                         double**, double**, double**, int**, int**, int*);

void   VF_BSP_CreateTree(VFtopology*);
void   VF_BSP_PartitionNode(BinNodePtr, int, int, int);
int    VF_BSP_FindPartition(BinNodePtr, int);
int    VF_BSP_IsFacetInNode(Facet*, double, int);
void   VF_BSP_DeleteTree(BinNodePtr);
Facet* VF_FirstFacetOfLinkList(FacetList*);
Facet* VF_NextFacetOfLinkList(FacetList*);

double VF_FacetArea(Facet*);
double VF_SurfArea(Poly*, int);
int    VF_IsFacetPlanar(Facet*);
double VF_PolyArea(Poly*);
double VF_GenericPolyArea(Poly*);
void   VF_PolyNormal(Poly*, Vector*);
void   VF_PolyNormal_Aux(Poly*, Vector*, int);
void   VF_GenericPolyNormal(Poly*, Vector*);
int    VF_IsPolyPlanar(Poly*);
void   VF_PolyCenter(Poly*, Vector*);
void   VF_FacetToPoly(Facet*, Poly*);
void   VF_FacetToPolyStack(Facet*, POLYstack*, int*);
void   VF_SubPoly(int, int, int, int, Poly*, Poly*);
void   VF_PrintPoly(Poly*, char*);
void   VF_PrintVector(Vector*, char*);
int    VF_SharedEdge(Facet*, Facet*, Poly*, Poly*);
int    VF_SharedPolyEdge(Poly*, Poly*);
void   VF_UV_to_XYZ(Poly*, double, double, Vector*);
Point  VF_UVtoXYZ(double, double, Poly*);
double VF_GetPolyJacobian(Poly*, double, double);
void   VF_GetShapeFunction(Poly*, double, double, double*, double*);
void   VF_TransformPoly(Poly*, Poly*, Matrix*);
double VF_QuadAspectRatio(Poly*);
double VF_TriAspectRatio(Poly*, double*, double*, double*);


void   VF_ComputeFacetBoundingBox(Box*, Facet*);
void   VF_ComputePolyBoundingBox(Box*, Poly*);
void   VF_ExtentsBoundingBox(Box*, Box*, Box*, Box*);



int    VF_BackFaceCullPolys(Poly*, Poly*);
int    VF_BackFaceCullViewPoint(ViewPort*, Poly*);
int    VF_BehindPoly(Poly*, Poly*);
int    VF_BehindPlane(Plane*, Poly*);
int    VF_BehindViewPoint(ViewPort*, Poly*);
int    VF_BehindAndBackFaceCullViewPoint(ViewPort*, Poly*);
int    VF_VoxelBehindPlane(Box*, Plane*);
int    VF_VoxelBehindFacet(Box*, Facet*);
int    VF_VoxelBehindPoly(Box*, Poly*);


int    VF_ClipToFrustum(Poly*, int);
void   VF_ClipToFrustumPlane(Poly*, Poly*, double d[]);
void   VF_ClipToPolyPlane(Poly*, Poly*, Poly *poly);
int    VF_ClipBoxToFrustum(ViewPort*, Box*, int);
void   VF_PolyScan(Poly *p0, Poly*, Window*, int, int, Hemicube*);


void   VF_MatrixCalcHemicube(void);
void   VF_AllocHemicube(Hemicube*);
void   VF_FreeHemicube(Hemicube*);
void   VF_JitterHemicube(Vector*, Vector*, Vector*, Poly*);

void   VF_HemicubeSub(int, int, Poly*, ViewPort*, 
                      double*, double*, int, int);
void   VF_HemicubeProjectRow(int, int, Poly*, ViewPort*, 
                             int, double, double*, double*);
void   VF_ProjectOntoHemicube(Poly*, Poly*, int, ViewPort*, int, int, Hemicube*);
double VF_FindMinSeperationDist(Poly*, int, int);
void   VF_SetView(ViewPort*, Poly*);
void   VF_SetViewPort(ViewPort*, int);



void   VF_MatrixCalcPairwise(void);
void   VF_InitCandidates(Facet *facet);
int    VF_FindCandidates(Poly *poly_i, Poly *poly_j,
                         CandidateList **candidates);
void   VF_SamplePoly(Poly *poly, Sampling *sampling, Point samples[]);
void   VF_SampleRandom(Poly *poly, Sampling *sampling, Point samples[]);
void   VF_SampleUniform(Poly *poly, Sampling *sampling, Point samples[]);
void   VF_SampleJitter(Poly *poly, Sampling *sampling, Point samples[]);
void   VF_SampleHalton(Poly *poly, Sampling *sampling, Point samples[]);
void   VF_SetupSampling(Sampling *sampling, Point uv_samples[]);
void   VF_SetupSampleRandom(Sampling *sampling, Point uv_samples[]);
void   VF_SetupSampleUniform(Sampling *sampling, Point uv_samples[]);
void   VF_SetupSampleJitter(Sampling *sampling, Point uv_samples[]);
void   VF_SetupSampleHalton(Sampling *sampling, Point uv_samples[]);
void   VF_ShuffleSamples(int ns, Point samples[]);
void   VF_SamplesUVtoXYZ(Poly *poly, Sampling *sampling, 
                    Point uv_samples[], Point xyz_samples[]);
void   VF_ConvertUVtoXYZ(Poly *poly, int ns, Point samples[]);
void   VF_CreateShaft(Shaft *shaft, Poly *poly_i, Poly *poly_j);
int    VF_ShaftCull(Box *bounds, Shaft *shaft, double tol);
double VF_Visibility(Poly *poly_i, Poly *poly_j, 
                     CandidateList *candidates, int nc);
int    VF_RayPatchListTest(Ray *ray, int sample_num, double tmax,
                    CandidateList *candidates, int ncandidates);
int    VF_RayIntersect(Ray *ray, double tmax, Poly *poly);
int    VF_InsideLine(Ray *ray, double t, Poly *poly);
int    VF_InsidePoly(Ray *ray, double t, Poly *poly);
double VF_CalcVF_ComputePair(Poly*, Poly*);
double VF_CalcVF_Unoccluded(Poly *poly_i, Poly *poly_j);
double VF_CalcVF_Occluded(Poly *poly_i, Poly *poly_j, CandidateList *candidates,
                          int ncandidates, double visibility);
double VF_CalcVF_Analytic(Poly *poly_i, Poly *poly_j);
double VF_CalcVF_Contour(Poly *poly_i, Poly *poly_j, int nseg_i, int nseg_j);
double VF_CalcVF_DoubleArea (Poly *poly_i, Poly *poly_j, 
                      CandidateList *candidates, int ncandidates, 
                      int nseg_i, int nseg_j);
double VF_CalcVF_Gauss(Poly *poly_i, Poly *poly_j, int nseg_i, int nseg_j);
double VF_CalcVF_MonteCarlo (Poly *poly_i, Poly *poly_j, 
                             CandidateList *candidates, 
                             int ncandidates);
double VF_CalcVF_Hemicube (Poly *poly_i, Poly *poly_j, 
                           CandidateList *candidates, 
                           int ncandidates);
double VF_CalcVF_Hemicube0(Poly *poly_i, Poly *poly_j);
double VF_CalcVF_HemicubeRow (Facet *facet);
double VF_CalcVF_Hottel(Poly *poly_i, Poly *poly_j);


void   VF_PartialEnclosure(void);

void   VF_RadSolveCleanUp(void);
void   VF_RadSolveCleanupCG(void);
void   VF_RadSolveCleanupGMRES(void);
void   VF_RadSolveCleanupAZTEC(void);

void   VF_RadSolveCleanpCG(void);
void   VF_RadSolveCleanpGMRES(void);
void   VF_RadSolveCleanpAZTEC(void);

void   VF_RadSolveCG(double*, double*, double*, double*, 
                     double, double, int, int*, int);

void   VF_RadSolveGMRES(double*, double*, double*, double*, 
                        double, double, int, int*, int);

void   VF_RadSolveAZTEC(double*, double*, double*, double*, 
                        double, double, int, int*, int, int, int);

void   VF_ComputeFluxes(double*, double*, double*, double*, int);

void   VF_RadSolve_Residual(double*, double*, double*, double*);
void   VF_RadSolve_MatVecMul(double*, double*, double*);
void   VF_RadSolve_Residual_Scale(double*, double*, double*, double*, double*);
void   VF_RadSolve_MatVecMul_Scale(double*, double*, double*, double*);

double ran2(long*);

#ifdef VF_INITIALIZE

MPI_Comm VFLIB_Comm;
int      VFLIB_Size=1;
int      VFLIB_Rank=0;
int      VFLIB_StagedIO=0;
int      VFLIB_Ncntrls=1;

#else

extern MPI_Comm VFLIB_Comm;
extern int VFLIB_Size;
extern int VFLIB_Rank;
extern int VFLIB_StagedIO;
extern int VFLIB_Ncntrls;

#endif

#ifdef __cplusplus
}
#endif

#endif
