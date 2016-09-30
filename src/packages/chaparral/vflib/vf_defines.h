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
@(#)    $RCSfile: vf_defines.h,v $
@(#)    $Revision: 1.10.2.6 $  $Date: 2006/12/06 22:05:50 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/vf_defines.h,v $
@(#)
@(#)--------------------------------------------------------------------------
*/

#ifndef _VF_DEFINES_H_

#define _VF_DEFINES_H_

#include "vf_mpi.h"

#include <math.h>
#if defined(WIN32) || defined(macosx)
#  include <float.h>
#else
#  include <values.h>
#endif

#define VF_VERSION     "3.2.4"
#define VF_DATE        "Dec 06, 2006"

#ifndef M_PI
#  define M_PI      3.14159265358979323846
#endif

#ifndef M_1_PI
#  define M_1_PI    0.31830988618379067154
#endif

#ifndef M_LN2
#  define M_LN2     0.69314718055994530942
#endif

#ifndef MAXINT
#  define MAXINT    2147483647 
#endif

#ifndef MAXDOUBLE
#  ifdef DBL_MAX
#    define MAXDOUBLE DBL_MAX
#  else
#    ifdef BIGNUM
#      define MAXDOUBLE BIGNUM
#    else
#      define MAXDOUBLE 1.7976931348623157e+308
#    endif
#  endif
#endif

#define SQ2         0.70710678118654752440

typedef struct Point3Struct {
  double x, y, z;
} Point;

typedef Point Vector;

typedef struct MatrixStruct { 
  double element[4][4];
} Matrix;

#ifndef TRUE
#  define TRUE      1
#endif

#ifndef FALSE
#  define FALSE     0
#endif

/* find minimum of a and b */
#define MIN(a,b)    (((a)<(b))?(a):(b))

/* find maximum of a and b */
#define MAX(a,b)    (((a)>(b))?(a):(b))

/* take sign of a, either -1, 0, or 1 */
#define ZSGN(a)     (((a)<0) ? -1 : (a)>0 ? 1 : 0)

/* take binary sign of a, either -1, or 1 if >= 0 */
#define SGN(a)      (((a)<0) ? -1 : 1)

/* shout if something that should be true isn't */
#define ASSERT(x) \
if (!(x)) fprintf(stderr," Assert failed: x\n");

/***********************************/
/* mocros for some vector routines */
/***********************************/

/*** 2D ***/

#define V2_Zero(a) \
        (a)->x = 0.0;  (a)->y = 0.0;  (a)->z = 0.0

#define V2_Init(a,b)     \
        (b)->x = (a)->x; \
        (b)->y = (a)->y

#define V2_Sub(a,b,c)             \
        (c)->x = (a)->x - (b)->x; \
        (c)->y = (a)->y - (b)->y

#define V2_Add(a,b,c)             \
        (c)->x = (a)->x + (b)->x; \
        (c)->y = (a)->y + (b)->y

#define V2_Scale(a,b)          \
        (a)->x = (a)->x * (b); \
        (a)->y = (a)->y * (b)

#define V2_Mul(a,b,c)             \
        (c)->x = (a)->x * (b)->x; \
        (c)->y = (a)->y * (b)->y
        
#define V2_Negate(a)      \
        (a)->x = -(a)->x; \
        (a)->y = -(a)->y

#define V2_Length(a) \
        sqrt((a)->x * (a)->x + (a)->y * (a)->y)

#define V2_Dot(a,b) \
        ((a)->x * (b)->x + (a)->y * (b)->y)

#define V2_Cross(a,b)             \
        ((a)->x * (b)->y - (a)->y * (b)->x)

#define V2_DistanceBetween2Points(a,b) \
        sqrt(((a)->x - (b)->x) * ((a)->x - (b)->x) + \
             ((a)->y - (b)->y) * ((a)->y - (b)->y))

#define V2_Normalize(v,len)     \
        {                       \
          len = V2_Length((v)); \
          if (len != 0.0) {     \
            (v)->x /= len;      \
            (v)->y /= len;      \
          }                     \
        }
        
/*** 3D ***/

#define V3_Zero(a) \
        (a)->x = 0.0;  (a)->y = 0.0;  (a)->z = 0.0

#define V3_Init(a,b)     \
        (b)->x = (a)->x; \
        (b)->y = (a)->y; \
        (b)->z = (a)->z

#define V3_Sub(a,b,c)             \
        (c)->x = (a)->x - (b)->x; \
        (c)->y = (a)->y - (b)->y; \
        (c)->z = (a)->z - (b)->z

#define V3_Add(a,b,c)             \
        (c)->x = (a)->x + (b)->x; \
        (c)->y = (a)->y + (b)->y; \
        (c)->z = (a)->z + (b)->z

#define V3_Scale(a,b)          \
        (a)->x = (a)->x * (b); \
        (a)->y = (a)->y * (b); \
        (a)->z = (a)->z * (b)

#define V3_Mul(a,b,c)             \
        (c)->x = (a)->x * (b)->x; \
        (c)->y = (a)->y * (b)->y; \
        (c)->z = (a)->z * (b)->z
        
#define V3_Negate(a)      \
        (a)->x = -(a)->x; \
        (a)->y = -(a)->y; \
        (a)->z = -(a)->z

#define V3_Length(a) \
        sqrt((a)->x * (a)->x + (a)->y * (a)->y + (a)->z * (a)->z)

#define V3_Dot(a,b) \
        ((a)->x * (b)->x + (a)->y * (b)->y + (a)->z * (b)->z)

#define V3_Cross(a,b,c)             \
        (c)->x = ((a)->y*(b)->z) - ((a)->z*(b)->y); \
        (c)->y = ((a)->z*(b)->x) - ((a)->x*(b)->z); \
        (c)->z = ((a)->x*(b)->y) - ((a)->y*(b)->x)

#define V3_DistanceBetween2Points(a,b) \
        sqrt(((a)->x - (b)->x) * ((a)->x - (b)->x) + \
             ((a)->y - (b)->y) * ((a)->y - (b)->y) + \
             ((a)->z - (b)->z) * ((a)->z - (b)->z))

#define V3_Normalize(v,len)     \
        {                       \
          len = V3_Length((v)); \
          if (len != 0.0) {     \
            (v)->x /= len;      \
            (v)->y /= len;      \
            (v)->z /= len;      \
          }                     \
        }
        
/**********************************/
/* additional one-argument macros */
/**********************************/

/* round a to nearest integer towards 0 */
#define FLOOR(a)                ((a)>0 ? (int)(a) : -(int)(-a))

/* round a to nearest integer away from 0 */
#define CEILING(a) \
((a)==(int)(a) ? (a) : (a)>0 ? 1+(int)(a) : -(1+(int)(-a)))


/*******************************/
/* additional useful constants */
/*******************************/

#define BIG_NUMBER      1000000.0
#define SMALL_NUMBER    0.00000001

#define RANDOM    ran2(&(encl->seed))
#define RAND(l,h) (((h)-(l))+(l))*RANDOM
#define RANI(l,h) ((int)(RAND(0,(h)-(l)+1)+(l)))

#define VF_INVALID            -1

#define VF_OUTPUT_HISTOGRAM  256
#define VF_OUTPUT_NONE         0
#define VF_OUTPUT_SUMMARY      1
#define VF_OUTPUT_VERBOSE      2
#define VF_OUTPUT_DEBUG_0      3
#define VF_OUTPUT_DEBUG_1      4
#define VF_OUTPUT_DEBUG_2      5
#define VF_OUTPUT_DEBUG_3      6
#define VF_OUTPUT_DEBUG_4      7

#define VF_FILE_READ           0
#define VF_PAIRWISE            1
#define VF_HEMICUBE            2
#define VF_HYBRID              3

#define VF_2Daxisym            1
#define VF_2Dplanar            2
#define VF_3D                  3

#define VF_LINE                2
#define VF_TRI                 3
#define VF_QUAD                4
   
#define VF_RANDOM_SAMPLE       0
#define VF_UNIFORM_SAMPLE      1
#define VF_JITTER_SAMPLE       2
#define VF_HALTON_SAMPLE       3

#define VF_SYMMETRIC_NONE      0
#define VF_SYMMETRIC_SUB       1
#define VF_SYMMETRIC_ADD       2
#define VF_SYMMETRIC_AVG       3

#define VF_BSP_FRONT           0
#define VF_BSP_BACK            1
#define VF_BSP_SPLIT           2
#define VF_BSP_SPLITON         3

#define VF_BSP_XAXIS           1
#define VF_BSP_YAXIS           2
#define VF_BSP_ZAXIS           4

#define VF_RAY_TEST  1.0e-10

#define VF_FRUSTUM_CLIP_UNKNOWN  0     /* entity not yet clipped by frustum */
#define VF_FRUSTUM_CLIP_OUT      1     /* entity entirely outside frustum   */
#define VF_FRUSTUM_CLIP_PARTIAL  2     /* entity partially inside frustum   */
#define VF_FRUSTUM_CLIP_IN       3     /* entity entirely inside frustum    */

#define VF_POLY_CLIP_OUT      0     /* polygon entirely outside box */
#define VF_POLY_CLIP_PARTIAL  1     /* polygon partially inside */
#define VF_POLY_CLIP_IN       2     /* polygon entirely inside box */

#define VF_POLY_NMAX         20

#define VF_FACET_MASK_PLANAR        1
#define VF_FACET_MASK_SPLIT         2
#define VF_FACET_MASK_DOIT          4
#define VF_FACET_MASK_ABOVE         8
#define VF_FACET_MASK_TOP          16
#define VF_FACET_MASK_LEFT         32
#define VF_FACET_MASK_RIGHT        64
#define VF_FACET_MASK_FRONT       128
#define VF_FACET_MASK_BACK        256
#define VF_FACET_MASK_HYBRID      512

#define VF_MATRIX_EXCL_VIRT 1 /* exclude the virtual vf when getting the matrix */
#define VF_MATRIX_EXCL_DIAG 2 /* exclude the diagonal when getting the matrix */

/*==========================================*/
/* DEFINE TABLE OF EXTERNAL STORAGE OPTIONS */
/*==========================================*/
#define VF_ASCII       0   /* External storage in ASCII text format          */
#define VF_BINARY      1   /* External storage in machine dependent binary   */
#define VF_XDRFMT      2   /* External storage in external data rep. format  */

#define VF_RADSOLVE_CG           1
#define VF_RADSOLVE_GMRES        2
#define VF_RADSOLVE_AZTEC_CG     3
#define VF_RADSOLVE_AZTEC_GMRES  4

#define VF_Newv(a)     VF_Malloc_void      ((a),__FILE__,__LINE__)
#define VF_Newc(a)     VF_Malloc_char      ((a),__FILE__,__LINE__)
#define VF_Newi(a)     VF_Malloc_int       ((a),__FILE__,__LINE__)
#define VF_Newl(a)     VF_Malloc_long      ((a),__FILE__,__LINE__)
#define VF_Newf(a)     VF_Malloc_float     ((a),__FILE__,__LINE__)
#define VF_Newd(a)     VF_Malloc_double    ((a),__FILE__,__LINE__)
#define VF_NewvPtr(a)  VF_MallocPtr_void   ((a),__FILE__,__LINE__)
#define VF_NewcPtr(a)  VF_MallocPtr_char   ((a),__FILE__,__LINE__)
#define VF_NewiPtr(a)  VF_MallocPtr_int    ((a),__FILE__,__LINE__)
#define VF_NewlPtr(a)  VF_MallocPtr_long   ((a),__FILE__,__LINE__)
#define VF_NewfPtr(a)  VF_MallocPtr_float  ((a),__FILE__,__LINE__)
#define VF_NewdPtr(a)  VF_MallocPtr_double ((a),__FILE__,__LINE__)
#define VF_ReNewv(a,b) VF_Realloc_void     ((a),(b),__FILE__,__LINE__)
#define VF_ReNewc(a,b) VF_Realloc_char     ((a),(b),__FILE__,__LINE__)
#define VF_ReNewi(a,b) VF_Realloc_int      ((a),(b),__FILE__,__LINE__)
#define VF_ReNewl(a,b) VF_Realloc_long     ((a),(b),__FILE__,__LINE__)
#define VF_ReNewf(a,b) VF_Realloc_float    ((a),(b),__FILE__,__LINE__)
#define VF_ReNewd(a,b) VF_Realloc_double   ((a),(b),__FILE__,__LINE__)
#define VF_Free(a)     if (a!=NULL) {free(a); a = NULL;}


#define VF_GlobalSumInt(a)               VF_GlobalSum_Int((a),__FILE__,__LINE__)
#define VF_GlobalSumUInt64(a)            VF_GlobalSum_UInt64((a),__FILE__,__LINE__)
#define VF_GlobalSumFloat(a)             VF_GlobalSum_Float((a),__FILE__,__LINE__)
#define VF_GlobalSumDouble(a)            VF_GlobalSum_Double((a),__FILE__,__LINE__)
#define VF_GlobalMinInt(a)               VF_GlobalMin_Int((a),__FILE__,__LINE__)
#define VF_GlobalMinFloat(a)             VF_GlobalMin_Float((a),__FILE__,__LINE__)
#define VF_GlobalMinDouble(a)            VF_GlobalMin_Double((a),__FILE__,__LINE__)
#define VF_GlobalMaxInt(a)               VF_GlobalMax_Int((a),__FILE__,__LINE__)
#define VF_GlobalMaxFloat(a)             VF_GlobalMax_Float((a),__FILE__,__LINE__)
#define VF_GlobalMaxDouble(a)            VF_GlobalMax_Double((a),__FILE__,__LINE__)
#define VF_SumInt(a,b)                   VF_Sum_Int((a),(b),__FILE__,__LINE__)
#define VF_SumFloat(a,b)                 VF_Sum_Float((a),(b),__FILE__,__LINE__)
#define VF_SumDouble(a,b)                VF_Sum_Double((a),(b),__FILE__,__LINE__)
#define VF_MinInt(a,b)                   VF_Min_Int((a),(b),__FILE__,__LINE__)
#define VF_MinFloat(a,b)                 VF_Min_Float((a),(b),__FILE__,__LINE__)
#define VF_MinDouble(a,b)                VF_Min_Double((a),(b),__FILE__,__LINE__)
#define VF_MaxInt(a,b)                   VF_Max_Int((a),(b),__FILE__,__LINE__)
#define VF_MaxFloat(a,b)                 VF_Max_Float((a),(b),__FILE__,__LINE__)
#define VF_MaxDouble(a,b)                VF_Max_Double((a),(b),__FILE__,__LINE__)
#define VF_BroadcastChar(a,b,c)          VF_Broadcast_Char((a),(b),(c),__FILE__,__LINE__)
#define VF_BroadcastInt(a,b,c)           VF_Broadcast_Int((a),(b),(c),__FILE__,__LINE__)
#define VF_BroadcastFloat(a,b,c)         VF_Broadcast_Float((a),(b),(c),__FILE__,__LINE__)
#define VF_BroadcastDouble(a,b,c)        VF_Broadcast_Double((a),(b),(c),__FILE__,__LINE__)
#define VF_AllgatherInt(a,b,c)           VF_Allgather_Int((a),(b),(c),__FILE__,__LINE__)
#define VF_AllgatherFloat(a,b,c)         VF_Allgather_Float((a),(b),(c),__FILE__,__LINE__)
#define VF_AllgatherDouble(a,b,c)        VF_Allgather_Double((a),(b),(c),__FILE__,__LINE__)
#define VF_AllgathervInt(a,b,c,d,e)      VF_Allgatherv_Int((a),(b),(c),(d),(e),__FILE__,__LINE__)
#define VF_AllgathervDouble(a,b,c,d,e)   VF_Allgatherv_Double((a),(b),(c),(d),(e),__FILE__,__LINE__)
#define VF_ExchangeInt(a)                VF_Exchange_Int((a),__FILE__,__LINE__)
#define VF_ExchangeFloat(a)              VF_Exchange_Float((a),__FILE__,__LINE__)
#define VF_ExchangeDouble(a)             VF_Exchange_Double((a),__FILE__,__LINE__)
#define VF_GatherVFs(a,b)                VF_Gather_VFs((a),(b),__FILE__,__LINE__)
#define VF_GatherExodusData(a,b,c,d,e)   VF_Gather_ExodusData((a),(b),(c),(d),(e),__FILE__,__LINE__)


/*=======================*/
/* POLY STACK OPERATIONS */
/*=======================*/
#define VF_POLY_STACKSIZE  50

#define VF_InitPolyStack(a)        (a)->cnt = 0;

#define VF_PushPolyStack(a,b)      (a)->poly[(a)->cnt] = b; ((a)->cnt)++;

#define VF_PopPolyStack(a,b)       ((a)->cnt)--; (b) = (a)->poly[(a)->cnt];

#define VF_PopPolyStackPtr(a,b)    ((a)->cnt)--; (b) = &((a)->poly[(a)->cnt]);

#define VF_PushSplitPolyStack(a,b) (a)->poly[(a)->cnt].np     = 4;\
                                (a)->poly[(a)->cnt].p[0]   = (b)->p[0];\
                                (a)->poly[(a)->cnt].p[1]   = (b)->p[1];\
                                (a)->poly[(a)->cnt].p[2]   = (b)->p[2];\
                                (a)->poly[(a)->cnt].p[3]   = (b)->p[3];\
                                (a)->poly[(a)->cnt].d      = (b)->d;\
                                (a)->poly[(a)->cnt].normal = (b)->normal;\
                                ((a)->cnt)++;\
                                (a)->poly[(a)->cnt].np     = 3;\
                                (a)->poly[(a)->cnt].p[0]   = (b)->p[3];\
                                (a)->poly[(a)->cnt].p[1]   = (b)->p[4];\
                                (a)->poly[(a)->cnt].p[2]   = (b)->p[0];\
                                (a)->poly[(a)->cnt].d      = (b)->d;\
                                (a)->poly[(a)->cnt].normal = (b)->normal;\
                                ((a)->cnt)++;

/*===========================*/
/* BSP TREE STACK OPERATIONS */
/*===========================*/
#define VF_BSP_STACKSIZE  100

#define VF_BSP_InitStack(a)    (a)->stack[0].node = NULL; (a)->stackPtr = 1;

#define VF_BSP_PushStack(a,b)  (a)->stack[(a)->stackPtr].node = b; ((a)->stackPtr)++;

#define VF_BSP_PopStack(a,b)   ((a)->stackPtr)--; (b) = (a)->stack[(a)->stackPtr].node;

/*===========================*/
/* DEFINE VARIOUS DATA TYPES */
/*===========================*/
typedef struct IPointStruct {
  int x, y, z;
} IPoint;

typedef struct RayStruct {
  int    vis;
  Point  O;
  Vector D;
} Ray;

typedef struct PlaneStruct {  
  Vector normal;
  double d;
} Plane;

typedef struct BoxStruct {
  double xmin,ymin,zmin;
  double xmax,ymax,zmax;
} Box;

typedef struct FacetStruct {
  int    index;
  int    patch_proc;
  int    patch_gid;
  int    patch_local_index;
  int    patch_global_index;
  int    num_vertices;
  int    *vertex_list;
  int    mask;
  int    sector;
  Vector normal;
  double d;
  Box    bounds;  
} Facet;

struct FacetLink {
  char   type;
  Facet* facet;
  struct FacetLink *next;
};

typedef struct {
  struct FacetLink *head;
  struct FacetLink *tail;
  struct FacetLink *current;
  int    length;
} FacetList;

typedef struct BinNode {
  int       axis;               /* 0=x, 1=y, 2=z                      */
  double    d;                  /* plane distance                     */
  Box       bounds;             /* bounding box for the node          */
  int       behind;             /* behind patch                       */
  int       clipped;            /* frustum clipped                    */
  int       depth;              /* node depth                         */
  FacetList members;            /* list of enclosed patches           */
  struct    BinNode *child[2];  /* pointers to children nodes, if any */
} BinNode, *BinNodePtr;

typedef struct {
  Box        bounds;        /* extent of the entire bin tree */
  FacetList  members;       /* list of all of the patches */
  int        MaxDepth;      /* max allowed depth of the tree */
  int        MaxListLength; /* max primitive allowed in a leaf node */
  BinNodePtr root;          /* root of the entire bin tree */
  struct     FacetLink *facetLinkList;
} BinTree;

typedef struct {
  BinNodePtr node;
} BSPstackElem;

typedef struct {
  int          stackPtr;
  BSPstackElem stack[VF_BSP_STACKSIZE];
} BSPstack, *BSPstackPtr;

typedef struct PolyStruct {
  int    facet;
  int    sub_facet;
  int    np;
  int    mask;
  Point  p[VF_POLY_NMAX];
  Vector normal;
  double d; 
} Poly, *PolyPtr;

typedef struct PolyStackStruct {
  int      cnt;
  Poly     poly[VF_POLY_STACKSIZE];
} POLYstack, *POLYstackPtr;


typedef struct PolyTreeStruct {
  int    level;
  Poly   parent;
  double area;
  int    nchildren;
  struct PolyTreeStruct *child;   
} PolyTree;

typedef struct PolyBoxStruct {
  double x0, x1;
  double y0, y1;
  double z0, z1;
} PolyBox;

typedef struct WindowStruct {
  int x0, y0;
  int x1, y1;
} Window;

typedef struct ViewPortStruct {
  Vector view_point;
  Vector view_normal;
  Vector up_vector;
  Vector n[5];
  Vector u[5];
  double d;
  Window window;
  Matrix xform;
  Plane  frustum_planes0[2][5];
  Plane  frustum_planes1[5];
} ViewPort;

typedef struct ShaftStruct {
  Box   box_i;
  Box   box_j;
  Box   extent;
  int   nplanes;
  Plane plane[8];
} Shaft;

typedef struct CandidatesStruct {
  Poly poly;
  int  facet_num;
  int  poly_num;
  int  cached;
  int  check_it;
} CandidateList;

typedef struct HemicubeStruct {
  int    sub_divide;
  double min_separation;
  double min_distance;
  double min_equiv_radius;
  int    imax;
  int    jmax;
  int    resolution;
  int    nholes;
  float  *side_deltaVF;
  float  *top_deltaVF;
  Vector *side_dir;
  Vector *top_dir;
  Vector *dir;
  double *zbuffer;
  int    *ibuffer;
  double *vf;
} Hemicube;

typedef struct SampleStruct {
  int method;
  int index;
  int n;
} Sampling;

typedef struct MonteCarloStruct {
  Sampling sampling;
  Point    *uv_samples;
  double   tol1;
  double   tol2;
} MonteCarlo;

typedef struct VisStruct {
  Sampling sampling;
  Point    *sampleuv_i;
  Point    *sampleuv_j;
  Ray      *ray;
} VisCheck;

typedef struct AdaptiveStruct {
  int        num_numerical;
  int        num_analytic;
  int        num_zero;
  int        num_total;
  int        num_numerical_row;
  int        num_analytic_row;
  int        num_zero_row;
  int        num_total_row;
  MonteCarlo montecarlo;
  VisCheck   visibility;
  Hemicube   hemicube;
  int        hc_mode;
} Adaptive;
 
typedef struct TopologyStruct {
  int      geom;
  int      nonblocking;
  int      nnodes;
  int      nfacets_g;
  int      nfacets_l;
  int      nfacets_base;
  int      nrotations;
  int      nsections;
  int      xmirror;
  int      ymirror;
  int      zmirror;
  int      vertex_offset;
  int      surf_offset;
  Box      bounds;
  double   theta;
  double   spatial_tol;
  double   *x, *y, *z;
  int      *connect;
  Facet    *facets;
    
  int      bsp_depth;
  int      bsp_length;
  BinTree  BSP_tree;
  BSPstack *BSP_stack;
    
  int      debug_level;
} VFtopology;

typedef struct SparseArrayStruct {
  int   cnt;
  int   size;
  float *data;
  int   *index;
} VFsparse_array;

typedef struct VFrowStruct {           /* A VIEWFACTOR ROW */
  int            host_gid;
  int            global_index;
  int            local_index;
  double         area;
  double         raw_rowsum;
  double         sym_rowsum;
  double         smooth_rowsum;
  float          diagonal;
  VFsparse_array array0;
  VFsparse_array array1;
} VFrow;

typedef struct VFcommpartnerStruct {
  int proc;
  int cnt;
  int offset;
  int *global_index;
} VFcommPartner;

typedef struct VFcommplanStruct {
  int npartners;
  VFcommPartner *partners;
} VFcommPlan;

typedef struct VFenclStruct {          /* AN ENCLOSURE */
  int        initialized;
  char*      id;
  int        nonblocking;
  int        partial;
  double     asink;
  int        npatches_g;
  int        npatches_l;
  int        host_npatches;
  int        *host2vflib_map;
  int        vf_method;
  int        sym_method;
  int        smoothed;
  VFrow      *row;
  VFtopology *topology;
  Hemicube   hemicube;
  Adaptive   adaptive;
  long       seed;
  int        debug_level;
#ifdef AZTEC
  int        *data_org;
  int        *update_index;
  int        *external;
  int        *external_index;
#endif
  double     time_vf_init;
  double     time_vf_calc;
  double     time_vf_symmetry;
  double     time_vf_smooth;
  double     time_vf_solve;
  VFcommPlan comm_plan;
} VFenclosure;

#endif
