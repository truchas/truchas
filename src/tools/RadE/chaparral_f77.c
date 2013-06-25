#include <stdlib.h>
#include <string.h>

#include "vf_api.h"

#include <FortranCInterface_names.h>

#define VF_SETUP            TR_ROUTINE_GLOBAL_(vf_setup,VF_SETUP)
#define VF_SETNUMENCLOSURES TR_ROUTINE_GLOBAL_(vf_setnumenclosures,VF_SETNUMENCLOSURES)
#define VF_SETMAXSURFACES   TR_ROUTINE_GLOBAL_(vf_setmaxsurfaces,VF_SETMAXSURFACES)
#define VF_CLEANUP          TR_ROUTINE_GLOBAL_(vf_cleanup,VF_CLEANUP)
#define VF_DEFINEENCLOSURE  TR_ROUTINE_GLOBAL_(vf_defineenclosure,VF_DEFINEENCLOSURE)
#define VF_RANDOMIZESURFACESON TR_ROUTINE_GLOBAL_(vf_randomizesurfaceson,VF_RANDOMIZESURFACESON)
#define VF_RANDOMIZESURFACESOFF TR_ROUTINE_GLOBAL_(vf_randomizesurfacesoff,VF_RANDOMIZESURFACESOFF)
#define VF_DEFINETOPOLOGY   TR_ROUTINE_GLOBAL_(vf_definetopology,VF_DEFINETOPOLOGY)
#define VF_RESETTOPOLOGY    TR_ROUTINE_GLOBAL_(vf_resettopology,VF_RESETTOPOLOGY)
#define VF_CALCHEMICUBE     TR_ROUTINE_GLOBAL_(vf_calchemicube,VF_CALCHEMICUBE)
#define VF_CALCPAIRWISE     TR_ROUTINE_GLOBAL_(vf_calcpairwise,VF_CALCPAIRWISE)
#define VF_JITTERON         TR_ROUTINE_GLOBAL_(vf_jitteron,VF_JITTERON)
#define VF_JITTEROFF        TR_ROUTINE_GLOBAL_(vf_jitteroff,VF_JITTEROFF)
#define VF_SMOOTHMATRIX     TR_ROUTINE_GLOBAL_(vf_smoothmatrix,VF_SMOOTHMATRIX)
#define VF_OUTPUTMATRIXSUMMARYBANNER TR_ROUTINE_GLOBAL_(vf_outputmatrixsummarybanner,VF_OUTPUTMATRIXSUMMARYBANNER)
#define VF_GETROWCOUNTS     TR_ROUTINE_GLOBAL_(vf_getrowcounts,VF_GETROWCOUNTS)
#define VF_GETMATRIX        TR_ROUTINE_GLOBAL_(vf_getmatrix,VF_GETMATRIX)

void VF_SETUP()
{
  VF_Setup(0,1,MPI_COMM_WORLD);
}

void VF_SETNUMENCLOSURES(int *nenclosures)
{
  VF_SetNumEnclosures(*nenclosures);
}

void VF_SETMAXSURFACES(int *max_surfaces)
{
  VF_SetMaxSurfaces(*max_surfaces);
}

void VF_CLEANUP()
{
  VF_CleanUp();
}

int VF_DEFINEENCLOSURE(char *name, int *nonblocking,
		        int *partial, double *asink, 
		        int *npatches, int global_ids[], 
		        int *debug_level, int len)
{
  int encl;
  char *enclID;
  enclID = (char*)malloc((len+1)*sizeof(char));
  strncpy(enclID,name,len);
  enclID[len] = '\0';
  encl = VF_DefineEnclosure(enclID, *nonblocking, *partial, *asink,
                            *npatches, global_ids, *debug_level);
  free(enclID);
  return encl;
}

void VF_RANDOMIZESURFACESON()
{
  VF_RandomizeSurfacesOn();
}

void VF_RANDOMIZESURFACESOFF()
{
  VF_RandomizeSurfacesOff();
}

void VF_DEFINETOPOLOGY(int* enclosure, int* geom_type, 
		       int* nfacets, int* nnodes, 
		       double *x, double *y, double *z, int *c, 
		       int *f2p_map, int *nrotations,
		       int *x_mirror, int *y_mirror, int *z_mirror,
		       int *bsp_depth, int *bsp_length,
		       double *spatial_tol,
		       int *debug_level)
{
  VF_DefineTopology(*enclosure, *geom_type, *nfacets, *nnodes, x, y, z, c, 
                    1, f2p_map, *nrotations, *x_mirror, *y_mirror, *z_mirror,
                    *bsp_depth, *bsp_length, *spatial_tol, *debug_level);
}

void VF_RESETTOPOLOGY(int *enclosure)
{
  VF_ResetTopology(*enclosure);
}

void VF_JITTERON ()
{
  VF_JitterOn();
}

void VF_JITTEROFF ()
{
  VF_JitterOff();
}

void VF_CALCHEMICUBE(int *encl, int *hc_sub_divide, int *hc_resolution, double *hc_min_sep)
{
  VF_CalcHemicube(*encl,*hc_sub_divide,*hc_resolution,*hc_min_sep);
}

void VF_CALCPAIRWISE(int *encl,
                     int *vis_nsamples, int *vis_sampling,
                     int *mc_nsamples, int *mc_sampling,
                     double *mc_tol1,  double *mc_tol2)
{
  VF_CalcPairwise(*encl, *vis_nsamples, *vis_sampling,
                  *mc_nsamples, *mc_sampling, *mc_tol1, *mc_tol2);
}

void VF_SMOOTHMATRIX(int *encl, double *wt, double *tol, int *max_iter,
                     int *symmetric, int *output)
{
  VF_SmoothMatrix(*encl, *wt, *tol, *max_iter, *symmetric, *output);
}

void VF_OUTPUTMATRIXSUMMARYBANNER()
{
  VF_OutputMatrixSummaryBanner();
}

void VF_GETROWCOUNTS(int *encl, int *mode, int count[])
{
  VF_GetRowCounts(*encl, *mode, count);
}

void VF_GETMATRIX(int *encl, int vf_cnt[], int vf_index[], float vf_data[],
                  float vf_diag[], float vf_virt[])
{
  VF_GetMatrix(*encl, vf_cnt, vf_index, vf_data, vf_diag, vf_virt);
}
