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
@(#)    $RCSfile: VF_OutputBanner.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_OutputBanner.c,v $
@(#)
@(#)    DESCRIPTION:  Top-level routine to calculate viewfactors.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>

#include "vf.h"

void VF_OutputInitBanner(void)
{   
  if (VFLIB_Rank==0) {
    printf("\n\n");
    printf("  ********************************");
    printf("**********************************\n");
    printf("     C H A P A R R A L  --  Version %s  --  %s\n",
           VF_VERSION,VF_DATE);
    printf("  ********************************");
    printf("**********************************\n");
    printf("\n");
    printf("   Initializing for:\n");
    printf("      number of processors  = %d\n",VFLIB_Size);
    printf("      number of enclosures  = %d \n",VF_NumEnclosures());
    printf("      max number of patches = %d\n",VF_MaxPatches());
  }
}

void VF_OutputCalcStartBanner(void)
{   
  int         maxDepth, maxLength, minLength, Nleafs;
  BinNodePtr  currentNode;
  BSPstackPtr stack;
  BinTree     *BSPTree;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();

  double dsize = VF_GetMemoryUsage();
  VF_GlobalSumDouble(&dsize);
  if (enclosure->debug_level>=VF_OUTPUT_SUMMARY && VFLIB_Rank==0) {
    if (enclosure->vf_method!=VF_FILE_READ) {
      maxDepth    = 0;
      maxLength   = 0;
      minLength   = topology->nfacets_g+10;
      Nleafs      = 0;
      stack       = topology->BSP_stack;
      BSPTree     = &(topology->BSP_tree);
      currentNode = BSPTree->root;
      VF_BSP_InitStack(stack);
      while (currentNode != NULL) {
        if (currentNode->child[0] == NULL) {
          if (currentNode->members.length>0) Nleafs++;
          maxDepth  = MAX(maxDepth,currentNode->depth);
          maxLength = MAX(maxLength,currentNode->members.length);
          minLength = MIN(minLength,currentNode->members.length);
        } else {
          VF_BSP_PushStack(stack, currentNode->child[0]);
          VF_BSP_PushStack(stack, currentNode->child[1]);
        }
        VF_BSP_PopStack(stack, currentNode);
      }
    }
    printf("\n\n");
    printf("  ********************************");
    printf("**********************************\n");
    printf("   V I E W F A C T O R    C A L C U L A T I O N\n");
    printf("  ********************************");
    printf("**********************************\n");
    printf("\n");
    if (enclosure->vf_method!=VF_FILE_READ) {
      printf("   Calculating viewfactors for enclosure <%s>\n",enclosure->id);
      printf("     enclosure geometry:    ");
      switch (topology->geom) {
      case VF_2Daxisym:
        printf("axisymmetric\n");
        break;
      case VF_2Dplanar:
        printf("planar\n");
        break;
      case VF_3D:
        printf("3D\n");
        break;
      }
    } else {
      printf("   Reading viewfactors for enclosure <%s>\n",enclosure->id);
    }
    printf("     enclosure type:        ");
    if (enclosure->partial) {
      printf("partial (area=%g)",enclosure->asink);
    } else {
      printf("full");
    }
    if (enclosure->nonblocking) {
      printf(", nonblocking\n");
    } else {
      printf(", blocking\n");
    }
    if (enclosure->vf_method!=VF_FILE_READ && topology->geom==1) {
      printf("     # of rotations:        %d\n",topology->nrotations);
    }
    printf("     # of patches:          %d\n",enclosure->npatches_g);
    if (enclosure->vf_method!=VF_FILE_READ) {
      printf("     # of facets:           %d\n",topology->nfacets_g);
      printf("     # of nodes:            %d\n",topology->nnodes);
      printf("     spatial tolerance:     %g\n",topology->spatial_tol);
      printf("     BSP target max depth:  %d\n",topology->bsp_depth);
      printf("     BSP target min length: %d\n",topology->bsp_length);
      printf("     BSP Tree num leafs:    %d\n",Nleafs);
      printf("     BSP Tree max depth:    %d\n",maxDepth);
      printf("     BSP Tree min length:   %d\n",minLength);
      printf("     BSP Tree max length:   %d\n",maxLength);
      printf("     output level:          %d\n",enclosure->debug_level);
    }
    printf("\n");
    printf("   Data segment memory size = %.2fMb\n\n",dsize);
    switch (enclosure->vf_method) {
    case VF_PAIRWISE:
      printf("   Calling VF_Adaptive()...\n");
      printf("     visibility samples           = %d\n",
             enclosure->adaptive.visibility.sampling.n);
      switch (enclosure->adaptive.visibility.sampling.method) {
      case VF_RANDOM_SAMPLE:
        printf("     visibility sampling method   = random\n");
        break;
      case VF_UNIFORM_SAMPLE:
        printf("     visibility sampling method   = uniform\n");
        break;
      case VF_JITTER_SAMPLE:
        printf("     visibility sampling method   = jittered\n");
        break;
      case VF_HALTON_SAMPLE:
        printf("     visibility sampling method   = halton\n");
        break;
      default:
        printf("     visibility sampling method   = unknown\n");
        break; 
      }
      printf("     monte-carlo samples          = %d\n",
             enclosure->adaptive.montecarlo.sampling.n);
      switch (enclosure->adaptive.montecarlo.sampling.method) {
      case VF_RANDOM_SAMPLE:
        printf("     monte-carlo sampling method  = random\n");
        printf("     monte-carlo convergence tol1 = %f\n",enclosure->adaptive.montecarlo.tol1);
        printf("     monte-carlo convergence tol2 = %f\n",enclosure->adaptive.montecarlo.tol2);
        break;
      case VF_UNIFORM_SAMPLE:
        printf("     monte-carlo sampling method  = uniform\n");
        break;
      case VF_JITTER_SAMPLE:
        printf("     monte-carlo sampling method  = jittered\n");
        break;
      case VF_HALTON_SAMPLE:
        printf("     monte-carlo sampling method  = halton\n");
        printf("     monte-carlo convergence tol1 = %f\n",enclosure->adaptive.montecarlo.tol1);
        printf("     monte-carlo convergence tol2 = %f\n",enclosure->adaptive.montecarlo.tol2);
        break;
      default:
        printf("     monte-carlo sampling method  = unknown\n");
        break;
      }
      break;
    case VF_HEMICUBE:
      printf("   Calling VF_Hemicube()...\n");
      printf("     resolution       = %d\n",enclosure->hemicube.resolution);
      printf("     max subdivisions = %d\n",enclosure->hemicube.sub_divide);
      printf("     min seperation   = %g\n",enclosure->hemicube.min_separation);
      break;
    case VF_FILE_READ:
      printf("   Calling VF_Read()...");
    }
    printf("\n");
    fflush(stdout);
  }
}

void VF_OutputCalcEndBanner(void)
{
  VFenclosure *enclosure=VF_CurrentEnclosure();

  if (enclosure->debug_level>0) {
    if (enclosure->vf_method==VF_HEMICUBE) {
      VF_MinDouble(&(enclosure->hemicube.min_distance),0);
      VF_MaxInt(&(enclosure->hemicube.imax),0);
      VF_MaxInt(&(enclosure->hemicube.jmax),0);
    }
    double dsize = VF_GetMemoryUsage();
    VF_GlobalSumDouble(&dsize);
    if (VFLIB_Rank==0) {
      if (enclosure->vf_method==VF_HEMICUBE) {
        printf("     Minimum effective surface radius    = %g\n",
               enclosure->hemicube.min_equiv_radius);
        printf("     Minimum seperation distance         = %g\n",
               enclosure->hemicube.min_distance);
        printf("     Maximum desired surface subdivision = %d, %d\n",
               enclosure->hemicube.imax,enclosure->hemicube.jmax);
        printf("     Actual maximum surface subdivision  = %d\n",
               enclosure->hemicube.sub_divide);
      } else if (enclosure->vf_method==VF_PAIRWISE) {
        printf("     N total      = %d\n",enclosure->adaptive.num_total);
        printf("     N zero       = %d\n",enclosure->adaptive.num_zero);
        printf("     N analytic   = %d\n",enclosure->adaptive.num_analytic);
        printf("     N numerical  = %d\n",enclosure->adaptive.num_numerical);
      }
      printf("\n");
      printf("   Data segment memory size = %.2fMb\n",dsize);
      printf("\n");
      printf("   Elapsed time = %.2f\n",enclosure->time_vf_calc);
      printf("\n\n");
      fflush(stdout);
    }
  }
}

void VF_OutputSmoothingBanner(double wt, double tol, int max_iter)
{   
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();

  if (enclosure->debug_level>=VF_OUTPUT_SUMMARY && VFLIB_Rank==0) {
    printf("\n\n");
    printf("  ********************************");
    printf("**********************************\n");
    printf("   V I E W F A C T O R    M A T R I X    S M O O T H I N G\n");
    printf("  ********************************");
    printf("**********************************\n");
    printf("\n");
    printf("   Smoothing of viewfactor matrix for enclosure <%s>\n",
           enclosure->id);
    printf("     PCG Solver, wt = %.1f\n",wt);
    printf("     Max iterations = %d\n",max_iter);
    printf("     Tolerance      = %g\n",tol);
    printf("\n");
    fflush(stdout);
  }
}

void VF_OutputMatrixSummaryBanner(void)
{
  int         i, debug_level;
  double      density;
  double      min0, max0, sum0, sd0, tsum0;
  double      min1, max1, sum1, sd1, tsum1;
  double      min2, max2, sum2, sd2, tsum2;
  VFrow       *row;
  VFenclosure *enclosure=VF_CurrentEnclosure();

  debug_level = enclosure->debug_level;
  if (debug_level>0) {
    if (VFLIB_Rank==0) {
      printf("\n\n");
      printf("  ********************************");
      printf("**********************************\n");
      printf("   V I E W F A C T O R    M A T R I X    S U M M A R Y\n");
      printf("  ********************************");
      printf("**********************************\n");
      printf("\n");
      printf("   Viewfactor matrix summary for enclosure <%s>\n",enclosure->id);
      printf("\n");
    }
    /*===========================================*/
    /* FIND THE DENSITY OF THE VIEWFACTOR MATRIX */
    /*===========================================*/
    for (density=0.0, i=0; i<enclosure->npatches_l; i++) {
      density += (double)enclosure->row[i].array0.cnt;
      density += (double)enclosure->row[i].array1.cnt;
    }
    VF_GlobalSumDouble(&density);
    density /= (double)enclosure->npatches_g*(double)enclosure->npatches_g;
    VF_RowsumStats(&min0, &max0, &sum0, &sd0, 
                   &min1, &max1, &sum1, &sd1,
                   &min2, &max2, &sum2, &sd2);
    row   = enclosure->row;
    tsum0 = tsum1 = tsum2 =0.0;
    for (i=0; i<enclosure->npatches_l; i++, row++) {
      tsum0 += row->raw_rowsum;
      tsum1 += row->sym_rowsum;
      tsum2 += row->smooth_rowsum;
    } 
    VF_GlobalSumDouble(&tsum0);
    VF_GlobalSumDouble(&tsum1);
    VF_GlobalSumDouble(&tsum2);
    if (enclosure->vf_method==VF_HEMICUBE) {
      VF_MinDouble(&(enclosure->hemicube.min_distance),0);
      VF_MinInt(&(enclosure->hemicube.imax),0);
      VF_MinInt(&(enclosure->hemicube.jmax),0);
    }
    double dsize = VF_GetMemoryUsage();
    VF_GlobalSumDouble(&dsize);
    if (VFLIB_Rank==0) {
      printf("   Target rowsum                  = %d\n",enclosure->npatches_g);
      printf("   Raw rowsum total               = %f\n",tsum0);
      printf("   Raw rowsum error min           = %11.4e\n",min0);
      printf("   Raw rowsum error max           = %11.4e\n",max0);
      printf("   Raw rowsum error mean          = %11.4e +/- %10.4e\n",
             sum0,sd0);
      if (enclosure->smoothed) {
        printf("   Smoothed rowsum total          = %f\n",tsum2);
        printf("   Smoothed rowsum error min      = %11.4e\n",min2);
        printf("   Smoothed rowsum error max      = %11.4e\n",max2);
        printf("   Smoothed rowsum error mean     = %11.4e +/- %10.4e\n",sum2,sd2);
      }
      printf("\n");
      printf("   Viewfactor matrix is %.2f%% dense\n\n",density*100.0);
      printf("   Total elapsed time = %.2f (sec)\n",enclosure->time_vf_init+enclosure->time_vf_calc+
                                                    enclosure->time_vf_symmetry+enclosure->time_vf_smooth);
      printf("     (Initialization) = %.2f\n",enclosure->time_vf_init);
      printf("     (Calculation)    = %.2f\n",enclosure->time_vf_calc);
      printf("     (Smoothing)      = %.2f\n",enclosure->time_vf_symmetry+enclosure->time_vf_smooth);
      printf("\n");
      printf("   Data segment memory size  = %.2fMb\n",dsize);
      printf("\n\n");
      fflush(stdout);
    }
  }
}
