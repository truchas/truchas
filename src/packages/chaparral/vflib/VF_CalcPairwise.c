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
@(#)    $RCSfile: VF_CalcPairwise.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcPairwise.c,v $
@(#)
@(#)    DESCRIPTION:  Top-level routine to calculate viewfactors.
@(#)
@(#)    encl          enclosure number
@(#)    vis_nsamples  number of samples to use for visibility checks
@(#)    vis_sampling  type of sampling to use for visibility checks
@(#)    mc_nsamples   number of samples to use for QMC calculations
@(#)    mc_sampling   type of sampling to use for QMC calculations
@(#)    mc_tol1       convergence tol1 to use for QMC calculations
@(#)    mc_tol2       convergence tol2 to use for QMC calculations
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vf.h"

void VF_CalcPairwise (int encl, 
                      int vis_nsamples, int vis_sampling, 
                      int mc_nsamples, int mc_sampling, 
                      double mc_tol1, double mc_tol2)
{
  int         ns, vis_ns1, vis_ns2, mc_ns1, mc_ns2;
  double      clock0, clock1;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_GetEnclosure(encl);
    
  enclosure->vf_method = VF_PAIRWISE;
  if (!topology->nonblocking) {
    if (vis_sampling==VF_UNIFORM_SAMPLE || vis_sampling==VF_JITTER_SAMPLE) {
      ns      = (int)(sqrt((double)vis_nsamples)+0.5);
      vis_ns1 = ns*ns;
      vis_ns2 = 2*vis_ns1;
    } else {
      vis_ns1 = vis_nsamples;
      vis_ns2 = vis_nsamples;
    }
    if (mc_sampling==VF_UNIFORM_SAMPLE || mc_sampling==VF_JITTER_SAMPLE) {
      ns     = (int)(sqrt((double)mc_nsamples)+0.5);
      mc_ns1 = ns*ns;
      mc_ns2 = 2*mc_ns1;
    } else {
      mc_ns1 = mc_nsamples;
      mc_ns2 = mc_nsamples;
    }
    enclosure->adaptive.visibility.sampling.method = vis_sampling;
    enclosure->adaptive.visibility.sampling.index  = 11;
    enclosure->adaptive.visibility.sampling.n      = vis_ns1;
    enclosure->adaptive.visibility.ray             = (Ray   *)VF_Newv(vis_ns2*sizeof(Ray));
    enclosure->adaptive.visibility.sampleuv_i      = (Point *)VF_Newv(vis_ns2*sizeof(Point));
    enclosure->adaptive.visibility.sampleuv_j      = (Point *)VF_Newv(vis_ns2*sizeof(Point));
    enclosure->adaptive.montecarlo.sampling.method = mc_sampling;
    enclosure->adaptive.montecarlo.sampling.index  = 11;
    enclosure->adaptive.montecarlo.sampling.n      = mc_ns1;
    enclosure->adaptive.montecarlo.tol1            = mc_tol1;
    enclosure->adaptive.montecarlo.tol2            = mc_tol2;
    enclosure->adaptive.montecarlo.uv_samples      = (Point *)VF_Newv(mc_ns2*sizeof(Point));
    VF_SetupSampling(&(enclosure->adaptive.montecarlo.sampling),enclosure->adaptive.montecarlo.uv_samples);
  } else {
    if (vis_sampling==VF_UNIFORM_SAMPLE || vis_sampling==VF_JITTER_SAMPLE) {
      ns      = (int)(sqrt((double)vis_nsamples)+0.5);
      vis_ns1 = ns*ns;
      vis_ns2 = 2*vis_ns1;
    } else {
      vis_ns1 = vis_nsamples;
      vis_ns2 = vis_nsamples;
    }
    if (mc_sampling==VF_UNIFORM_SAMPLE || mc_sampling==VF_JITTER_SAMPLE) {
      ns     = (int)(sqrt((double)mc_nsamples)+0.5);
      mc_ns1 = ns*ns;
      mc_ns2 = 2*mc_ns1;
    } else {
      mc_ns1 = mc_nsamples;
      mc_ns2 = mc_nsamples;
    }
    enclosure->adaptive.visibility.sampling.method = vis_sampling;
    enclosure->adaptive.visibility.sampling.n      = vis_ns1;
    enclosure->adaptive.montecarlo.sampling.method = mc_sampling;
    enclosure->adaptive.montecarlo.sampling.n      = mc_ns1;
    enclosure->adaptive.montecarlo.tol1            = mc_tol1;
    enclosure->adaptive.montecarlo.tol2            = mc_tol2;
  }
  VF_OutputCalcStartBanner();
  clock0 = VF_Clock();
  VF_MatrixCalcPairwise();
  clock1 = VF_Clock();
  if (!topology->nonblocking) {
    VF_Free(enclosure->adaptive.visibility.ray);
    VF_Free(enclosure->adaptive.visibility.sampleuv_i);
    VF_Free(enclosure->adaptive.visibility.sampleuv_j);
    VF_Free(enclosure->adaptive.montecarlo.uv_samples);
  }
  enclosure->time_vf_calc = clock1-clock0;
  VF_OutputCalcEndBanner();
}
