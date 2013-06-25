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
@(#)    $RCSfile: VF_SamplePoly.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_SamplePoly.c,v $
@(#)
@(#)    DESCRIPTION:  Implement various polygon sampling methods for the
@(#)    visibility check, Classical MonteCarlo and Quasi-MonteCarlo
@(#)    integration routines.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "vf.h"

double VF_Halton(int i, int base);

void VF_SamplePoly(Poly *poly, Sampling *sampling, Point samples[])
{
  switch (sampling->method) {
  case VF_RANDOM_SAMPLE:
    VF_SampleRandom(poly, sampling, samples);
    break;
  case VF_UNIFORM_SAMPLE:
    VF_SampleUniform(poly, sampling, samples);
    break;
  case VF_JITTER_SAMPLE:
    VF_SampleJitter(poly, sampling, samples);
    break;
  case VF_HALTON_SAMPLE:
    VF_SampleHalton(poly, sampling, samples);
    break;
  }
}

void VF_SampleRandom(Poly *poly, Sampling *sampling, Point samples[])
{
  int    i;
  double u, v, t=1.0e-3;
  VFenclosure *encl=VF_CurrentEnclosure();
    
  for (i=0; i<sampling->n; i++) {
    u          = RAND(0.0+t,1.0-t);
    v          = RAND(0.0+t,1.0-t);
    samples[i] = VF_UVtoXYZ(u, v, poly);
  }
}

void VF_SampleUniform(Poly *poly, Sampling *sampling, Point samples[])
{
  int    i, j, k, m, n;
  double width, L2, L3, L2sum, L3sum;
  
  switch (poly->np) {
  case VF_LINE:
    m     = sampling->n;
    width = 1.0/(double)(m);
    /*===============================*/
    /* INITIALIZE TO UNIFORM PATTERN */
    /*===============================*/
    for (i=0; i<m; i++) {
      samples[i].x = i*width+0.5*width;
      samples[i].y = 0.0;
    }
    break;
  case VF_TRI:
    m     = (int)(sqrt((double)(sampling->n))+0.5);
    width = 1.0/(double)(m);
    /*===============================*/
    /* INITIALIZE TO UNIFORM PATTERN */
    /*===============================*/
    for (k=0, i=0; i<m; i++) {
      n  = i*2+1;
      L2 = i*width;
      L3 = 0.0;
      for (j=0; j<n; j++) {
        if (j%2) {
          L2sum = 2.0*L2+(L2-width);
          L3sum = 2.0*L3+(L3-width);
        } else {
          L2sum = 2.0*L2+(L2+width);
          L3sum = 2.0*L3+(L3+width);
        }
        L2sum       /= 3.0;
        L3sum       /= 3.0;
        samples[k].x = L2sum;
        samples[k].y = L3sum;
        if (j%2) {
          L2 -= width;
        } else {
          L3 += width;
        }
        k++;
      }
    }
    break;
  case VF_QUAD:
    m     = (int)(sqrt((double)(sampling->n))+0.5);
    width = 1.0/(double)(m);
    /*===============================*/
    /* INITIALIZE TO UNIFORM PATTERN */
    /*===============================*/
    for (i=0; i<m; i++) {
      for (j=0; j<m; j++) {
        samples[i*m+j].x = i*width+0.5*width;
        samples[i*m+j].y = j*width+0.5*width;
      }
    }
    break;
  default:
    printf("VF_SampleUniform():  poly->np = %d\n",poly->np);
    VF_Exit(1);
    break;
  }
  VF_ShuffleSamples(sampling->n, samples);
  VF_ConvertUVtoXYZ(poly, sampling->n, samples);
}

void VF_SampleJitter(Poly *poly, Sampling *sampling, Point samples[])
{
  int    i, j, k, m, n;
  double width, L2, L3, L2sum, L3sum;
  VFenclosure *encl=VF_CurrentEnclosure();
    
  switch (poly->np) {
  case VF_LINE:
    m     = sampling->n;
    width = 1.0/(double)(m);
    /*=========================================*/
    /* INITIALIZE TO STRATIFIED JITTER PATTERN */
    /*=========================================*/
    for (i=0; i<m; i++) {
      samples[i].x = i*width+RAND(0.0,width);
      samples[i].y = 0.0;
    }
    break;
  case VF_TRI:
    m     = (int)(sqrt((double)(sampling->n))+0.5);
    width = 1.0/(double)(m);
    /*===============================*/
    /* INITIALIZE TO UNIFORM PATTERN */
    /*===============================*/
    for (k=0, i=0; i<m; i++) {
      n  = i*2+1;
      L2 = i*width;
      L3 = 0.0;
      for (j=0; j<n; j++) {
        if (j%2) {
          L2sum = 2.0*L2+(L2-width);
          L3sum = 2.0*L3+(L3-width);
        } else {
          L2sum = 2.0*L2+(L2+width);
          L3sum = 2.0*L3+(L3+width);
        }
        L2sum       /= 3.0;
        L3sum       /= 3.0;
        L2sum       += RAND(0.0,width*2.0/3.0)-width/3.0;
        L3sum       += RAND(0.0,width*2.0/3.0)-width/3.0;
        samples[k].x = L2sum;
        samples[k].y = L3sum;
        if (j%2) {
          L2 -= width;
        } else {
          L3 += width;
        }
        k++;
      }
    }
    break;
  case VF_QUAD:
    m     = (int)(sqrt((double)(sampling->n))+0.5);
    width = 1.0/(double)(m);
    /*=========================================*/
    /* INITIALIZE TO STRATIFIED JITTER PATTERN */
    /*=========================================*/
    for (i=0; i<m; i++) {
      for (j=0; j<m; j++) {
        samples[i*m+j].x = i*width+RAND(0.0,width);
        samples[i*m+j].y = j*width+RAND(0.0,width);
      }
    }
    break;
  default:
    printf("VF_SampleJitter():  poly->np = %d\n",poly->np);
    VF_Exit(1);
    break;
  }
  VF_ShuffleSamples(sampling->n, samples);
  VF_ConvertUVtoXYZ(poly, sampling->n, samples);
}

void VF_SampleHalton(Poly *poly, Sampling *sampling, Point samples[])
{
  int i;
    
  for (i=0; i<sampling->n; i++) {
    /*
    samples[i].x = VF_Halton3(sampling->index);
    samples[i].y = VF_Halton5(sampling->index);
    */
    samples[i].x = VF_Halton(sampling->index,2);
    samples[i].y = VF_Halton(sampling->index,3);
    sampling->index++;
  }
  VF_ShuffleSamples(sampling->n, samples);
  VF_ConvertUVtoXYZ(poly, sampling->n, samples);
}

void VF_SetupSampling(Sampling *sampling, Point uv_samples[])
{
  switch (sampling->method) {
  case VF_RANDOM_SAMPLE:
    VF_SetupSampleRandom(sampling, uv_samples);
    break;
  case VF_UNIFORM_SAMPLE:
    VF_SetupSampleUniform(sampling, uv_samples);
    break;
  case VF_JITTER_SAMPLE:
    VF_SetupSampleJitter(sampling, uv_samples);
    break;
  case VF_HALTON_SAMPLE:
    VF_SetupSampleHalton(sampling, uv_samples);
    break;
  }
}

void VF_SetupSampleRandom(Sampling *sampling, Point uv_samples[])
{
  int    i;
  double t=1.0e-3;
  VFenclosure *encl=VF_CurrentEnclosure();
    
  for (i=0; i<sampling->n; i++) {
    uv_samples[i].x = RAND(0.0+t,1.0-t);
    uv_samples[i].y = RAND(0.0+t,1.0-t);
  }
}

void VF_SetupSampleUniform(Sampling *sampling, Point uv_samples[])
{
  int    i, j, m, n, nq, nt;
  double width, L2, L3, L2sum, L3sum;
  VFtopology *topology=VF_CurrentTopology();;
    
  switch (topology->geom) {
  case VF_2Dplanar:
    m     = sampling->n;
    width = 1.0/(double)(m);
    /*===============================*/
    /* INITIALIZE TO UNIFORM PATTERN */
    /*===============================*/
    for (i=0; i<m; i++) {
      uv_samples[i].x = i*width+0.5*width;
    }
    break;
  default:
    m     = (int)(sqrt((double)(sampling->n))+0.5);
    nq    = m*m;
    width = 1.0/(double)(m);
    /*===============================*/
    /* INITIALIZE TO UNIFORM PATTERN */
    /*===============================*/
    for (i=0; i<m; i++) {
      for (j=0; j<m; j++) {
        uv_samples[i*m+j].x = i*width+0.5*width;
        uv_samples[i*m+j].y = j*width+0.5*width;
      }
    }
        
    m     = (int)(sqrt((double)(sampling->n))+0.5);
    width = 1.0/(double)(m);
    /*===============================*/
    /* INITIALIZE TO UNIFORM PATTERN */
    /*===============================*/
    for (nt=0, i=0; i<m; i++) {
      n  = i*2+1;
      L2 = i*width;
      L3 = 0.0;
      for (j=0; j<n; j++) {
        if (j%2) {
          L2sum = 2.0*L2+(L2-width);
          L3sum = 2.0*L3+(L3-width);
        } else {
          L2sum = 2.0*L2+(L2+width);
          L3sum = 2.0*L3+(L3+width);
        }
        L2sum              /= 3.0;
        L3sum              /= 3.0;
        uv_samples[nt+nq].x = L2sum;
        uv_samples[nt+nq].y = L3sum;
        if (j%2) {
          L2 -= width;
        } else {
          L3 += width;
        }
        nt++;
      }
    }
    break;
  }
}

void VF_SetupSampleJitter(Sampling *sampling, Point uv_samples[])
{
  int    i, j, m, n, nq, nt;
  double width, L2, L3, L2sum, L3sum;
  VFtopology *topology=VF_CurrentTopology();
  VFenclosure *encl=VF_CurrentEnclosure();
    
  switch (topology->geom) {
  case VF_2Dplanar:
    m     = sampling->n;
    width = 1.0/(double)(m);
    /*=========================================*/
    /* INITIALIZE TO STRATIFIED JITTER PATTERN */
    /*=========================================*/
    for (i=0; i<m; i++) {
      uv_samples[i].x = i*width+RAND(0.0,width);
    }
    break;
  default:
    m     = (int)(sqrt((double)(sampling->n))+0.5);
    nq    = m*m;
    width = 1.0/(double)(m);
    /*=========================================*/
    /* INITIALIZE TO STRATIFIED JITTER PATTERN */
    /*=========================================*/
    for (i=0; i<m; i++) {
      for (j=0; j<m; j++) {
        uv_samples[i*m+j].x = i*width+RAND(0.0,width);
        uv_samples[i*m+j].y = j*width+RAND(0.0,width);
      }
    }
    
    m     = (int)(sqrt((double)(sampling->n))+0.5);
    width = 1.0/(double)(m);
    /*===============================*/
    /* INITIALIZE TO UNIFORM PATTERN */
    /*===============================*/
    for (nt=0, i=0; i<m; i++) {
      n  = i*2+1;
      L2 = i*width;
      L3 = 0.0;
      for (j=0; j<n; j++) {
        if (j%2) {
          L2sum = 2.0*L2+(L2-width);
          L3sum = 2.0*L3+(L3-width);
        } else {
          L2sum = 2.0*L2+(L2+width);
          L3sum = 2.0*L3+(L3+width);
        }
        L2sum              /= 3.0;
        L3sum              /= 3.0;
        L2sum              += RAND(0.0,width*2.0/3.0)-width/3.0;
        L3sum              += RAND(0.0,width*2.0/3.0)-width/3.0;
        uv_samples[nt+nq].x = L2sum;
        uv_samples[nt+nq].y = L3sum;
        if (j%2) {
          L2 -= width;
        } else {
          L3 += width;
        }
        nt++;
      }
    }
    break;
  }
}

void VF_SetupSampleHalton(Sampling *sampling, Point uv_samples[])
{
  int    i;
    
  for (i=0; i<sampling->n; i++) {
    uv_samples[i].x = VF_Halton(sampling->index,2);
    uv_samples[i].y = VF_Halton(sampling->index,3);
    sampling->index++;
  }
}

void VF_ShuffleSamples(int ns, Point samples[])
{
  int   i, j;
  Point tmp;
  VFenclosure *encl=VF_CurrentEnclosure();
    
  for (i=ns-1; i>1; i--) {
    j          = RANI(0,i-1);
    tmp        = samples[j];
    samples[j] = samples[i];
    samples[i] = tmp;
  }
}

void VF_SamplesUVtoXYZ(Poly *poly, Sampling *sampling, 
                       Point uv_samples[], Point xyz_samples[])
{
  int    i;
  double u, v;
  Point  *samples;
    
  samples = uv_samples;
  if ((sampling->method==VF_UNIFORM_SAMPLE || 
       sampling->method==VF_JITTER_SAMPLE) && 
       poly->np==3) samples = &uv_samples[sampling->n];
  VF_ShuffleSamples(sampling->n, samples);
  for (i=0; i<sampling->n; i++) {
    u              = samples[i].x;
    v              = samples[i].y;
    xyz_samples[i] = VF_UVtoXYZ(u, v, poly);
  }
}

void VF_ConvertUVtoXYZ(Poly *poly, int ns, Point samples[])
{
  int    i;
  double u, v;
    
  for (i=0; i<ns; i++) {
    u          = samples[i].x;
    v          = samples[i].y;
    samples[i] = VF_UVtoXYZ(u, v, poly);
  }
}

double VF_Halton(int i, int base)
{
  double h=0.0, f, factor;
    
  f = factor = 1.0/(double)base;
  while (i>0) {
    h      += (double)(i%base) * factor;
    i      /= base;
    factor *= f;
  }
  return h;
}
