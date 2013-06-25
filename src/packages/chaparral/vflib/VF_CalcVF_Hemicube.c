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
@(#)    $RCSfile: VF_CalcVF_Hemicube.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcVF_Hemicube.c,v $
@(#)
@(#)    DESCRIPTION:  Use monte-carlo method to numerically 
@(#)    evaluate the viewfactor between two surfaces.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vf.h"

double VF_CalcVF_Hemicube (Poly *poly_i, Poly *poly_j, 
                           CandidateList *candidates, 
                           int ncandidates)
{
  int         i, is, j, js, k, nnn, *iptr;
  int         hc_size, top_size, side_size;
  int         face, hc_res, hc_res2, iseg=1, jseg=1;
  float       *delta_ptr;
  double      phi, phi1, phi2, r, patch_area=0.0, poly_area;
  double      sub_area, new_area, new_diameter, eff_diameter;
  double      *fptr, ff=0.0, vf, aspect;
  Plane       view_plane;
  Point       center_i, center_j;
  Vector      Phi1, Phi2, Phi3, Phi4;
  Poly        *poly_jj, poly;
  ViewPort    view;
  VFenclosure *enclosure=VF_CurrentEnclosure();
  Hemicube    *hemicube=&(enclosure->hemicube);

  poly.np     = poly_i->np;
  poly.d      = poly_i->d;
  poly.normal = poly_i->normal;
  hc_res      = hemicube->resolution;
  hc_res2     = hc_res/2;
  top_size    = hc_res*hc_res;
  side_size   = hc_res*hc_res2;
  VF_PolyCenter(poly_i, &center_i);
  VF_PolyCenter(poly_j, &center_j);
  eff_diameter = 2.0*sqrt(VF_PolyArea(poly_i)/M_PI);
  r            = V3_DistanceBetween2Points(&center_i,&center_j);
  if (r<eff_diameter*hemicube->min_separation) {
    if (poly_i->np == 4) {
      new_diameter = r/hemicube->min_separation;
      aspect       = VF_QuadAspectRatio(poly_i);
      if (aspect>=1.0) {
        sub_area = VF_PolyArea(poly_i)/aspect;
        new_area = M_PI*0.25*new_diameter*new_diameter;
        jseg     = (int)(sqrt((double)(sub_area/new_area))+0.5);
        iseg     = (int)((double)jseg*aspect);
      } else {
        aspect   = 1.0/aspect;
        sub_area = VF_PolyArea(poly_i)/aspect;
        new_area = M_PI*0.25*new_diameter*new_diameter;
        iseg     = (int)(sqrt((double)(sub_area/new_area))+0.5);
        jseg     = (int)((double)iseg*aspect);
      }
      iseg = MAX(1,iseg);
      jseg = MAX(1,jseg);
      iseg = MIN(hemicube->sub_divide,iseg);
      jseg = MIN(hemicube->sub_divide,jseg);
    } else {
      new_diameter = r/hemicube->min_separation;
      iseg         = (int)(eff_diameter/new_diameter+0.5);
      iseg         = MIN(hemicube->sub_divide,iseg);
      jseg         = iseg;
    }
  }
  if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
    printf("           (iseg, jseg) = (%d, %d)\n",iseg,jseg);
  }
  for (is=0; is<iseg; is++) {
    if (poly_i->np == 4) {
      nnn = jseg;
    } else {
      nnn         = 2*(is+1)-1;
      phi         = (double)(is)/(double)(iseg);
      Phi1.x      = (poly_i->p[1].x-poly_i->p[0].x)*
                    phi+poly_i->p[0].x;
      Phi1.y      = (poly_i->p[1].y-poly_i->p[0].y)*
                    phi+poly_i->p[0].y;
      Phi1.z      = (poly_i->p[1].z-poly_i->p[0].z)*
                    phi+poly_i->p[0].z;
      Phi2.x      = (poly_i->p[2].x-poly_i->p[0].x)*
                    phi+poly_i->p[0].x;
      Phi2.y      = (poly_i->p[2].y-poly_i->p[0].y)*
                    phi+poly_i->p[0].y;
      Phi2.z      = (poly_i->p[2].z-poly_i->p[0].z)*
                    phi+poly_i->p[0].z;
      phi         = (double)(is+1)/(double)(iseg);
      Phi3.x      = (poly_i->p[1].x-poly_i->p[0].x)*
                    phi+poly_i->p[0].x;
      Phi3.y      = (poly_i->p[1].y-poly_i->p[0].y)*
                    phi+poly_i->p[0].y;
      Phi3.z      = (poly_i->p[1].z-poly_i->p[0].z)*
                    phi+poly_i->p[0].z;
      Phi4.x      = (poly_i->p[2].x-poly_i->p[0].x)*
                    phi+poly_i->p[0].x;
      Phi4.y      = (poly_i->p[2].y-poly_i->p[0].y)*
                    phi+poly_i->p[0].y;
      Phi4.z      = (poly_i->p[2].z-poly_i->p[0].z)*
                    phi+poly_i->p[0].z;
      poly.p[2].x = Phi1.x;
      poly.p[2].y = Phi1.y;
      poly.p[2].z = Phi1.z;
      poly.p[1].x = Phi3.x;
      poly.p[1].y = Phi3.y;
      poly.p[1].z = Phi3.z;
    }
    for (js=0; js<nnn; js++) {
      if (poly_i->np == 4) {
        phi1 = (double)(is)/(double)(iseg);
        phi2 = (double)(js)/(double)(jseg);     
        VF_UV_to_XYZ(poly_i, phi1, phi2, &(poly.p[0]));
        phi1 = (double)(is+1)/(double)(iseg);
        phi2 = (double)(js)/(double)(jseg);     
        VF_UV_to_XYZ(poly_i, phi1, phi2, &(poly.p[1]));
        phi1 = (double)(is+1)/(double)(iseg);
        phi2 = (double)(js+1)/(double)(jseg);     
        VF_UV_to_XYZ(poly_i, phi1, phi2, &(poly.p[2]));
        phi1 = (double)(is)/(double)(iseg);
        phi2 = (double)(js+1)/(double)(jseg);     
        VF_UV_to_XYZ(poly_i, phi1, phi2, &(poly.p[3]));
      } else {
        phi1 = (double)(js)/(double)(nnn);
        phi2 = (double)(js)/(double)(iseg);
        if (js%2 == 0) {
          phi         = (double)(1+js/2)/(double)(is+1);
          poly.p[0]   = poly.p[2];
          poly.p[2].x = (Phi4.x-Phi3.x)*phi+Phi3.x;
          poly.p[2].y = (Phi4.y-Phi3.y)*phi+Phi3.y;
          poly.p[2].z = (Phi4.z-Phi3.z)*phi+Phi3.z;
        } else {
          phi         = (double)((1+js)/2)/(double)(is);
          poly.p[1]   = poly.p[2];
          poly.p[2].x = (Phi2.x-Phi1.x)*phi+Phi1.x;
          poly.p[2].y = (Phi2.y-Phi1.y)*phi+Phi1.y;
          poly.p[2].z = (Phi2.z-Phi1.z)*phi+Phi1.z;
        }
      }
      hc_size       = top_size;
      delta_ptr     = hemicube->top_deltaVF;
      hemicube->dir = hemicube->top_dir;
      VF_SetView(&view, &poly);
      view.window.x1 = hc_res-1;
      view.window.y1 = hc_res-1;
      for (vf=0.0, face=0; face<5; face++) {
        /*==================================*/
        /* SET THE VIEWPORT FOR THE CURRENT */ 
        /* FACE AND DEFINE THE VIEW PLANE   */
        /*==================================*/
        VF_SetViewPort(&view, face);
        view_plane.normal = view.view_normal;
        view_plane.d      = -V3_Dot(&(view_plane.normal), &(view.view_point));
        /*===================================*/
        /* REINITIALIZE THE HEMICUBE BUFFERS */
        /*===================================*/
        iptr = hemicube->ibuffer;
        fptr = hemicube->zbuffer;
        for (i=0; i<hc_size; i++) {
          *iptr++ = -1;
          *fptr++ = 0.0;
        }
        /*=======================*/
        /* TRAVERSE THE BSP TREE */
        /*=======================*/
        VF_ProjectOntoHemicube(&poly, poly_j, 10, &view,
                               VF_FRUSTUM_CLIP_PARTIAL, 
                               face, hemicube);
        for (k=0; k<ncandidates; k++) {
          poly_jj = &(candidates[k].poly);
          VF_ProjectOntoHemicube(&poly, poly_jj, 0, &view,
                                 VF_FRUSTUM_CLIP_PARTIAL, 
                                 face, hemicube);
        }
        iptr = hemicube->ibuffer;
        for (i=0; i<hc_size; i++) {
          j = *iptr++;
          if (j>1) vf += *delta_ptr;
          delta_ptr++;
        }
        view.window.y1 = hc_res2-1;
        hc_size        = side_size;
        delta_ptr      = hemicube->side_deltaVF;
        hemicube->dir  = hemicube->side_dir;
      }
      poly_area   = VF_PolyArea(&poly);
      patch_area += poly_area;
      ff += poly_area*vf;
    }
  }
  ff /= patch_area;
  return (ff);
}

double VF_CalcVF_Hemicube0 (Poly *poly_i, Poly *poly_j)
{
  int         i, j, *iptr;
  int         hc_size, top_size, side_size;
  int         face, hc_res, hc_res2;
  float       *delta_ptr;
  double      *fptr, vf=0.0;
  Plane       view_plane;
  ViewPort    view;
  VFenclosure *enclosure=VF_CurrentEnclosure();
  Hemicube    *hemicube=&(enclosure->hemicube);

  hc_res        = hemicube->resolution;
  hc_res2       = hc_res/2;
  top_size      = hc_res*hc_res;
  side_size     = hc_res*hc_res2;
  hc_size	= top_size;
  delta_ptr	= hemicube->top_deltaVF;
  hemicube->dir = hemicube->top_dir;
  VF_SetView(&view, poly_i);
  view.window.x1 = hc_res-1;
  view.window.y1 = hc_res-1;
  for (face=0; face<5; face++) {
    /*==================================*/
    /* SET THE VIEWPORT FOR THE CURRENT */ 
    /* FACE AND DEFINE THE VIEW PLANE	*/
    /*==================================*/
    VF_SetViewPort(&view, face);
    view_plane.normal = view.view_normal;
    view_plane.d      = -V3_Dot(&(view_plane.normal), &(view.view_point));
    /*===================================*/
    /* REINITIALIZE THE HEMICUBE BUFFERS */
    /*===================================*/
    iptr = hemicube->ibuffer;
    fptr = hemicube->zbuffer;
    for (i=0; i<hc_size; i++) {
      *iptr++ = -1;
      *fptr++ = 0.0;
    }
    VF_ProjectOntoHemicube(poly_i, poly_j, 1, &view,
        		   VF_FRUSTUM_CLIP_PARTIAL, 
        		   face, hemicube);
    iptr = hemicube->ibuffer;
    for (i=0; i<hc_size; i++) {
      j = *iptr++;
      if (j>0) vf += *delta_ptr;
      delta_ptr++;
    }
    view.window.y1 = hc_res2-1;
    hc_size	   = side_size;
    delta_ptr	   = hemicube->side_deltaVF;
    hemicube->dir  = hemicube->side_dir;
  }
  return (vf);
}
