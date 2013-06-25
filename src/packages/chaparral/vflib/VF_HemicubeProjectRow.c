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
@(#)    $RCSfile: VF_HemicubeProjectRow.c,v $
@(#)    $Revision: 1.5 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_HemicubeProjectRow.c,v $
@(#)
@(#)    DESCRIPTION:  Use hemicube method to calculate one row of
@(#)    the viewfactor matrix.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

#define DO_CLIP_TO_FRUSTUM  1
#define DO_FRONT_TO_BACK    0

#define NO_PPM_DEBUG
#ifdef PPM_DEBUG
void WritePPM(int dx, int dy, int raster[], int maxval);
#endif

void VF_HemicubeProjectRow(int facet_i, int subfacet_i, 
                           Poly *poly_i, ViewPort *view, 
                           int scale_it, double area_i, 
                           double VF[], double vf[])
{
  int         i, j, *iptr, nfaces=5;
  int         hc_size, top_size, side_size;
  int         face, hc_res, hc_res2;
  int         patch_j, facet_j, subfacet_j;
  int         inview, nsub_polys_j, mholes=0;
  float       *delta_ptr;
  double      *fptr, dist, *VFptr;
  Plane       view_plane;
  Point       center;
  Poly        *poly_j;
  POLYstack   PolyStack;
  BinNodePtr  currentNode;
  BSPstackPtr stack;
  Facet       *facet;
  VFtopology  *topology=VF_CurrentTopology();
  VFenclosure *enclosure=VF_CurrentEnclosure();

  if (topology->geom==VF_2Dplanar) nfaces=3;
  VF_InitPolyStack(&PolyStack);
  VF_PolyCenter(poly_i,&center);
  stack     = topology->BSP_stack;
  hc_res    = enclosure->hemicube.resolution;
  hc_res2   = hc_res/2;
  top_size  = hc_res*hc_res;
  side_size = hc_res*hc_res2;
  if (topology->geom==VF_2Dplanar) {
    top_size  = hc_res;
    side_size = hc_res2;
  } else {
    top_size  = hc_res*hc_res;
    side_size = hc_res*hc_res2;
  }
  hc_size                 = top_size;
  delta_ptr               = enclosure->hemicube.top_deltaVF;
  enclosure->hemicube.dir = enclosure->hemicube.top_dir;
  /*
  for (patch=mesh->patches, n=0; n<mesh->npatches_g; n++, patch++) {
    patch->mask &= MASK_PLANAR|MASK_SPLIT;
  }
  */
  /*===================================================*/
  /* TAG ALL THE BSP NODES THAT ARE INFRONT/BEHIND THE */
  /* PATCH AND COMPUTE THE TRANSFORMED BOUNDING BOX OF */
  /* THE BSP NODES THAT ARE IN FRONT OF THE PATCH      */
  /*===================================================*/
  VF_SetViewPort(view, 0);
  view_plane.normal = view->view_normal;
  view_plane.d      = -V3_Dot(&(view_plane.normal), &view->view_point);
  VF_BSP_InitStack(stack);
  currentNode = topology->BSP_tree.root;
  while (currentNode != NULL) {
    currentNode->behind = VF_VoxelBehindPlane(&(currentNode->bounds), 
                                              &view_plane);
    if (!currentNode->behind) {
      if (currentNode->child[0]!=NULL) {
        VF_BSP_PushStack(stack, currentNode->child[0]);
      }
      if (currentNode->child[1]!=NULL) {
        VF_BSP_PushStack(stack, currentNode->child[1]);
      }
    }
    VF_BSP_PopStack(stack, currentNode);
  }
    
  for (VFptr=vf, i=0; i<enclosure->npatches_g; i++) *VFptr++ = 0.0;
  for (face=0; face<nfaces; face++) {
    mholes = 0;
    if (enclosure->debug_level>=VF_OUTPUT_DEBUG_1) {
      printf("       Face %d\n",face);
      fflush(stdout);
    }
    /*==================================*/
    /* SET THE VIEWPORT FOR THE CURRENT */ 
    /* FACE AND DEFINE THE VIEW PLANE   */
    /*==================================*/
    VF_SetViewPort(view, face);
    view_plane.normal = view->view_normal;
    view_plane.d      = -V3_Dot(&(view_plane.normal), &view->view_point);
    /*===================================*/
    /* REINITIALIZE THE HEMICUBE BUFFERS */
    /*===================================*/
    iptr = enclosure->hemicube.ibuffer;
    fptr = enclosure->hemicube.zbuffer;
    for (i=0; i<hc_size; i++) {
      *iptr++ = -1;
      *fptr++ = 0.0;
    }
    /*=======================*/
    /* TRAVERSE THE BSP TREE */
    /*=======================*/
    VF_BSP_InitStack(stack);
    currentNode = topology->BSP_tree.root;
    currentNode->clipped = VF_FRUSTUM_CLIP_UNKNOWN;
    while (currentNode != NULL) {
      inview = VF_FRUSTUM_CLIP_OUT;
      if (!currentNode->behind) {
        if (DO_CLIP_TO_FRUSTUM) {
          if (currentNode->clipped==VF_FRUSTUM_CLIP_IN) {
            inview = VF_FRUSTUM_CLIP_IN;
          } else {
            inview = VF_ClipBoxToFrustum(view,&(currentNode->bounds),face);
          }
        } else {
          if (face==0) {
            inview = VF_FRUSTUM_CLIP_PARTIAL;
          } else {
            if (!VF_VoxelBehindPlane(&(currentNode->bounds), &view_plane)) {
              inview = VF_FRUSTUM_CLIP_PARTIAL;
            }
          }
        }
      }
      if (inview!=VF_FRUSTUM_CLIP_OUT) {
        if (currentNode->child[0]==NULL && currentNode->child[1]==NULL) {
          /*======================================*/
          /* NO MORE BSP NODES BELOW THIS ONE SO  */
          /* PROJECT ALL THE PATCHES IN THIS NODE */
          /*======================================*/
          facet = VF_FirstFacetOfLinkList(&currentNode->members);
          while (facet) {
            if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2) {
              printf("         Projecting facet %d (gid=%d)\n",
                     facet->index,facet->patch_global_index);
            }
            facet_j = facet->index;
            patch_j = facet->patch_global_index;
            VF_FacetToPolyStack(facet, &PolyStack, &nsub_polys_j);
            for (subfacet_j=0; subfacet_j<nsub_polys_j; subfacet_j++) {
              VF_PopPolyStackPtr(&PolyStack, poly_j);
              if (enclosure->debug_level>=VF_OUTPUT_DEBUG_2 && 
                nsub_polys_j>1) {
                printf("           Sub poly_j %d of %d\n",
                       subfacet_j,nsub_polys_j);
              }
              if (facet_i==facet_j && subfacet_i==subfacet_j) continue;
              VF_ProjectOntoHemicube(poly_i, poly_j, patch_j, 
                		     view, inview, face, 
                		     &(enclosure->hemicube));
            }
            facet = VF_NextFacetOfLinkList(&currentNode->members);
          }
        } else {
          if (currentNode->clipped==VF_FRUSTUM_CLIP_IN) {
            if (currentNode->child[0]!=NULL) {
              currentNode->child[0]->clipped = VF_FRUSTUM_CLIP_IN;
            }
            if (currentNode->child[1]!=NULL) {
              currentNode->child[1]->clipped = VF_FRUSTUM_CLIP_IN;
            }
          } else {
            if (currentNode->child[0]!=NULL) {
              currentNode->child[0]->clipped = VF_FRUSTUM_CLIP_UNKNOWN;
            }
            if (currentNode->child[1]!=NULL) {
              currentNode->child[1]->clipped = VF_FRUSTUM_CLIP_UNKNOWN;
            }
          }
          switch (currentNode->axis) {
          case VF_BSP_XAXIS:
            dist = SGN(view->view_normal.x);
            break;
          case VF_BSP_YAXIS:
            dist = SGN(view->view_normal.y);
            break;
          case VF_BSP_ZAXIS:
            dist = SGN(view->view_normal.z);
            break;
          }
          if (dist>0) {
            if (currentNode->child[1]!=NULL) {
              VF_BSP_PushStack(stack, currentNode->child[1]); /* far  */
            }
            if (currentNode->child[0]!=NULL) {
              VF_BSP_PushStack(stack, currentNode->child[0]); /* near */
            }
          } else {
            if (currentNode->child[0]!=NULL) {
              VF_BSP_PushStack(stack, currentNode->child[0]); /* far  */
            }
            if (currentNode->child[1]!=NULL) {
              VF_BSP_PushStack(stack, currentNode->child[1]); /* near */
            }
          }
        }
      }
      VF_BSP_PopStack(stack, currentNode);
    }
    /*============================================*/
    /* SUM DELTA-FORMFACTORS FROM THIS PROJECTION */
    /*============================================*/
    iptr = enclosure->hemicube.ibuffer;
    if (scale_it) {
      for (i=0; i<hc_size; i++) {
        j = *iptr++;
        if (j>=0) {
          vf[j] += *delta_ptr;
        } else {
          enclosure->hemicube.nholes |= 2;
          mholes++;
        }
        delta_ptr++;
      }
    } else {
      for (i=0; i<hc_size; i++) {
        j = *iptr++;
        if (j>=0) {
          vf[j] += *delta_ptr;
        } else {
          enclosure->hemicube.nholes |= 1;
          mholes++;
        }
        delta_ptr++;
      }
    }
#ifdef PPM_DEBUG
    if (mholes!=0) {
      printf("patch %d, face %d, has %d holes!!\n",patch_ii,face,mholes);
      if (face==0) {
        WritePPM(hc_res,hc_res,enclosure->hemicube.ibuffer,enclosure->npatches_g);
      } else {
        WritePPM(hc_res,hc_res/2,enclosure->hemicube.ibuffer,enclsoure->npatches_g);
      }
    }
#endif
    view->window.y1         = hc_res2-1;
    hc_size                 = side_size;
    delta_ptr               = enclosure->hemicube.side_deltaVF;
    enclosure->hemicube.dir = enclosure->hemicube.side_dir;
  }
  if (scale_it) {
    for (i=0; i<enclosure->npatches_g; i++) {
      if (vf[i]>0.0) VF[i] += vf[i]*area_i;
    }
  } else {
    for (i=0; i<enclosure->npatches_g; i++) {
      if (vf[i]>0.0) VF[i] += vf[i];
    }
  }
}

#ifdef PPM_DEBUG
void WritePPM(int dx, int dy, int raster[], int maxval)
{
  FILE   *fp;
  char   filename[132];
  int    x, y, n, nn;
  int    r, g, b;
  int    minv, maxv;
  float  val, rr, gg, bb;
  static int cnt=0;
    
  maxv = -1;
  minv = maxval+1;
  for (y=0; y<dy; y++) {
    for (x=0; x<dx; x++, n++) {
      nn = (dy-y)*dx+x;
      if (raster[nn]>=0) {
        minv = MIN(minv,raster[nn]);
        maxv = MAX(maxv,raster[nn]);
      }
    }
  }
  sprintf(filename,"chaparral.%d.ppm",cnt); cnt++;
  fp = fopen(filename,"w");
  fprintf(fp,"P3\n");
  fprintf(fp,"%d %d\n255\n",dx,dy);
  for (n=1, y=0; y<dy; y++) {
    for (x=0; x<dx; x++, n++) {
      nn = (dy-1-y)*dx+x;
      if (raster[nn]<0) {
        r = g = b = 0;
      } else {
        val = raster[nn]/maxval;
        val = (raster[nn]-minv)/(maxv-minv);
        if (val<0.5) {
          rr = 0.0;
          gg = 2.0*val;
          bb = 1.0-gg;
        } else {
          rr = 1.0-gg;
          gg = 2.0*val-1.0;
          bb = 0.0;
        }
        r = 255*rr;
        b = 255*gg;
        g = 255*bb;
      }
      fprintf(fp," %3d %3d %3d",r,g,b);
      if (n%6==0) fprintf(fp,"\n");
    }
  }
  if (n%6!=1) fprintf(fp,"\n");
  fclose(fp);
}
#endif

