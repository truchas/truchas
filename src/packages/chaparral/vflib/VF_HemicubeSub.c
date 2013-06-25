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
@(#)    $RCSfile: VF_HemicubeSub.c,v $
@(#)    $Revision: 1.3 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_HemicubeSub.c,v $
@(#)
@(#)    DESCRIPTION:  Use hemicube method to calculate the viewfactor matrix.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

void VF_HemicubeSub(int facet_i, int subfacet_i, Poly *poly_i, 
                    ViewPort *view, double VF[], double vf[], int iseg, int jseg)
{   
  int	      k, is, js, nnn;
  double      *VFptr,rowsum, area;
  double      phi, phi1, phi2;
  Point       Phi1, Phi2, Phi3, Phi4;
  Poly        poly;
  VFenclosure *enclosure=VF_CurrentEnclosure();
  
  if (enclosure->debug_level>=VF_OUTPUT_DEBUG_1 && VFLIB_Rank==0) {
    printf("	   MULTI-LEVEL Pass for patch %d -------------------\n",
           facet_i);
    printf("	     Recalculating viewfactors with iseg=%d and jseg=%d\n",
           iseg,jseg);
  }
  poly.np     = poly_i->np;
  poly.d      = poly_i->d;
  poly.normal = poly_i->normal;
  for (is=0; is<iseg; is++) {
    if (poly_i->np == 4 || poly_i->np == 2) {
      nnn = jseg;
    } else {
      nnn	  = 2*(is+1)-1;
      phi	  = (double)(is)/(double)(iseg);
      Phi1.x	  = (poly_i->p[1].x-poly_i->p[0].x)*
    		    phi+poly_i->p[0].x;
      Phi1.y	  = (poly_i->p[1].y-poly_i->p[0].y)*
    		    phi+poly_i->p[0].y;
      Phi1.z	  = (poly_i->p[1].z-poly_i->p[0].z)*
    		    phi+poly_i->p[0].z;
      Phi2.x	  = (poly_i->p[2].x-poly_i->p[0].x)*
    		    phi+poly_i->p[0].x;
      Phi2.y	  = (poly_i->p[2].y-poly_i->p[0].y)*
    		    phi+poly_i->p[0].y;
      Phi2.z	  = (poly_i->p[2].z-poly_i->p[0].z)*
    		    phi+poly_i->p[0].z;
      phi	  = (double)(is+1)/(double)(iseg);
      Phi3.x	  = (poly_i->p[1].x-poly_i->p[0].x)*
    		    phi+poly_i->p[0].x;
      Phi3.y	  = (poly_i->p[1].y-poly_i->p[0].y)*
    		    phi+poly_i->p[0].y;
      Phi3.z	  = (poly_i->p[1].z-poly_i->p[0].z)*
    		    phi+poly_i->p[0].z;
      Phi4.x	  = (poly_i->p[2].x-poly_i->p[0].x)*
    		    phi+poly_i->p[0].x;
      Phi4.y	  = (poly_i->p[2].y-poly_i->p[0].y)*
    		    phi+poly_i->p[0].y;
      Phi4.z	  = (poly_i->p[2].z-poly_i->p[0].z)*
    		    phi+poly_i->p[0].z;
      poly.p[2].x = Phi1.x;
      poly.p[2].y = Phi1.y;
      poly.p[2].z = Phi1.z;
      poly.p[1].x = Phi3.x;
      poly.p[1].y = Phi3.y;
      poly.p[1].z = Phi3.z;
    }
    for (js=0; js<nnn; js++) {
      if (poly_i->np == 4 || poly_i->np == 2) {
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
    	  phi	      = (double)(1+js/2)/(double)(is+1);
    	  poly.p[0]   = poly.p[2];
    	  poly.p[2].x = (Phi4.x-Phi3.x)*phi+Phi3.x;
    	  poly.p[2].y = (Phi4.y-Phi3.y)*phi+Phi3.y;
    	  poly.p[2].z = (Phi4.z-Phi3.z)*phi+Phi3.z;
    	} else {
    	  phi	      = (double)((1+js)/2)/(double)(is);
    	  poly.p[1]   = poly.p[2];
    	  poly.p[2].x = (Phi2.x-Phi1.x)*phi+Phi1.x;
    	  poly.p[2].y = (Phi2.y-Phi1.y)*phi+Phi1.y;
    	  poly.p[2].z = (Phi2.z-Phi1.z)*phi+Phi1.z;
    	}
      }
      area = VF_PolyArea(&poly);
      VF_SetView(view, &poly);
      VF_HemicubeProjectRow(facet_i, subfacet_i, 
    			    &poly, view, 1, area, VF, vf);
      if (enclosure->debug_level>=VF_OUTPUT_DEBUG_1 && VFLIB_Rank==0) {
    	for (rowsum=0.0, VFptr=vf, k=0; k<enclosure->npatches_g; k++) {
    	  rowsum += *VFptr++;
    	}
    	printf("	 iseg = %d   jseg = %d\n",is, js);
    	printf("	   ViewPnt = (%g, %g, %g)\n",
    	       view->view_point.x, 
    	       view->view_point.y, 
    	       view->view_point.z);
    	printf("	   rowsum  = %g (%g)\n",rowsum,1.0-rowsum);
      }
    } /* END OF SUB-ELEMENT LOOP */
  } /* END OF SUB-ELEMENT LOOP */
}
