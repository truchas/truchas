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
@(#)    $RCSfile: VF_CalcVF_Hottel.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcVF_Hottel.c,v $
@(#)
@(#)    DESCRIPTION:  Use Hottel's cross string method to calculate the 
@(#)    viewfactor bewtween two 2D planar surfaces.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

double 
VF_CalcVF_Hottel(Poly *poly_i, Poly *poly_j)
{
  double L1, L3, L4, L5, L6, vf=0.0;
  Point  r;

  V2_Sub(&(poly_i->p[1]), &(poly_i->p[0]), &r);
  L1 = V2_Length(&r);
  V2_Sub(&(poly_j->p[0]), &(poly_i->p[0]), &r);
  L5 = V2_Length(&r);
  V2_Sub(&(poly_j->p[1]), &(poly_i->p[0]), &r);
  L4 = V2_Length(&r);
  V2_Sub(&(poly_j->p[0]), &(poly_i->p[1]), &r);
  L3 = V2_Length(&r);
  V2_Sub(&(poly_j->p[1]), &(poly_i->p[1]), &r);
  L6 = V2_Length(&r);
  vf = ((L5+L6)-(L3+L4))/(2.0*L1);
  return (vf);
}
