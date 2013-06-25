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
@(#)    $RCSfile: VF_CalcVF_Gauss.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_CalcVF_Gauss.c,v $
@(#)
@(#)    DESCRIPTION:  Use gauss quadrature to calculate the viewfactor
@(#)    between 2 unoccluded polygons.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vf.h"

void gauss_coeff(double x[], double w[], int n);
void get_shape(Poly *poly_i, double u, double v, 
               double phi[], double *jacobian);

double 
VF_CalcVF_Gauss(Poly *poly_i, Poly *poly_j, int nseg_i, int nseg_j)
{
  int     i, n, ni1, ni2, nj1, nj2;
  double  r, x_i, y_i, z_i, x_j, y_j, z_j, u_i, v_i, u_j, v_j, cos_i, cos_j;
  double  vf=0.0, *gp_i, *gp_j, *wt_i, *wt_j, phi[4], jacobian_i, jacobian_j;
  static  int    *offset=NULL,max_pnts=20;
  static  double *gauss_pnt, *gauss_wts;

  if (offset==NULL) {
    offset    = VF_Newi(max_pnts);
    offset[0] = 0;
    for (n=1, i=1; i<max_pnts; i++) {
      offset[i] = offset[i-1]+i;
      n        += i+1;
    }
    gauss_pnt = VF_Newd(n);
    gauss_wts = VF_Newd(n);
    for (i=0; i<max_pnts; i++) {
      n = offset[i];
      gauss_coeff(&gauss_pnt[n], &gauss_wts[n], i+1);
    }
  }
  if (nseg_i>max_pnts || nseg_j>max_pnts) {
    printf("Maximum gauss points exceeded!! (%d || %d) > %d\n",
           nseg_i,nseg_j,max_pnts);
    VF_Exit(0);
  }
    
  gp_i = &gauss_pnt[offset[nseg_i-1]];
  gp_j = &gauss_pnt[offset[nseg_j-1]];
  wt_i = &gauss_wts[offset[nseg_i-1]];
  wt_j = &gauss_wts[offset[nseg_j-1]];
  for (ni1=0; ni1<nseg_i; ni1++) {
    u_i = gp_i[ni1];
    for (ni2=0; ni2<nseg_i; ni2++) {
      v_i = gp_i[ni2];
      get_shape(poly_i,u_i,v_i,phi,&jacobian_i);
      for (x_i=0.0, y_i=0.0, z_i=0.0, n=0; n<poly_i->np; n++) {
        x_i += phi[n]*poly_i->p[n].x;
        y_i += phi[n]*poly_i->p[n].y;
        z_i += phi[n]*poly_i->p[n].z;
      }
      for (nj1=0; nj1<nseg_j; nj1++) {
        u_j = gp_j[nj1];
        for (nj2=0; nj2<nseg_j; nj2++) {
          v_j = gp_j[nj2];
          get_shape(poly_j,u_j,v_j,phi,&jacobian_j);
          for (x_j=0.0, y_j=0.0, z_j=0.0, n=0; n<poly_j->np; n++) {
                        x_j += phi[n]*poly_j->p[n].x;
                        y_j += phi[n]*poly_j->p[n].y;
                        z_j += phi[n]*poly_j->p[n].z;
          }
          r     = sqrt((x_j-x_i)*(x_j-x_i)+
                       (y_j-y_i)*(y_j-y_i)+
                       (z_j-z_i)*(z_j-z_i));
          cos_i = ((x_j-x_i)*poly_i->normal.x+
                   (y_j-y_i)*poly_i->normal.y+
                    (z_j-z_i)*poly_i->normal.z)/r;
          cos_j = ((x_i-x_j)*poly_j->normal.x+
                   (y_i-y_j)*poly_j->normal.y+
                   (z_i-z_j)*poly_j->normal.z)/r;
          if (cos_i>0.0 && cos_j>0.0) {
            vf += wt_i[ni1]*wt_i[ni2]*wt_j[nj1]*wt_j[nj2]*
                  jacobian_i*jacobian_j*cos_i*cos_j/(r*r);
          }
        }
      }
    }
  }
  vf /= M_PI*VF_PolyArea(poly_i);
  return (vf);
}

void get_shape(Poly *poly, double u, double v, double phi[], double *jacobian)
{
  int    n;
  double dpdu[4], dpdv[4];
  double dxdu=0.0, dydu=0.0, dzdu=0.0;
  double dxdv=0.0, dydv=0.0, dzdv=0.0;
  double e1_dot_e1, e2_dot_e2, e1_dot_e2;

  /*==========================*/
  /* EVALUATE SHAPE FUNCTIONS */
  /*==========================*/
  phi[0] = 0.25*(1.0-u)*(1.0-v);
  phi[1] = 0.25*(1.0+u)*(1.0-v);
  phi[2] = 0.25*(1.0+u)*(1.0+v);
  phi[3] = 0.25*(1.0-u)*(1.0+v);
  /*===============================================*/
  /* DERIVATIVES OF SHAPE FUNCTIONS W.R.T. U AND V */
  /*===============================================*/
  dpdu[0] = -0.25*(1.0-v);
  dpdu[1] = -dpdu[0];
  dpdu[2] =  0.25*(1.0+v);
  dpdu[3] = -dpdu[2];
  dpdv[0] = -0.25*(1.0-u);
  dpdv[1] = -0.25*(1.0+u);
  dpdv[2] = -dpdv[1];
  dpdv[3] = -dpdv[0];
  /*=====================================*/
  /* TREAT TRIANGLE AS DEGENERATIVE QUAD */
  /*=====================================*/
  if (poly->np==3) {
    phi[2]  += phi[3];
    dpdu[2] += dpdu[3];
    dpdv[2] += dpdv[3];
    phi[3]   = 0.0;
    dpdu[3]  = 0.0;
    dpdv[3]  = 0.0;
  }
  /*===========================================*/
  /* DERIVATIVES OF X, Y, AND Z W.R.T. U AND V */
  /*===========================================*/
  for (n=0; n<poly->np; n++) {
    dxdu += dpdu[n]*poly->p[n].x;
    dydu += dpdu[n]*poly->p[n].y;
    dzdu += dpdu[n]*poly->p[n].z;
    dxdv += dpdv[n]*poly->p[n].x;
    dydv += dpdv[n]*poly->p[n].y;
    dzdv += dpdv[n]*poly->p[n].z;
  }
  /*=====================================*/
  /* COMPUTE DETERMINANT OF THE JACOBIAN */
  /*=====================================*/
  e1_dot_e1 = dxdu*dxdu+dydu*dydu+dzdu*dzdu;
  e2_dot_e2 = dxdv*dxdv+dydv*dydv+dzdv*dzdv;
  e1_dot_e2 = dxdu*dxdv+dydu*dydv+dzdu*dzdv;
  *jacobian = sqrt(e1_dot_e1*e2_dot_e2-e1_dot_e2*e1_dot_e2);
}

#define EPS 3.0e-11

void gauss_coeff(double x[], double w[], int n)
{
  int    i, j, m;
  double x1 = -1.0e0;
  double x2 =  1.0e0;
  double z1,z,xm,xl,pp,p3,p2,p1;

  m  = (n+1)/2;
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  for (i=1; i<=m; i++)  {
    z = cos(M_PI*(i-0.25)/(n+0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (j=1; j<=n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp = n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z  = z1-p1/pp;
    } while (fabs(z-z1) > EPS);
    x[i-1] = xm-xl*z;
    x[n-i] = xm+xl*z;
    w[i-1] = 2.0*xl/((1.0-z*z)*pp*pp);
    w[n-i] = w[i-1];
  }
}

#undef EPS
