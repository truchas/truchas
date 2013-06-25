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
@(#)    $RCSfile: VF_PolyScan.c,v $
@(#)    $Revision: 1.5 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_PolyScan.c,v $
@(#)
@(#)    DESCRIPTION:  Point-sampled scan conversion of convex polygons based
@(#)    on the routines from the section on "Generic Convex Polygon Scan 
@(#)    Conversion and Clipping" by Paul Heckbert from "Graphics Gems",
@(#)    Academic Press, 1990.  It scan converts a polygon, checking
@(#)    and setting the z-buffer at each pixel. Polygon must be counter-
@(#)    clockwise.  Polygon is clipped in 2-D to the screen space window.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdio.h>
#include <math.h>

#include "vf.h"

 
void VF_PolyScan(Poly *p0, Poly *p, Window *win, int surface_num, int face, Hemicube *hc)
{
  int    x, lx, rx;
  int    y, ly, ry;
  int    i, bot, rem, nonblocking;
  int    li, li_old, ri, ri_old;
  int    index, hr, *ibuffer;
  double dx, dy, frac, ymin;
  double *zbuffer, zval, d, ax, ay, az, Ax, Ay, Az;
  double scalex, scaley, offsety, ds, ds2;
  Point  *p1, *p2, l, r, dl, dr;
  Vector *normal;
  VFtopology* topology=VF_CurrentTopology();
    
  if (p->np>VF_POLY_NMAX) {
    fprintf(stdout, "poly_scan: too many vertices: %d\n", p->np);
    return;
  }
  nonblocking = topology->nonblocking;
  hr          = hc->resolution;
  ibuffer     = hc->ibuffer;
  zbuffer     = hc->zbuffer;
  normal      = &(p0->normal);
  d           = p0->d;
  ds          = 2.0/hr;
  ds2         = 1.0/hr;
  ax          = -normal->x/d;
  ay          = -normal->y/d;
  az          = -normal->z/d;
  if (face==0) {
    scalex  = scaley = (2.0-ds)/(hr-1);
    offsety = ds2-1.0;
  } else {
    scalex  = (2.0-ds)/(hr-1);
    scaley  = (1.0-ds)/(hr/2-1);
    offsety = ds2;
  }
  Ax = ax * scalex;
  Ay = ay * scaley;
  Az = ax*(ds2-1.0) + ay*offsety - az;
  if (topology->geom==VF_2Dplanar) {
    if (face==0) {
      p1  = &p->p[0];
      p2  = &p->p[1];
      l.x = p1->x;
      r.x = p2->x;
      lx  = ceil(l.x-0.5);
      rx  = floor(r.x-0.5);
      if (nonblocking) {
        if (lx<=rx) {
          for (index=lx, x=lx; x<=rx; x++, index++) {
            zbuffer[index] = 1.0;
            ibuffer[index] = surface_num;
          }
        }
      } else {
        if (lx<=rx) {
          zval = Ax*lx + Az;
          for (index=lx, x=lx; x<=rx; x++, index++, zval+=Ax) {
            if (zval>zbuffer[index]) {
              zbuffer[index] = zval;
              ibuffer[index] = surface_num;
            }
          }
        }
      }
    } else {
      if (face==1) {
        p1 = &p->p[0];
        p2 = &p->p[1];
      } else {
        p1 = &p->p[1];
        p2 = &p->p[0];
      }
      l.y = p1->y;
      r.y = p2->y;
      ly = ceil(l.y-0.5);
      ry = floor(r.y-0.5);
      if (nonblocking) {
        if (ly<=ry) {
          for (index=ly, y=ly; y<=ry; y++, index++) {
            zbuffer[index] = 1.0;
            ibuffer[index] = surface_num;
          }
        }
      } else {
        if (ly<=ry) {
          zval = Ay*ly + Az;
          for (index=ly, y=ly; y<=ry; y++, index++, zval+=Ay) {
            if (zval>zbuffer[index]) {
              zbuffer[index] = zval;
              ibuffer[index] = surface_num;
            }
          }
        }
      }
    }
  } else {
    /*====================*/
    /* find bottom vertex */
    /*====================*/
    ymin = MAXDOUBLE;
    for (i=0; i<p->np; i++) {
      if (p->p[i].y < ymin) {
        ymin = p->p[i].y;
        bot  = i;
      }
    }
    y   = (int)ceil(ymin-0.5);  /* current scan line            */
    li  = bot;                  /* left vertex indices          */
    ri  = bot;                  /* right vertex indices         */
    ly  = y-1;                  /* lower end of left edges      */
    ry  = y-1;                  /* lower end of right edges     */
    rem = p->np;                /* number of vertices remaining */
    /*=================================================*/
    /* scan in y, activating new edges on left & right */
    /* as scan line passes over new vertices           */
    /*=================================================*/
    while (rem>0) {
      li_old = -1;
      while (ly<=y && rem>0) {
        /* advance left edge by     */
        /* stepping cw up left edge */
        rem--;
        i = li-1;
        if (i<0) i = p->np-1;
        ly     = floor(p->p[i].y+0.5);
        li_old = li;
        li     = i;
      }
      if (li_old>=0) {
        /*=============================================*/
        /* incrementalize_y: put intersection of line  */
        /* Y=y+.5 with edge between points li_old and  */
        /* li in l, put change with respect to y in dl */
        /*=============================================*/
        p1 = &p->p[li_old];
        p2 = &p->p[li];
        dx = p2->x-p1->x;
        dy = p2->y-p1->y;
        if (dy==0.0) dy = 1.0;
        frac = y + 0.5 - p1->y;
        dl.x = dx/dy;
        l.x  = p1->x+dl.x*frac;
      }
      ri_old = -1;
      while (ry<=y && rem>0) {
        /* advance right edge by      */
        /* stepping ccw up right edge */
        rem--;
        i = ri+1;
        if (i>=p->np) i = 0;
        ry     = floor(p->p[i].y+0.5);
        ri_old = ri;
        ri     = i;
      }
      if (ri_old>=0) {
        /*=============================================*/
        /* incrementalize_y: put intersection of line  */
        /* Y=y+.5 with edge between points ri_old and  */
        /* ri in r, put change with respect to y in dr */
        /*=============================================*/
        p1 = &p->p[ri_old];
        p2 = &p->p[ri];
        dx = p2->x-p1->x;
        dy = p2->y-p1->y;
        if (dy==0.0) dy = 1.0;
        frac = y + 0.5 - p1->y;
        dr.x = dx/dy;
        r.x  = p1->x+dr.x*frac;
      }
      /*===================================*/
      /* do ScanXs till end of l or r edge */
      /*===================================*/
      if (nonblocking) {
        while (y<ly && y<ry) {
          /*================================================*/
          /* process scanline by sampling polygon at Y=y+.5 */
          /*================================================*/
          lx = ceil(l.x-0.5);
          rx = floor(r.x-0.5);
          if (lx<=rx) {
            index = y*hr+lx;
            for (x=lx; x<=rx; x++, index++) {
              zbuffer[index] = 1.0;
              ibuffer[index] = surface_num;
            }
          }
          l.x += dl.x;
          r.x += dr.x;
          y++;
        }
      } else {
        while (y<ly && y<ry) {
          /*================================================*/
          /* process scanline by sampling polygon at Y=y+.5 */
          /*================================================*/
          lx = ceil(l.x-0.5);
          rx = floor(r.x-0.5);
          if (lx<=rx) {
            index = y*hr+lx;
            zval  = Ax*lx + Ay*y + Az;
            for (x=lx; x<=rx; x++, index++, zval+=Ax) {
              if (zval>zbuffer[index]) {
                zbuffer[index] = zval;
                ibuffer[index] = surface_num;
              }
            }
          }
          l.x += dl.x;
          r.x += dr.x;
          y++;
        }
      }
    }
  }
}

