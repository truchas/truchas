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
@(#)    $RCSfile: VF_BoundingBox.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_BoundingBox.c,v $
@(#)
@(#)    DESCRIPTION:  Bounding box utilities for the adaptive method.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vf.h"

void VF_ComputeFacetBoundingBox(Box *box, Facet *facet)
{
    int    n, vertex, vertex0, vertex1;
    double theta,cos0,cos1,sin0,sin1;
    VFtopology *topology = VF_CurrentTopology();
    
    switch (topology->geom) {
    case VF_2Daxisym:
        theta = topology->theta*(double)(facet->sector);
        cos0  = cos(theta);
        sin0  = sin(theta);
        if (facet->sector+1 == topology->nrotations) {
            theta = 0.0;
        }
        cos1      = cos(theta);
        sin1      = sin(theta);
        vertex0   = facet->vertex_list[0];
        vertex1   = facet->vertex_list[1];
        box->xmin = cos0*topology->x[vertex0];
        box->ymin = topology->y[vertex0];
        box->zmin = sin0*topology->x[vertex0];
        box->xmax = box->xmin;
        box->ymax = box->ymin;
        box->zmax = box->zmin;
        box->xmin = MIN(box->xmin,cos0*topology->x[vertex1]);
        box->ymin = MIN(box->ymin,topology->y[vertex1]);
        box->zmin = MIN(box->zmin,sin0*topology->x[vertex1]);
        box->xmax = MAX(box->xmax,cos0*topology->x[vertex1]);
        box->ymax = MAX(box->ymax,topology->y[vertex1]);
        box->zmax = MAX(box->zmax,sin0*topology->x[vertex1]);
        if (facet->num_vertices==4) {
            box->xmin = MIN(box->xmin,cos1*topology->x[vertex0]);
            box->ymin = MIN(box->ymin,topology->y[vertex0]);
            box->zmin = MIN(box->zmin,sin1*topology->x[vertex0]);
            box->xmax = MAX(box->xmax,cos1*topology->x[vertex0]);
            box->ymax = MAX(box->ymax,topology->y[vertex0]);
            box->zmax = MAX(box->zmax,sin1*topology->x[vertex0]);
            box->xmin = MIN(box->xmin,cos1*topology->x[vertex1]);
            box->ymin = MIN(box->ymin,topology->y[vertex1]);
            box->zmin = MIN(box->zmin,sin1*topology->x[vertex1]);
            box->xmax = MAX(box->xmax,cos1*topology->x[vertex1]);
            box->ymax = MAX(box->ymax,topology->y[vertex1]);
            box->zmax = MAX(box->zmax,sin1*topology->x[vertex1]);
        } else {
            if (fabs((double)(topology->x[vertex0])) < topology->spatial_tol) {
                vertex = vertex1;
            }
            if (fabs((double)(topology->x[vertex1])) < topology->spatial_tol) {
                vertex = vertex0;
            }
            box->xmin = MIN(box->xmin,cos1*topology->x[vertex]);
            box->ymin = MIN(box->ymin,topology->y[vertex]);
            box->zmin = MIN(box->zmin,sin1*topology->x[vertex]);
            box->xmax = MAX(box->xmax,cos1*topology->x[vertex]);
            box->ymax = MAX(box->ymax,topology->y[vertex]);
            box->zmax = MAX(box->zmax,sin1*topology->x[vertex]);
        }
        break;
    case VF_2Dplanar:
        vertex0   = facet->vertex_list[0];
        vertex1   = facet->vertex_list[1];
        box->xmin = topology->x[vertex0];
        box->ymin = topology->y[vertex0];
        box->xmax = box->xmin;
        box->ymax = box->ymin;
        box->xmin = MIN(box->xmin,topology->x[vertex1]);
        box->ymin = MIN(box->ymin,topology->y[vertex1]);
        box->xmax = MAX(box->xmax,topology->x[vertex1]);
        box->ymax = MAX(box->ymax,topology->y[vertex1]);
        box->zmin = 0.0;
        box->zmax = 0.0;
        break;
    case VF_3D:
        vertex    = facet->vertex_list[0];
        box->xmin = topology->x[vertex];
        box->ymin = topology->y[vertex];
        box->zmin = topology->z[vertex];
        box->xmax = box->xmin;
        box->ymax = box->ymin;
        box->zmax = box->zmin;
        for (n=1; n<facet->num_vertices; n++) {
            vertex    = facet->vertex_list[n];
            box->xmin = MIN(box->xmin,topology->x[vertex]);
            box->ymin = MIN(box->ymin,topology->y[vertex]);
            box->zmin = MIN(box->zmin,topology->z[vertex]);
            box->xmax = MAX(box->xmax,topology->x[vertex]);
            box->ymax = MAX(box->ymax,topology->y[vertex]);
            box->zmax = MAX(box->zmax,topology->z[vertex]);
        }
        break;
    }
}

void VF_ComputePolyBoundingBox(Box *box, Poly *poly)
{
    int  n;
    VFtopology *topology = VF_CurrentTopology();

    switch (topology->geom) {
    case VF_2Dplanar:
        box->xmin = poly->p[0].x;
        box->ymin = poly->p[0].y;
        box->xmax = box->xmin;
        box->ymax = box->ymin;
        box->xmin = MIN(box->xmin,poly->p[1].x);
        box->ymin = MIN(box->ymin,poly->p[1].y);
        box->xmax = MAX(box->xmax,poly->p[1].x);
        box->ymax = MAX(box->ymax,poly->p[1].y);
        box->zmin = 0.0;
        box->zmax = 0.0;
        break;
    case VF_2Daxisym:
    case VF_3D:
        box->xmin = poly->p[0].x;
        box->ymin = poly->p[0].y;
        box->zmin = poly->p[0].z;
        box->xmax = box->xmin;
        box->ymax = box->ymin;
        box->zmax = box->zmin;
        for (n=1; n<poly->np; n++) {
            box->xmin = MIN(box->xmin,poly->p[n].x);
            box->ymin = MIN(box->ymin,poly->p[n].y);
            box->zmin = MIN(box->zmin,poly->p[n].z);
            box->xmax = MAX(box->xmax,poly->p[n].x);
            box->ymax = MAX(box->ymax,poly->p[n].y);
            box->zmax = MAX(box->zmax,poly->p[n].z);
        }
        break;
    }
}

void VF_ExtentsBoundingBox(Box *box_i, Box *box_j, 
                           Box *extents, Box *minmax)
{
    minmax->xmin = 0.0;
    minmax->ymin = 0.0;
    minmax->zmin = 0.0;
    minmax->xmax = 0.0;
    minmax->ymax = 0.0;
    minmax->zmax = 0.0;
    /*=========================*/
    /* CHECK MIN X COORDINATES */
    /*=========================*/
    if (box_i->xmin < box_j->xmin) {
        extents->xmin = box_i->xmin;
        minmax->xmin  = -1.0;
    } else {
        extents->xmin = box_j->xmin;
        minmax->xmin  = 1.0;
    }  
    /*=========================*/
    /* CHECK MIN Y COORDINATES */
    /*=========================*/
    if (box_i->ymin < box_j->ymin) {
        extents->ymin = box_i->ymin;
        minmax->ymin  = -1.0;
    } else { 
        extents->ymin = box_j->ymin;
        minmax->ymin  = 1.0;
    }  
    /*=========================*/
    /* CHECK MIN Z COORDINATES */
    /*=========================*/
    if (box_i->zmin < box_j->zmin) {
        extents->zmin = box_i->zmin;
        minmax->zmin  = -1.0;
    } else {
        extents->zmin = box_j->zmin;
        minmax->zmin  = 1.0;
    }  

    /*=========================*/
    /* CHECK MAX X COORDINATES */
    /*=========================*/
    if (box_i->xmax > box_j->xmax) {
        extents->xmax = box_i->xmax;
        minmax->xmax  = -1.0;
    } else { 
        extents->xmax = box_j->xmax;
        minmax->xmax  = 1.0;
    }  
    /*=========================*/
    /* CHECK MAX Y COORDINATES */
    /*=========================*/
    if (box_i->ymax > box_j->ymax) {
        extents->ymax = box_i->ymax;
        minmax->ymax  = -1.0;
    } else { 
        extents->ymax = box_j->ymax;
        minmax->ymax  = 1.0;
    }  
    /*=========================*/
    /* CHECK MAX Z COORDINATES */
    /*=========================*/
    if (box_i->zmax > box_j->zmax) {
        extents->zmax = box_i->zmax;
        minmax->zmax  = -1.0;
    } else {
        extents->zmax = box_j->zmax;
        minmax->zmax  = 1.0;
    }
}
