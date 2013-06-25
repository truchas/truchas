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
@(#)    $RCSfile: VF_BSPtree.c,v $
@(#)    $Revision: 1.2.4.1 $  $Date: 2006/12/05 19:44:48 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_BSPtree.c,v $
@(#)
@(#)    DESCRIPTION:  Build a BSP tree to partition the facets.  This is
@(#)    based on the BSP code in GraphicsGems III.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "vf.h"

static VFtopology *Topology;
/*******************************************************************************
 INITIALIZE AND START THE BUILDING OF BSP TREE.
 ******************************************************************************/
void VF_BSP_CreateTree(VFtopology* topology)
{
  int     n;
  BinTree *BSPTree;
  static int* map=NULL;
  static int  map_size=0;
  
  /*VF_GetINTbuffer_ptr(&map);*/
  Topology                = topology;
  BSPTree                 = &(topology->BSP_tree);
  BSPTree->bounds         = topology->bounds;
  BSPTree->MaxDepth       = topology->bsp_depth;
  BSPTree->MaxListLength  = topology->bsp_length;
  BSPTree->members.head   = NULL;
  BSPTree->members.tail   = NULL;
  BSPTree->members.length = topology->nfacets_g;
  BSPTree->facetLinkList  = (struct FacetLink*)VF_Newv(topology->nfacets_g*
                                                       sizeof(struct FacetLink));
  if (map_size==0) {
    map = VF_Newi(topology->nfacets_g);
  } else if (map_size<topology->nfacets_g) {
    map = VF_ReNewi(map,topology->nfacets_g);
  }
  for (n=0; n<topology->nfacets_g;  n++) {
    map[n] = topology->facets[n].index;
  }
  for (n=0; n<topology->nfacets_g;  n++) {
    BSPTree->facetLinkList[n].type  = 0;
    BSPTree->facetLinkList[n].next  = NULL;
    BSPTree->facetLinkList[n].facet = &(topology->facets[map[n]]);
    if (n==0) {
      BSPTree->members.head = &BSPTree->facetLinkList[n];
      BSPTree->members.tail = &BSPTree->facetLinkList[n];
    } else {
      BSPTree->members.tail->next = &(BSPTree->facetLinkList[n]);
      BSPTree->members.tail       = &(BSPTree->facetLinkList[n]);
    } 
  }
  /*===================================================================*/
  /* START BUILDING THE BSP TREE BY SUBDIVIDING ALONG THE LONGEST AXIS */
  /*===================================================================*/
  BSPTree->root                          = (BinNodePtr)VF_Newv(sizeof(BinNode));
  BSPTree->root->bounds                  = BSPTree->bounds;
  BSPTree->root->members                 = BSPTree->members;
  BSPTree->root->child[0]                = NULL;
  BSPTree->root->child[1]                = NULL;
  VF_BSP_PartitionNode(BSPTree->root, 0, BSPTree->MaxDepth, BSPTree->MaxListLength);
  topology->BSP_stack = (BSPstackPtr)VF_Newv(sizeof(BSPstack));
}

/*******************************************************************************
 Builds the BSP tree by subdividing along the center of x, y, or z bounds, once
 each time this function is called. This function calls itself recursively until
 either the tree is deeper than MaxDepth or all of the tree leaves contains less
 than MaxListLength of objects.

 Entry:
   node          - node currently working on
   depth         - current tree depth
   MaxDepth      - Max allowed tree depth
   MaxList       - Max allowed object list length of a leave node
 ******************************************************************************/
void VF_BSP_PartitionNode(BinNodePtr node, int depth, int MaxDepth, int MaxList)
{
  int    v, mask, doit=0;
  double xmin_f, ymin_f, zmin_f;
  double xmax_f, ymax_f, zmax_f;
  double xmin_b, ymin_b, zmin_b;
  double xmax_b, ymax_b, zmax_b;
  Poly   poly;
  Facet  *facet;
  struct FacetLink *currentLink, *nextLink;

  node->depth    = depth;
  node->child[0] = NULL;
  node->child[1] = NULL;
  if ((node->members.length > MaxList) && (depth < MaxDepth)) {
    doit = 1;
    mask = VF_BSP_XAXIS | VF_BSP_YAXIS | VF_BSP_ZAXIS;
    if (!VF_BSP_FindPartition(node,mask)) {
      mask ^= node->axis; 
      if (!VF_BSP_FindPartition(node,mask)) {
        mask ^= node->axis; 
        if (!VF_BSP_FindPartition(node,mask)) {
          doit = 0;
        }
      }
    }
  }
  if (doit) {
    xmin_f = ymin_f = zmin_f =  MAXDOUBLE;
    xmax_f = ymax_f = zmax_f = -MAXDOUBLE;
    xmin_b = ymin_b = zmin_b =  MAXDOUBLE;
    xmax_b = ymax_b = zmax_b = -MAXDOUBLE;
    currentLink = node->members.head;
    while (currentLink) {
      nextLink = currentLink->next;
      facet    = currentLink->facet;
      switch (VF_BSP_IsFacetInNode(facet, node->d, node->axis)) {
      case VF_BSP_FRONT:
      case VF_BSP_SPLIT:
        /* REMOVE LINK FROM CURRENT NODE */
        currentLink->next = NULL;
        /* AND ADD IT TO THE FRONT CHILD */
        if (node->child[0] == NULL) {
          node->child[0]                 = (BinNodePtr)VF_Newv(sizeof(BinNode));
          node->child[0]->behind         = 0;
          node->child[0]->bounds         = node->bounds;
          node->child[0]->members.head   = currentLink;
          node->child[0]->members.tail   = NULL;
          node->child[0]->members.length = 0;
        } else {
          node->child[0]->members.tail->next = currentLink;
        }
        node->child[0]->members.tail = currentLink;
        node->child[0]->members.length++;
        facet = currentLink->facet;
        VF_FacetToPoly(facet, &poly);
        for (v=0; v<poly.np; v++) {
          xmin_f = MIN(xmin_f,poly.p[v].x);
          ymin_f = MIN(ymin_f,poly.p[v].y);
          zmin_f = MIN(zmin_f,poly.p[v].z);
          xmax_f = MAX(xmax_f,poly.p[v].x);
          ymax_f = MAX(ymax_f,poly.p[v].y);
          zmax_f = MAX(zmax_f,poly.p[v].z);
        }
        break;
      case VF_BSP_BACK:
        /* REMOVE LINK FROM CURRENT NODE */
        currentLink->next = NULL;
        /* AND ADD IT TO THE BACK CHILD  */
        if (node->child[1] == NULL) {
          node->child[1]                 = (BinNodePtr)VF_Newv(sizeof(BinNode));
          node->child[1]->behind         = 0;
          node->child[1]->bounds         = node->bounds;
          node->child[1]->members.head   = currentLink;
          node->child[1]->members.tail   = NULL;
          node->child[1]->members.length = 0;
        } else {
          node->child[1]->members.tail->next = currentLink;
        }
        node->child[1]->members.tail = currentLink;
        node->child[1]->members.length++;
        facet = currentLink->facet;
        VF_FacetToPoly(facet, &poly);
        for (v=0; v<poly.np; v++) {
          xmin_b = MIN(xmin_b,poly.p[v].x);
          ymin_b = MIN(ymin_b,poly.p[v].y);
          zmin_b = MIN(zmin_b,poly.p[v].z);
          xmax_b = MAX(xmax_b,poly.p[v].x);
          ymax_b = MAX(ymax_b,poly.p[v].y);
          zmax_b = MAX(zmax_b,poly.p[v].z);
        }
        break;
      }
      currentLink = nextLink;
    }
    if (node->child[0]==NULL || node->child[1]==NULL) {
      printf("BSPtree fatal error!!\n");
      VF_Exit(-1);
    }
    if (node->child[0]!=NULL) {
      node->child[0]->bounds.xmin = xmin_f;
      node->child[0]->bounds.ymin = ymin_f;
      node->child[0]->bounds.zmin = zmin_f;
      node->child[0]->bounds.xmax = xmax_f;
      node->child[0]->bounds.ymax = ymax_f;
      node->child[0]->bounds.zmax = zmax_f;
    }
    if (node->child[1]!=NULL) {
      node->child[1]->bounds.xmin = xmin_b;
      node->child[1]->bounds.ymin = ymin_b;
      node->child[1]->bounds.zmin = zmin_b;
      node->child[1]->bounds.xmax = xmax_b;
      node->child[1]->bounds.ymax = ymax_b;
      node->child[1]->bounds.zmax = zmax_b;
    }
    node->members.head   = NULL;
    node->members.tail   = NULL;
    node->members.length = 0;
    if (node->child[0]!=NULL) 
      VF_BSP_PartitionNode(node->child[0], depth+1, MaxDepth, MaxList);
    if (node->child[1]!=NULL) 
      VF_BSP_PartitionNode(node->child[1], depth+1, MaxDepth, MaxList);
  }
}

int VF_BSP_FindPartition(BinNodePtr node, int mask)
{
  int    i, majorAxis=0, cnt0, cnt1, min_diff, status=1;
  double dx, dy, dz, dmax=0.0, sum=0.0;
  double davg, d1, d2, dnew;
  Facet  *facet;
  Poly   poly;
  Point  center;
  struct FacetLink *currentLink;
    
  /*=================*/
  /* FIND MAJOR AXIS */
  /*=================*/
  if (VF_BSP_XAXIS&mask) {
    dx = node->bounds.xmax - node->bounds.xmin;
    if (dx>dmax) {
      dmax      = dx;
      majorAxis = VF_BSP_XAXIS;
    }         
  }
  if (VF_BSP_YAXIS&mask) {
    dy = node->bounds.ymax - node->bounds.ymin;
    if (dy>dmax) {
      dmax      = dy;
      majorAxis = VF_BSP_YAXIS;
    }
  } 
  if (VF_BSP_ZAXIS&mask) {
    dz = node->bounds.zmax - node->bounds.zmin;
    if (dz>dmax) {
      dmax      = dz;
      majorAxis = VF_BSP_ZAXIS;
    }
  }
  node->axis = majorAxis;
  /* FOR FIRST GUESS, USE MEAN CENTER POSITION */
  currentLink = node->members.head;
  while (currentLink)  {
    facet = currentLink->facet;
    VF_FacetToPoly(facet, &poly);
    VF_PolyCenter(&poly, &center);
    switch (majorAxis) {
    case VF_BSP_XAXIS:
      sum += center.x;
      break;
    case VF_BSP_YAXIS:
      sum += center.y;
      break;
    case VF_BSP_ZAXIS:
      sum += center.z;
      break;
    }
    currentLink = currentLink->next;
  }
  node->d = sum/(node->members.length);
    
  davg = node->d;
  cnt0 = cnt1 = 0;
  currentLink = node->members.head;
  while (currentLink) {
    facet = currentLink->facet;
    switch (VF_BSP_IsFacetInNode(facet, davg, majorAxis)) {
    case VF_BSP_FRONT:
    case VF_BSP_SPLIT:
      cnt0++;
      break;
    case VF_BSP_BACK:
      cnt1++;
      break;
    }
    currentLink = currentLink->next;
  }
  min_diff = abs(cnt0-cnt1);
  if (min_diff<5) return status;
    
  if (cnt0<cnt1) {
    d1 = davg;
    switch (majorAxis) {
    case VF_BSP_XAXIS:
      d2 = node->bounds.xmax;
      break;
    case VF_BSP_YAXIS:
      d2 = node->bounds.ymax;
      break;
    case VF_BSP_ZAXIS:
      d2 = node->bounds.zmax;
      break;
    }
  } else {
    d2 = davg;
    switch (majorAxis) {
    case VF_BSP_XAXIS:
      d1 = node->bounds.xmin;
      break;
    case VF_BSP_YAXIS:
      d1 = node->bounds.ymin;
      break;
    case VF_BSP_ZAXIS:
      d1 = node->bounds.zmin;
      break;
    }
  }
  dnew = 0.5*(d1+d2);
  for (i=0; i<15; i++) {
    cnt0 = cnt1 = 0;
    currentLink = node->members.head;
    while (currentLink) {
      facet = currentLink->facet;
      switch (VF_BSP_IsFacetInNode(facet, dnew, majorAxis)) {
      case VF_BSP_FRONT:
      case VF_BSP_SPLIT:
        cnt0++;
        break;
      case VF_BSP_BACK:
        cnt1++;
        break;
      }
      currentLink = currentLink->next;
    }
    if (abs(cnt0-cnt1)<min_diff) {
      min_diff = abs(cnt0-cnt1);
      node->d  = dnew;
    }
    if (min_diff<5) {
      break;
    }
    if (cnt0<cnt1) {
      d1 = dnew;
    } else {
      d2 = dnew;
    }
    dnew = 0.5*(d1+d2);
  }
  if (cnt0==0 || cnt1==0) status = 0;
  return status;
}

int VF_BSP_IsFacetInNode(Facet *facet, double d, int axis)
{
  int    v, on=0, front=0, back=0;
  Poly   poly;
        
  VF_FacetToPoly(facet, &poly);
  for (v=0; v<poly.np; v++){
    switch (axis) {
    case VF_BSP_XAXIS:
      if (poly.p[v].x == d) on++;
      if (poly.p[v].x <  d) front++;
      if (poly.p[v].x >  d) back++;
      break;
    case VF_BSP_YAXIS:
      if (poly.p[v].y == d) on++;
      if (poly.p[v].y <  d) front++;
      if (poly.p[v].y >  d) back++;
      break;
    case VF_BSP_ZAXIS:
      if (Topology->geom!=VF_2Dplanar) {
        if (poly.p[v].z == d) on++;
        if (poly.p[v].z <  d) front++;
        if (poly.p[v].z >  d) back++;
      }
      break;
    }
  }
  if (on==poly.np) {
    return VF_BSP_SPLIT;
  } else if (front+on==poly.np) {
    return VF_BSP_FRONT;
  } else if (back+on==poly.np) {
    return VF_BSP_BACK;
  } else {
    return VF_BSP_SPLIT;
  }
}

void VF_BSP_DeleteTree(BinNodePtr node)
{
  struct FacetLink *currentLink, *nextLink;
    
  if (node->child[0]!=NULL) {
    VF_BSP_DeleteTree(node->child[0]);
    node->child[0] = NULL;
  }
  if (node->child[1]!=NULL) {
    VF_BSP_DeleteTree(node->child[1]);
    node->child[1] = NULL;
  }
  currentLink = node->members.head;
  if (currentLink!=NULL) nextLink = currentLink->next;
  while (currentLink) {
    if (currentLink->type) VF_Free(currentLink);
    currentLink = nextLink;
    if (currentLink!=NULL) nextLink = currentLink->next;
  }
  VF_Free(node);
}

Facet* VF_FirstFacetOfLinkList(FacetList *facetList)
{
  if (facetList->head==NULL) return(NULL);
  facetList->current = facetList->head;
  return (facetList->current->facet);
}

Facet* VF_NextFacetOfLinkList(FacetList *facetList)
{
  if (facetList->current->next == NULL) return(NULL);
  facetList->current = facetList->current->next;
  return (facetList->current->facet);
}
