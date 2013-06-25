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
@(#)    $RCSfile: VF_FindCandidates.c,v $
@(#)    $Revision: 1.2 $  $Date: 2005/09/08 16:41:25 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_FindCandidates.c,v $
@(#)
@(#)    DESCRIPTION:  Find any surfaces that may potentially occlude the
@(#)    view between 2 surfaces.  Use backface culling and shaft culling
@(#)    to test surfaces.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "vf.h"

#define LIST_INCR 500

void VF_InitCandidates(Facet *facet) 
{
  BinNodePtr  currentNode;
  BSPstackPtr stack;
  BinTree     *BSPtree;
  VFtopology  *topology=VF_CurrentTopology();

  stack       = topology->BSP_stack;
  BSPtree     = &(topology->BSP_tree);
  currentNode = BSPtree->root;
  VF_BSP_InitStack(stack);
  /*====================================================================*/
  /* TRAVERSE THE BSP TREE AND FIND WHICH NODES ARE BEHIND TARGET PLANE */
  /*====================================================================*/
  while (currentNode != NULL) {
    if (VF_VoxelBehindFacet(&(currentNode->bounds), facet)) {
      currentNode->behind = 1;
    } else {
      currentNode->behind = 0;
      if (currentNode->child[0] != NULL) {
        VF_BSP_PushStack(stack, currentNode->child[0]);
      }
      if (currentNode->child[1] != NULL) {
        VF_BSP_PushStack(stack, currentNode->child[1]);
      }
    }
    VF_BSP_PopStack(stack, currentNode);
  }
}

int VF_FindCandidates(Poly *poly_i, Poly *poly_j,
                      CandidateList **candidates)
{
  int         size, process, cnt, poly_cnt;
  int         ncandidates=0;
  double      spatial_tol;
  Point       center;
  Shaft       shaft;
  Poly        *poly;
  Facet       *facet;
  BinNodePtr  currentNode;
  BSPstackPtr stack;
  BinTree     *BSPtree;
  POLYstack   PolyStack;
  static      int list_size=0;
  static      CandidateList *list;
  VFtopology  *topology=VF_CurrentTopology();

  stack   = topology->BSP_stack;
  BSPtree = &(topology->BSP_tree);
  if (list_size==0) {
    list_size = LIST_INCR;
    size      = list_size*sizeof(CandidateList);
    list      = (CandidateList *)VF_Newv((unsigned)size);
  }
  VF_PolyCenter(poly_i, &center);
  VF_CreateShaft(&shaft, poly_i, poly_j);
  spatial_tol = topology->spatial_tol;
  VF_InitPolyStack(&PolyStack);

  /*=========================================================*/
  /* TRAVERSE THE BSP TREE AND FIND WHICH NODES INTERSECT    */
  /* THE CULLING SHAFT.  IF A LEAF NODE IS REACHED, TRAVERSE */
  /* ITS'S PATCHES TO FIND WHICH ONES GO INTO THE VISIBILITY */
  /* CANDIDATE LIST.                                         */
  /*=========================================================*/
  VF_BSP_InitStack(stack);
  currentNode = BSPtree->root;
  while (currentNode != NULL) {
    /*=================================================*/
    /* SKIP LEAF NODES THAT ARE NOT VISIBLE TO PATCH I */
    /*=================================================*/
    if (currentNode->behind) {
       VF_BSP_PopStack(stack, currentNode);
    /*=================================================*/
    /* SKIP LEAF NODES THAT ARE NOT VISIBLE TO PATCH J */
    /*=================================================*/
    } else if (VF_VoxelBehindPoly(&currentNode->bounds, poly_j)) {
      VF_BSP_PopStack(stack, currentNode);
    /*========================================*/
    /* SKIP LEAF NODES THAT DO NOT INTERSECT  */
    /* THE CULL SHAFT FROM PATCH I TO PATCH J */
    /*========================================*/
    } else if (!VF_ShaftCull(&currentNode->bounds, &shaft, spatial_tol)) {
      VF_BSP_PopStack(stack, currentNode);
    } else {
      if (currentNode->child[0] == NULL) {
        /*=======================================*/
        /* WE'RE AT A LEAF, PROCESS IT'S PATCHES */
        /*=======================================*/
        facet = VF_FirstFacetOfLinkList(&currentNode->members);
        while (facet) {
          /*================================*/
          /* CHECK IF PATCH INTERSECT SHAFT */
          /* FROM PATCH I TO PATCH J        */
          /*================================*/
          if (VF_ShaftCull(&(facet->bounds),&shaft, spatial_tol)) {
            /*=========================================*/
            /* CONVERT PATCH TO POLY(S) AND PERFORM    */
            /* MORE CHECKS OF EACH POLY FOR VISIBILITY */
            /*=========================================*/
            VF_FacetToPolyStack(facet, &PolyStack, &poly_cnt);
            for (cnt=0; cnt<poly_cnt; cnt++) {
              VF_PopPolyStackPtr(&PolyStack, poly);
              if (VF_BehindPoly(poly_i, poly)) {
                process = 0;
              } else if (VF_BehindPoly(poly_j, poly)) {
                process = 0;
              } else if (VF_BackFaceCullPolys(poly_i,poly)) {
                process = 0;
              } else {
                process = 1;
              }
              if (process) {
                /*====================*/
                /* ADD PATCH/POLY TO  */
                /* THE CANDIDATE LIST */
                /*====================*/
                if (ncandidates==list_size) {
                  list_size  = ncandidates;
                  list_size += LIST_INCR;
                  size       = list_size*sizeof(CandidateList);
                  list       = (CandidateList *)VF_ReNewv(list, size);
                }
                list[ncandidates].poly      = *poly;
                list[ncandidates].cached    = 0;
                list[ncandidates].check_it  = 1;
                list[ncandidates].facet_num = facet->index;
                list[ncandidates].poly_num  = PolyStack.cnt;
                ncandidates++;
              }
            }
          }
          facet = VF_NextFacetOfLinkList(&currentNode->members);
        }
        VF_BSP_PopStack(stack, currentNode);
      } else {
        VF_BSP_PushStack(stack, currentNode->child[1]);
        currentNode = currentNode->child[0];
      }
    }
  }
  *candidates = list;
  return (ncandidates);
}
