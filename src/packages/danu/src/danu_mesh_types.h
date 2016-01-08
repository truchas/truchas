/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* 
*
*
*/

#ifndef TOUT_MESH_TYPE_H
#define TOUT_MESH_TYPE_H

/* Mesh types */
typedef enum tmesh_t {
    INVALID_MESH = -3,
    UNSTRUCTURED_MESH,
    STRUCTURED_MESH
} tmesh_t;

#define UNSTRUCT_MESH_NAME   "UNSTRUCTURED"
#define STRUCTURED_MESH_NAME "STRUCTURED"

/* Element types */
typedef enum telem_t {
    INVALID_ELEM = 0,
    LINE_ELEM,
    TRI_ELEM,
    QUAD_ELEM,
    TET_ELEM,
    HEX_ELEM,
} telem_t;

/* Element Orders (number of vertices) */
#define LINE_ELEM_ORDER     2
#define TRI_ELEM_ORDER      3
#define QUAD_ELEM_ORDER     4
#define TET_ELEM_ORDER      4
#define HEX_ELEM_ORDER      8

#endif
