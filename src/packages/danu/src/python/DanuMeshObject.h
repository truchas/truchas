/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * Danu Mesh Objects
 *
 * Python Interface Ouput Objects
 *
 */

#ifndef DANU_PY_MESH_OBJ_H
#define DANU_PY_MESH_OBJ_H

#include "DanuH5Object.h"
#include "DanuOutputObject.h"


typedef struct {

  H5Obj  *h5;          /* Internal HDF5 Identifier Object */

  tmesh_t  mesh_type;  /* Mesh type INVALID, UNSTRUCTURED, STRUCTURED */

  telem_t  elem_type;  /* Mesh element type INVALID, LINE, TRI, QUAD, TET, HEX */
  
  int dim;             /* Mesh dimension */

} Mesh;


#endif

