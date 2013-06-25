/* *************************************************************************** *
*                                                                              *
*                                                                              *
*                             Copyright  (C) 20xx,                             *
*                      Los Alamos National Security, LLC                       *
*                                                                              *
*                             LA-CC-xxxxxx                                     *
*                                                                              *
* **************************************************************************** */

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

  Output *file;       /* File where mesh is located */

  tmesh_t  mesh_type;  /* Mesh type INVALID, UNSTRUCTURED, STRUCTURED */

  telem_t  elem_type;  /* Mesh element type INVALID, LINE, TRI, QUAD, TET, HEX */
  
  int dim;             /* Mesh dimension */

} Mesh;


#endif

