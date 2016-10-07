/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * Danu Output Objects
 *
 * Python Interface for Output Objects
 *
 */

/* System includes */

#include <hdf5.h>
#include <Python.h>

/* Danu includes */
#include <danu.h>

#include "DanuH5Object.h"
#include "DanuOutput.h"

#include "DanuMesh.h"


Mesh * allocate_mesh_object(hid_t mid, const char *meshname, tmesh_t mesh_type, telem_t elem_type )
{
    Mesh *m = DANU_MALLOC(Mesh,1);
    int dim;
    if ( m != NULL ) {
        m->h5   = allocate_h5_object(mid,meshname);
	m->mesh_type = mesh_type;
	m->elem_type = elem_type;
	if ( H5_RETURN_OK(mesh_get_dimension(mid,&dim))) {
	    m->dim = dim;
	}
    }
    return m;
}

void deallocate_mesh_object(Mesh *m)
{
    deallocate_h5_object(m->h5);
    DANU_FREE(m);
}

Mesh * construct_mesh_object(Output *fh, const char *meshname, tmesh_t mesh_type, telem_t elem_type )
{
  hid_t fid = GET_OUTPUT_OBJ_HID(fh);
  herr_t err;
  Mesh *m = NULL;
  hid_t mid;
  int exists;

  err = mesh_exists(fid,meshname,&exists);

  if ( exists ) {
    mid = mesh_open(fid,meshname);
    m = construct_mesh_object_by_id(mid);
  }
  else {
    err = mesh_create(fid, meshname,mesh_type,elem_type,&mid);
    if ( H5_ISA_VALID_ID(mid) ) {
      m = allocate_mesh_object(mid,meshname,mesh_type,elem_type);
    }
  }

  return m;
}

Mesh * construct_mesh_object_by_id(hid_t mid)
{
  Mesh *m = NULL;
  char meshname[128];
  herr_t err1, err2;
  tmesh_t m_type;
  telem_t m_elem;


  if ( H5_ISA_VALID_ID(mid) ) {
    err1 = mesh_get_type(mid,&m_type);
    err2 = mesh_get_elementtype(mid,&m_elem);
    if ( DANU_RETURN_SELECT(err1,err2) == DANU_SUCCESS ) {
      H5Iget_name(mid,meshname,128);
      m = allocate_mesh_object(mid,meshname,m_type,m_elem);
    }
    else {
      if ( err1 != DANU_SUCCESS ) {
        throw_exception("Failed to determine mesh type");
      }
      if ( err2 != DANU_SUCCESS ) {
        throw_exception("Failed to determine mesh element type");
      }
    }
  }
  else {
    throw_exception("Invalid HDF5 id");
  }

  return m;

}

void deconstruct_mesh_object(Mesh *m)
{
    hid_t mid = GET_MOBJECT_HID(m);
    if ( H5_ISA_VALID_ID(mid) ) {
        danu_group_close(mid);
    }
    deallocate_mesh_object(m);
}
