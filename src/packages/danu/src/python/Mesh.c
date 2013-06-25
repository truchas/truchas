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


Mesh * allocate_mesh_object(Output *fh, hid_t mid, const char *meshname, tmesh_t mesh_type, telem_t elem_type )
{
    Mesh *m = DANU_MALLOC(Mesh,1);
    int dim;
    if ( m != NULL ) {
	m->file = allocate_output_object(GET_OUTPUT_OBJ_HID(fh),fh->h5->name,fh->access);
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
    deallocate_output_object(m->file);
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
      if ( H5_ISA_VALID_ID(mid) ) {
	  err = mesh_get_type(mid,&mesh_type);
	  err = mesh_get_elementtype(mid,&elem_type);
      }
  }
  else {
      err = mesh_create(fid, meshname,mesh_type,elem_type,&mid);
  }

  if ( H5_ISA_VALID_ID(mid) ) {
    m = allocate_mesh_object(fh,mid,meshname,mesh_type,elem_type);
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
