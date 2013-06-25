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
* 
*
*
*/

#ifndef TOUT_MESH_H
#define TOUT_MESH_H

#include <hdf5.h>

#include "danu_mesh_types.h"

/* Public defines */
#define MESH_ROOT_GROUP_NAME "Meshes"

#define MESH_NODAL_DATA_NAME    "Nodal Coordinates"
#define MESH_CONNECT_DATA_NAME "Element Connectivity"

#define MESH_NNODES_ATTR_NAME   "Number of Nodes"
#define MESH_NELEM_ATTR_NAME    "Number of Elements"
#define MESH_TYPE_ATTR_NAME     "Mesh Type"
#define MESH_DIM_ATTR_NAME      "Dimension"
#define MESH_ETYPE_ATTR_NAME    "Element Type"
#define MESH_EORDER_ATTR_NAME   "Element Order"

#define UNSTRUCT_MESH_NAME   "UNSTRUCTURED"
#define STRUCTURED_MESH_NAME "STRUCTURED"

/* Element string names */
#define LINE_ELEM_NAME "LINE"
#define TRI_ELEM_NAME  "TRI"
#define QUAD_ELEM_NAME "QUAD"
#define TET_ELEM_NAME  "TET"
#define HEX_ELEM_NAME  "HEX"

/* Useful error checking macros */
#define DIM_IS_VALID(a)         ( ( ( (a) > 0 ) && ( (a) < 4 ) )  ? 1 : 0 )
#define DIM_IS_INVALID(a)       ( ( ( (a) < 0 ) || ( (a) > 4 ) )  ? 1 : 0 )

#define EORDER_IS_VALID(a)       ( (  (a) == LINE_ELEM_ORDER || \
                                      (a) == TRI_ELEM_ORDER  || \
                                      (a) == QUAD_ELEM_ORDER || \
                                      (a) == TET_ELEM_ORDER  || \
                                      (a) == HEX_ELEM_ORDER   ) ? 1 : 0 )
#define EORDER_IS_INVALID(a)     ( ! EORDER_IS_VALID((a)) )

#define ETYPE_IS_VALID(a)        ( ( ( (a) > INVALID_ELEM ) && ( (a) <= HEX_ELEM ) ) ? 1 : 0 )
#define ETYPE_IS_INVALID(a)      ( ! ETYPE_IS_VALID((a)) )


/* Prototypes */
herr_t mesh_create_root_group(hid_t fid);
hid_t  mesh_open_root_group(hid_t fid);

herr_t mesh_count(hid_t fid, int *nmeshes);
herr_t mesh_list(hid_t fid, int num, char **meshnames, int *num_found);
herr_t mesh_exists(hid_t fid, const char * meshname, int * exist); 

herr_t mesh_create(hid_t fid, const char *mesh, tmesh_t mtype, telem_t etype, hid_t *mid);
hid_t  mesh_open(hid_t fid, const char *mesh);

herr_t mesh_add_unstructured(hid_t fid, const char *name, int elemorder, int dim, hid_t *mid);

/* Mesh Coordinates */
herr_t mesh_create_coordinates(hid_t mid, int dim, int nnodes);
hid_t  mesh_open_coordinates(hid_t mid);

herr_t mesh_write_coordinates(hid_t mid, int nnodes, double * x, double * y, double * z);
herr_t mesh_write_coordinates_1d(hid_t mid, int nnodes, double * x );
herr_t mesh_write_coordinates_2d(hid_t mid, int nnodes, double * x, double *y );

herr_t mesh_read_coordinates(hid_t mid, double * x, double *y, double *z);
herr_t mesh_read_coordinates_byindex(hid_t mid, int idx, double *buff);
herr_t mesh_read_coordinates_1d(hid_t mid, double * x);
herr_t mesh_read_coordinates_2d(hid_t mid, double * x, double *y);

/* Mesh Connectivity */
hid_t  mesh_create_connectivity(hid_t mid, int nelem);
hid_t  mesh_open_connectivity(hid_t mid);
herr_t mesh_write_connectivity(hid_t mid, int nelem, const int *data);
herr_t mesh_read_connectivity(hid_t mid,  int *data);
herr_t mesh_connectivity_size(hid_t mid, int *nelem, int *elem_order);

/* Mesh Attributes */
herr_t mesh_get_type(hid_t mid, tmesh_t *type);
herr_t mesh_get_elementtype(hid_t mid, telem_t *type);
herr_t mesh_get_dimension(hid_t mid, int *dim);
herr_t mesh_get_nnodes(hid_t mid, int *nnodes);
herr_t mesh_get_nelem(hid_t mid, int *nelem);
herr_t mesh_get_elem_order(hid_t mid, int *order);

/* Mesh Wrapper functions */
hid_t mesh_create_hex_unstruct(hid_t fid, 
                               const char * mesh_name,
                               int nnodes,
                               double *x,
                               double *y,
                               double *z,
                               int nelem,
                               int *conn);

#endif
