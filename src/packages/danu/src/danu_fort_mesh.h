/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* fort_mesh.h
*
*  DANU FORTRAN/C interfaces for meshs
*
*/

#ifndef DANU_FORT_MESH_H
#define DANU_FORT_MESH_H

#include <danu_fort_hid.h>

/* Use these defines to improve the readability of the code */
/* Basic mesh control */
#if 0
#include "Fortran2C.h"
#define mesh_exists_f            FORTRAN_FUNC_GLOBAL_(mesh_exists_f, MESH_EXISTS_F)
#define mesh_count_f             FORTRAN_FUNC_GLOBAL_(mesh_count_f, MESH_COUNT_F)
#define mesh_list_f              FORTRAN_FUNC_GLOBAL_(mesh_list_f, MESH_LIST_F)
#define mesh_add_unstructured_f  FORTRAN_FUNC_GLOBAL_(mesh_add_unstructured_f, MESH_ADD_UNSTRUCTURED_F)
#define mesh_write_coordinates_f FORTRAN_FUNC_GLOBAL_(mesh_write_coordinates_f, MESH_WRITE_COORDINATES_F)
#endif


/* Function prototypes corresponding interface definition in module_iface.f90 */
void mesh_exists_f(const hid_t_ptr *fid, const char *meshname, const int *flen, int *flag, int *ierr);
void mesh_count_f(const hid_t_ptr *fid,  int *count, int *ierr);
void mesh_names_f(const hid_t_ptr *fid,  char *names, const int *flen, const int *num, int *ierr);

void mesh_add_unstructured_f(const hid_t_ptr *fid, 
                             const char *name,
                             int *flen,
                             int *elemorder,
                             int *dim,
                             hid_t_ptr *mid,
                             int *ierr);

/* Coordinates */
void mesh_write_coordinates_1d_f(const hid_t_ptr *mptr, int *nnodes, double *x, int *ierr);
void mesh_write_coordinates_2d_f(const hid_t_ptr *mptr, int *nnodes, double *x, double *y, int *ierr);
void mesh_write_coordinates_f(const hid_t_ptr *mptr, int *nnodes, double *x, double *y, double *z, int *ierr);

void mesh_read_coordinates_1d_f(const hid_t_ptr *mptr, double *x, int *ierr);
void mesh_read_coordinates_2d_f(const hid_t_ptr *mptr, double *x, double *y, int *ierr);
void mesh_read_coordinates_f(const hid_t_ptr *mptr, double *x, double *y, double *z, int *ierr);
void mesh_read_coordiantes_byindex_f(const hid_t_ptr *mptr, const int *fidx, double *buf, int *ierr); 

/* Connectivity */
void mesh_write_connectivity_f(const hid_t_ptr *mptr, const int *nelem, const int *data, int *ierr);
void mesh_read_connectivity_f(const hid_t_ptr *mptr, int *data, int *ierr);

/* Attributes */
void mesh_get_type_f(const hid_t_ptr *mptr, tmesh_t *type, int *ierr);
void mesh_get_elementtype_f(const hid_t_ptr *mptr, telem_t *elem_type, int *ierr);
void mesh_get_dimension_f(const hid_t_ptr *mptr, int *dim, int *ierr);
void mesh_get_nnodes_f(const hid_t_ptr *mptr, int *nnodes, int *ierr);
void mesh_get_nelem_f(const hid_t_ptr *mptr, int *nelem, int *ierr);
void mesh_get_elem_order_f(const hid_t_ptr *mptr, int *elem_order, int *ierr);

/* Wrappers */
void mesh_create_hex_unstruct_f(const hid_t_ptr *fptr,
                                const char *mname,
                                int *flen,
                                int *nnodes,
                                double *x,
                                double *y,
                                double *z,
                                int *nelem,
                                int *conn,
                                hid_t_ptr *mesh,
                                int *ierr);


                              

#endif /* DANU_FORT_MESH_H */

