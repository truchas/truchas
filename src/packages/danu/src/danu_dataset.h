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
* danu_dataset.h
*
*  DANU Dataset module
*
*/

#ifndef DANU_DATASET_H
#define DANU_DATASET_H

#include <hdf5.h>

#include <danu_slab.h>

hid_t  danu_dataset_create(hid_t loc, const char * name, hid_t type, int dim, const hsize_t * size, dbool_t extend);
hid_t  danu_dataset_open(hid_t loc, const char * name);
herr_t danu_dataset_close(hid_t id);

herr_t danu_dataset_read(hid_t id, dslab_t *slab, hid_t buf_type, int buf_dim, const hsize_t *buf_size, void *buf);   
herr_t danu_dataset_write(hid_t id, dslab_t *slab, hid_t buf_type,int buf_dim, const hsize_t * buf_size, const void *buf);

hbool_t danu_dataset_exists(hid_t loc, const char * name);

int     danu_dataset_rank(hid_t loc, const char *name);
herr_t  danu_dataset_dimensions(hid_t loc, const char *name, int size, hsize_t * dimensions);
herr_t  danu_dataset_max_dimensions(hid_t loc, const char *name, int size, hsize_t * max_dimensions);
herr_t  danu_dataset_type(hid_t loc, const char *dataname, int *typecode);

int     danu_dataset_rank2(hid_t id);
herr_t  danu_dataset_dimensions2(hid_t id, int size, hsize_t * dimensions);
herr_t  danu_dataset_max_dimensions2(hid_t id, int size, hsize_t * max_dimensions);
herr_t  danu_dataset_type2(hid_t id, int *typecode);

hid_t danu_dataset_create_double(hid_t loc, const char * name, int dim, const hsize_t*size);
hid_t danu_dataset_create_float(hid_t loc, const char * name, int dim, const hsize_t*size);
hid_t danu_dataset_create_int(hid_t loc, const char * name, int dim, const hsize_t*size);
hid_t danu_dataset_create_string(hid_t loc, const char * name, int dim, const hsize_t*size);

herr_t danu_data_write_double(hid_t loc, const char * name, int dim, const hsize_t *size, const double * buf);
herr_t danu_data_write_float(hid_t loc, const char * name, int dim, const hsize_t *size, const float * buf);
herr_t danu_data_write_int(hid_t loc, const char * name, int dim, const hsize_t *size, const int * buf);
herr_t danu_data_write_strings(hid_t loc, const char * name, int num, const char ** buf);

herr_t danu_data_write_double2(hid_t id, int dim, const hsize_t *size, const double * buf);
herr_t danu_data_write_float2(hid_t id,  int dim, const hsize_t *size, const float * buf);
herr_t danu_data_write_int2(hid_t id,    int dim, const hsize_t *size, const int * buf);
herr_t danu_data_write_strings2(hid_t id, int num, const char ** buf);

herr_t danu_data_read_double(hid_t loc, const char * name, int dim, const hsize_t *size, double * buf);
herr_t danu_data_read_float(hid_t loc, const char * name, int dim, const hsize_t *size, float * buf);
herr_t danu_data_read_int(hid_t loc, const char * name, int dim, const hsize_t *size, int * buf);
herr_t danu_data_read_strings(hid_t loc, const char * name, int num, char ** buf);

herr_t danu_data_read_double2(hid_t id, int dim, const hsize_t *size, double * buf);
herr_t danu_data_read_float2(hid_t id, int dim, const hsize_t *size, float * buf);
herr_t danu_data_read_int2(hid_t id, int dim, const hsize_t *size, int * buf);
herr_t danu_data_read_strings2(hid_t id, int num, char ** buf);

#endif

