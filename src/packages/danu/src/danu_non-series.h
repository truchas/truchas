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
* danu_sim.h
*
*  Output Simulation Non-series data
*
*/

#ifndef TOUT_NONSERIES_H
#define TOUT_NONSERIES_H

#include <hdf5.h>

/* Public define's */
#define NON_SERIES_GROUP_NAME "Non-series Data"

/* Public Functions */
herr_t data_create_group(hid_t sid);
hid_t  data_open_group(hid_t sid);

hid_t  data_open_dataset(hid_t sid, const char *name);

herr_t data_exists(hid_t sid, const char *dataname, int *exists);

herr_t data_count(hid_t sid, int *ndata);
herr_t data_rank(hid_t sid, const char *name, int *rank);
herr_t data_dimensions(hid_t sid, const char *name, int rank, int *dimensions);
herr_t data_type(hid_t sid, const char *name, int *type);
herr_t data_list(hid_t sid, int num, char **datanames, int * num_found);

herr_t data_write_int(hid_t sid, const char * name, int dim, const int * dimensions, int *data);
herr_t data_write_float(hid_t sid, const char * name, int dim, const int * dimensions, float *data);
herr_t data_write_double(hid_t sid, const char * name, int dim, const int * dimensions, double *data);
herr_t data_write_strings(hid_t sid, const char * name, int dim, const int * dimensions, char *data);

herr_t data_read_int(hid_t sid, const char * name, int dim, const int * dimensions, int *data);
herr_t data_read_float(hid_t sid, const char * name, int dim, const int * dimensions, float *data);
herr_t data_read_double(hid_t sid, const char * name, int dim, const int * dimensions, double *data);
herr_t data_read_strings(hid_t sid, const char * name, int dim, const int * dimensions, char *data);



#endif
