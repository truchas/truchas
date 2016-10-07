/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_series-dataset.h
*
*  Telluride Output Simulation Module 
*
*/

#ifndef TOUT_SERIES_DATA_H
#define TOUT_SERIES_DATA_H

#include <hdf5.h>

/* Public define's */

/* Public Functions */

hid_t simulation_open_data(hid_t nsid, const char *dataname);

herr_t simulation_data_exists(hid_t nsid, const char * dataname, int *exists);

int    simulation_data_rank(hid_t nsid, const char * dataname);

herr_t simulation_data_dimensions(hid_t nsid, 
                                  const char * dataname, 
                                  int size, hsize_t *dimensions);

herr_t simulation_data_type(hid_t nsid, const char *dataname, int *typecode);


herr_t simulation_data_count(hid_t nsid, int *nseriesdata);
herr_t simulation_data_list(hid_t nsid, int num,
                             char **datanames, int *num_found);

herr_t simulation_data_write_int(hid_t nsid, const char * dataname, 
                                  int dim, const int *size, const int *data);

herr_t simulation_data_read_int(hid_t nsid, const char * dataname, 
                                  int dim, const int *size, int *buffer);

herr_t simulation_data_write_double(hid_t nsid, const char * dataname, 
                                    int dim, const int *size, const double *data);

herr_t simulation_data_read_double(hid_t nsid, const char * dataname, 
                                   int dim, const int *size, double *buffer);

herr_t simulation_data_write_float(hid_t nsid, const char * dataname, 
                                    int dim, const int *size, const float *data);

herr_t simulation_data_read_float(hid_t nsid, const char * dataname, 
                                   int dim, const int *size, float *buffer);
#endif
