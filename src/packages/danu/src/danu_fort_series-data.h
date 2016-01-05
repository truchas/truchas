/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
*
*  DANU Fortran/C Series Data Interface
*
*/

#ifndef DANU_FORT_SERIES_DATA_H
#define DANU_FORT_SERIES_DATA_H

#include <danu_fort_hid.h>

// Query subroutines
void simulation_data_exists_f(const hid_t_ptr *nptr,
                              const char * dataname,
                              const int *len,
                              int *flag,
                              int *ierr);

void simulation_data_count_f(const hid_t_ptr *nptr,
                             int * count,
                             int *ierr);

void simulation_data_rank_f(const hid_t_ptr *nptr,
                          const char *dataname,
                          const int *flen,
                          int *rank,
                          int *ierr);

void simulation_data_dimensions_f(const hid_t_ptr *nptr,
                                const char *dataname,
                                const int *flen,
                                const int *size,
                                int *dimensions,
                                int *ierr);

void simulation_data_type_f(const hid_t_ptr *nptr,
                            const char *dataname,
                            const int *flen,
                            int *typecode,
                            int *ierr);

void simulation_data_write_int0_f(const hid_t_ptr * nptr,
                                  const char *data_name,
                                  const int *flen,
                                  const int *data,
                                  int *ierr);

void simulation_data_write_int_f(const hid_t_ptr *nptr,
                                 const char *data_name,
                                 const int *flen,
                                 const int *dim,
                                 const int *dimensions,
                                 const int *data,
                                 int *ierr);

void simulation_data_write_float0_f(const hid_t_ptr * nptr,
                                    const char *data_name,
                                    const int *flen,
                                    const float *data,
                                    int *ierr);

void simulation_data_write_float_f(const hid_t_ptr *nptr,
                                    const char *data_name,
                                    const int *flen,
                                    const int *dim,
                                    const int *dimensions,
                                    const float *data,
                                    int *ierr);

void simulation_data_write_double0_f(const hid_t_ptr * nptr,
                                    const char *data_name,
                                    const int *flen,
                                    const double *data,
                                    int *ierr);

void simulation_data_write_double_f(const hid_t_ptr *nptr,
                                    const char *data_name,
                                    const int *flen,
                                    const int *dim,
                                    const int *dimensions,
                                    const double *data,
                                    int *ierr);

void simulation_data_read_int0_f(const hid_t_ptr * nptr,
                                 const char *data_name,
                                 const int *flen,
                                 int *data,
                                 int *ierr);

void simulation_data_read_int_f(const hid_t_ptr *nptr,
                                const char *data_name,
                                const int *flen,
                                const int *dim,
                                const int *dimensions,
                                int *data,
                                int *ierr);

void simulation_data_read_float0_f(const hid_t_ptr * nptr,
                                   const char *data_name,
                                   const int *flen,
                                   float *data,
                                   int *ierr);

void simulation_data_read_float_f(const hid_t_ptr *nptr,
                                  const char *data_name,
                                  const int *flen,
                                  const int *dim,
                                  const int *dimensions,
                                  float *data,
                                  int *ierr);

void simulation_data_read_double0_f(const hid_t_ptr * nptr,
                                   const char *data_name,
                                   const int *flen,
                                   double *data,
                                   int *ierr);

void simulation_data_read_double_f(const hid_t_ptr *nptr,
                                  const char *data_name,
                                  const int *flen,
                                  const int *dim,
                                  const int *dimensions,
                                  double *data,
                                  int *ierr);


#endif

