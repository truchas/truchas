/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/


#ifndef IO_SCORPIO_H
#define IO_SCORPIO_H

#include "scorpio.h"

#include <danu_fort_hid.h>

struct TruchasScorpioFileHandle {
  int fhandle; // Scorpio file handle
  iogroup_t myIOgroup; // Internal Scorpio data
};

struct TruchasScorpioFileHandle * truchas_scorpio_create_handle(
        const char *filename, int numIOgroups);

void truchas_scorpio_free_handle(struct TruchasScorpioFileHandle *h);

void truchas_scorpio_write_dataset_1d_integer(
        struct TruchasScorpioFileHandle *h,
        const char *name,
        int *vector, int global_dim, int local_dim);

void truchas_scorpio_write_dataset_1d_double(struct TruchasScorpioFileHandle *h,
        const char *name,
        double *vector, int global_dim, int local_dim);

int truchas_scorpio_get_group_name(hid_t_ptr gid, char **cname);

void truchas_scorpio_str_free(char *cname);

hid_t truchas_scorpio_get_hdf5_handle(struct TruchasScorpioFileHandle *h);

#endif
