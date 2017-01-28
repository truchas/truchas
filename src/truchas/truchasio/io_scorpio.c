/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

#include <io_scorpio.h>

struct TruchasScorpioFileHandle * truchas_scorpio_create_handle(
    const char *filename, int numIOgroups)
{
  struct TruchasScorpioFileHandle *h
      = malloc(sizeof(struct TruchasScorpioFileHandle));
  iogroup_conf_t myIOconfig;
  int i;

  myIOconfig.numIOgroups = numIOgroups;
  myIOconfig.commIncoming = MPI_COMM_WORLD;

  scorpio_IOgroup_init(&myIOconfig, &h->myIOgroup);

  h->fhandle = scorpio_open_file(filename, &h->myIOgroup,
      SCORPIO_FILE_CREATE);
  if (h->fhandle == FAILURE) {
    fprintf(stderr, "Opening file for writing failed.\n");
    exit(-1);
  }
  return h;
}

void truchas_scorpio_free_handle(struct TruchasScorpioFileHandle *h)
{
  scorpio_close_file(h->fhandle, &h->myIOgroup);
  scorpio_IOgroup_cleanup(&h->myIOgroup);
  free(h);
}

void truchas_scorpio_write_dataset_1d_integer(
        struct TruchasScorpioFileHandle *h,
        const char *name,
        int *vector, int global_dim, int local_dim)
{
#define NDIMS 1
  int ndims = NDIMS;
  int globaldims[NDIMS];
  int localdims[NDIMS];
  globaldims[0] = global_dim;
  localdims[0] = local_dim;
  scorpio_write_dataset(vector, SCORPIO_INTEGER, ndims, globaldims, localdims,
      h->fhandle, (char *)name, &h->myIOgroup,
      SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE);
}

void truchas_scorpio_write_dataset_1d_double(struct TruchasScorpioFileHandle *h,
        const char *name,
        double *vector, int global_dim, int local_dim)
{
  const int ndims = 1;
  int globaldims[ndims];
  int localdims[ndims];
  globaldims[0] = global_dim;
  localdims[0] = local_dim;
  scorpio_write_dataset(vector, SCORPIO_DOUBLE, ndims, globaldims, localdims,
      h->fhandle, (char *)name, &h->myIOgroup,
      SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE);
}

void truchas_scorpio_write_dataset_2d_double(struct TruchasScorpioFileHandle *h,
        const char *name,
        double *vector, int global_dim, int local_dim1, int local_dim2)
{
  const int ndims = 2;
  int globaldims[ndims];
  int localdims[ndims];
  globaldims[0] = global_dim;
  globaldims[1] = local_dim1;
  localdims[0] = local_dim2;
  localdims[1] = local_dim1;
  scorpio_write_dataset(vector, SCORPIO_DOUBLE, ndims, globaldims, localdims,
      h->fhandle, (char *)name, &h->myIOgroup,
      SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE);
}

void truchas_scorpio_write_dataset_2d_byte(struct TruchasScorpioFileHandle *h,
        const char *name,
        char *vector, int global_dim, int local_dim1, int local_dim2)
{
  const int ndims = 2;
  int globaldims[ndims];
  int localdims[ndims];
  globaldims[0] = global_dim;
  globaldims[1] = local_dim1;
  localdims[0] = local_dim2;
  localdims[1] = local_dim1;
  scorpio_write_dataset(vector, SCORPIO_BYTE, ndims, globaldims, localdims,
      h->fhandle, (char *)name, &h->myIOgroup,
      SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE);
}

hid_t truchas_scorpio_get_hdf5_handle(struct TruchasScorpioFileHandle *h)
{
  if (h->myIOgroup.localrank != 0) {
    fprintf(stderr, "h->fid is only initialized on IO nodes.\n");
    exit(-1);
  }
  return h->myIOgroup.file[h->fhandle]->fid;
}

int truchas_scorpio_get_group_name(hid_t_ptr gid, char **cname)
{
  ssize_t size;
  size = 1 + H5Iget_name(gid->id, NULL, 0);
  *cname = (char *) malloc(size);
  size = H5Iget_name(gid->id, *cname, size);
  return size;
}

void truchas_scorpio_str_free(char *cname)
{
  free(cname);
}

// Danu helper functions

#include <danu_file.h>

herr_t output_file_initial_setup(hid_t fid);

void truchas_scorpio_handle_initial_setup(struct
    TruchasScorpioFileHandle *h)
{
  output_file_initial_setup(truchas_scorpio_get_hdf5_handle(h));
}

void truchas_scorpio_handle_file_close(struct TruchasScorpioFileHandle *h)
{
  danu_file_close(truchas_scorpio_get_hdf5_handle(h));
}

hid_t_ptr truchas_scorpio_hdf5_handle_danu_create(struct
    TruchasScorpioFileHandle *h)
{
  return create_hid_struct(truchas_scorpio_get_hdf5_handle(h));
}

void truchas_scorpio_hdf5_handle_danu_free(hid_t_ptr fid_ptr)
{
  destroy_hid_struct(fid_ptr);
}
