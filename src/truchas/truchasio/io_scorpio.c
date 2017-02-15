/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

#include <io_scorpio.h>
#include <assert.h>

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

void truchas_scorpio_write_attr_0d_integer(struct TruchasScorpioFileHandle *h,
        char *attr_name, int *attr_data, char *obj_name)
{
  scorpio_write_simple_attr(attr_name, attr_data, SCORPIO_INTEGER, h->fhandle, obj_name, &h->myIOgroup);
}

void truchas_scorpio_write_attr_0d_double(struct TruchasScorpioFileHandle *h,
        char *attr_name, double *attr_data, char *obj_name)
{
  scorpio_write_simple_attr(attr_name, attr_data, SCORPIO_DOUBLE, h->fhandle, obj_name, &h->myIOgroup);
}

void truchas_scorpio_write_attr_1d_double(struct TruchasScorpioFileHandle *h,
        char *attr_name, double *attr_data, int ndims, int *adims, char *obj_name)
{
  scorpio_write_attr(attr_name, attr_data, SCORPIO_DOUBLE, ndims, adims, h->fhandle, obj_name, &h->myIOgroup);
}

void truchas_scorpio_write_attr_0d_string(struct TruchasScorpioFileHandle *h,
        char *attr_name, char *attr_data, char *obj_name)
{
  scorpio_write_simple_attr(attr_name, attr_data, SCORPIO_STRING, h->fhandle, obj_name, &h->myIOgroup);
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

void truchas_scorpio_write_dataset_2d_integer(
        struct TruchasScorpioFileHandle *h,
        const char *name,
        int *vector, int global_dim, int local_dim1, int local_dim2)
{
  const int ndims = 2;
  int globaldims[ndims];
  int localdims[ndims];
  globaldims[0] = global_dim;
  globaldims[1] = local_dim1;
  localdims[0] = local_dim2;
  localdims[1] = local_dim1;
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

int truchas_scorpio_create_dataset_group(struct TruchasScorpioFileHandle *h, char *group_name)
{
  return scorpio_create_dataset_group(group_name, h->fhandle, &h->myIOgroup);
  /* returns the group id, but only on the ranks that write; other ranks return 0.
     Failure is captured by an internal assert */
  /* scorpio_close_dataset_group(group_id, h->fhandle, &h->myIOgroup); */
}

void truchas_scorpio_close_dataset_group(struct TruchasScorpioFileHandle *h, int gid)
{
  scorpio_close_dataset_group(gid, h->fhandle, &h->myIOgroup);
  /* ignores the return value (==0 success, !=0 error) which is not consistent across ranks */
}

void truchas_scorpio_create_link(struct TruchasScorpioFileHandle *h, char *target, int link_loc_id, char *link_name)
{
  scorpio_create_link(target, link_loc_id, link_name, h->fhandle, &h->myIOgroup);
}

void truchas_scorpio_handle_file_close(struct TruchasScorpioFileHandle *h)
{
  hid_t fid, *objects;
  ssize_t num_obj;
  int i;
  
  if (h->myIOgroup.localrank == 0) {
    fid = truchas_scorpio_get_hdf5_handle(h);

    /* Close all the datasets */
    num_obj = H5Fget_obj_count(fid,H5F_OBJ_DATASET);
    objects = (hid_t *) malloc(num_obj*sizeof(hid_t));
    if ( num_obj > 0 ) {
        H5Fget_obj_ids(fid,H5F_OBJ_DATASET,num_obj,objects);
        for(i=0;i<num_obj;i++)
            H5Dclose(objects[i]);
    }
    free(objects);

    /* Close all the groups */
    num_obj = H5Fget_obj_count(fid,H5F_OBJ_GROUP);
    objects = (hid_t *) malloc(num_obj*sizeof(hid_t));
    if ( num_obj > 0 ) {
        H5Fget_obj_ids(fid,H5F_OBJ_GROUP,num_obj,objects);
        for(i=0;i<num_obj;i++)
            H5Gclose(objects[i]);
    }
    free(objects);
  }
  scorpio_close_file(h->fhandle, &h->myIOgroup);
  scorpio_IOgroup_cleanup(&h->myIOgroup);
  free(h);
}

/* PROBE OUTPUT
 * This adapts Danu serial code from probes.c to the parallel structure of
 * scorpio_write_dataset.
 */
/*==============================================================================

Truchas probe data is migrated to the IO process (rank 0).  In fact I believe
it is replicated on all processes.  Danu would be called from rank 0 to write
the data.  Here we
==============================================================================*/

herr_t scorpio_write_probe_data(hid_t pid, hid_t hdf_type, int rank, int *dims, void *data, int fhandle, iogroup_t *myIOgroup)
{
  int ndims, ret, i;
  hsize_t *size, *count, *offset;
  hid_t sid, file_dataspace, mem_dataspace, xfer_plist;

  if (myIOgroup->localrank == 0) {

    sid = H5Dget_space(pid);
    assert(sid != FAILURE);
    ndims = H5Sget_simple_extent_ndims(sid);
    assert(rank == ndims);
    size = (hsize_t *) calloc(ndims, sizeof(hsize_t));
    H5Sget_simple_extent_dims(sid, size, NULL);
    H5Sclose(sid);
    for (i = 1; i < ndims; i++)
      assert(dims[i] == size[i]);

    count = (hsize_t *) calloc(ndims, sizeof(hsize_t));
    for (i = 0; i < ndims; i++)
      count[i] = dims[i];
    if (myIOgroup->globalrank != 0) count[0] = 0;

    offset = (hsize_t *) calloc(ndims, sizeof(hsize_t));
    offset[0] = size[0];
    for (i = 1; i < ndims; i++)
      offset[i] = 0;

    size[0] += dims[0];
    ret = H5Dextend(pid, size);
    assert(ret != FAILURE);

    file_dataspace = H5Dget_space(pid);
    ret = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    assert(ret != FAILURE);

    mem_dataspace = H5Screate_simple(ndims, count, NULL);

    xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    assert(xfer_plist != FAILURE);

#if !defined(SERIAL_HDF5)
    ret = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    assert(ret != FAILURE);
#endif

    ret = H5Dwrite(pid, hdf_type, mem_dataspace, file_dataspace, xfer_plist, data);

    H5Pclose(xfer_plist);
    H5Sclose(mem_dataspace);
    H5Sclose(file_dataspace);

    free(offset);
    free(count);
    free(size);
  }
  return 0;
}

void truchas_scorpio_write_probe_data_2d_double(struct TruchasScorpioFileHandle *h,
    hid_t pid, int *dims, double *data, int fhandle, iogroup_t *myIOgroup)
{
  int cdims[2];
  cdims[0] = dims[1];
  cdims[1] = dims[0];
  scorpio_write_probe_data(pid, H5T_NATIVE_DOUBLE, 2, cdims, data, h->fhandle, &h->myIOgroup);
}

hid_t scorpio_create_probe(hid_t sid, char *name, hid_t hdf_type, int rank, int *dims, void *data, int fhandle, iogroup_t *myIOgroup)
{
  int ret, i;
  hsize_t *size, *max_size, *chunk_size, *count, *start;
  hid_t pid, gid, space_id, file_dataspace, mem_dataspace, link_plist, xfer_plist;

  if (myIOgroup->localrank == 0) {

    gid = H5Gopen(sid, "Probes", H5P_DEFAULT);

    size = (hsize_t *) calloc(rank, sizeof(hsize_t));
    for (i = 0; i < rank; i++)
      size[i] = dims[i];

    max_size = (hsize_t *) calloc(rank, sizeof(hsize_t));
    max_size[0] = H5S_UNLIMITED;
    for (i = 1; i < rank; i++)
      max_size[i] = dims[i];

    chunk_size = (hsize_t *) calloc(rank, sizeof(hsize_t));
    chunk_size[0] = 1;
    for (i = 1; i < rank; i++)
      chunk_size[i] = dims[i];

    space_id = H5Screate_simple(rank, size, max_size);

    link_plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(link_plist, rank, chunk_size);
    pid = H5Dcreate2(gid, name, hdf_type, space_id, H5P_DEFAULT, link_plist, H5P_DEFAULT);
    H5Sclose(space_id);


    count = (hsize_t *) calloc(rank, sizeof(hsize_t));
    start = (hsize_t *) calloc(rank, sizeof(hsize_t));
    for (i = 0; i < rank; i++) {
      count[i] = dims[i];
      start[i] = 0;
    }
    if (myIOgroup->globalrank != 0) count[0] = 0;

    file_dataspace = H5Dget_space(pid);
    ret = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
    assert(ret != FAILURE);

    mem_dataspace = H5Screate_simple(rank, count, NULL);

    xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    assert(xfer_plist != FAILURE);

#if !defined(SERIAL_HDF5)
    ret = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    assert(ret != FAILURE);
#endif

    ret = H5Dwrite(pid, hdf_type, mem_dataspace, file_dataspace, xfer_plist, data);
    assert(ret != FAILURE);

    H5Pclose(xfer_plist);
    H5Sclose(mem_dataspace);
    H5Sclose(file_dataspace);
    H5Pclose(link_plist);
    H5Gclose(gid);

    free(size);
    free(max_size);
    free(chunk_size);
    free(start);
    free(count);

    return pid;
  } else {
    return H5I_INVALID_HID;
  }
}

hid_t truchas_scorpio_create_probe_2d_double(struct TruchasScorpioFileHandle *h,
    hid_t sid, char *name, int *dims, double *data)
{
  int cdims[2];
  cdims[0] = dims[1];
  cdims[1] = dims[0];
  return scorpio_create_probe(sid, name, H5T_NATIVE_DOUBLE, 2, cdims, data, h->fhandle, &h->myIOgroup);
}
