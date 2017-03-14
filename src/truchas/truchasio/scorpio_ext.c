/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

#include <assert.h>
#include "scorpio.h"

/* Wraps initialization of the iogroup and file opening into a single function */

void scorpio_open_file_ext(const char *filename, int groupSize,
    int *fhandle, iogroup_t **myIOgroup)
{
  iogroup_conf_t myIOconfig;

  myIOconfig.numIOgroups = 0;
  myIOconfig.preferredGroupSize = groupSize;
  myIOconfig.commIncoming = MPI_COMM_WORLD;

  *myIOgroup = malloc(sizeof(iogroup_t));
  scorpio_IOgroup_init(&myIOconfig, *myIOgroup);

  *fhandle = scorpio_open_file(filename, *myIOgroup, SCORPIO_FILE_CREATE);
  if (*fhandle == FAILURE) {
    fprintf(stderr, "Opening file for writing failed.\n");
    exit(-1);
  }
}

/* Wraps closing of open objects with file closing */

void scorpio_close_file_ext(int fhandle, iogroup_t *myIOgroup)
{
  iofile_t *currfile;
  hid_t *objects;
  ssize_t num_obj;
  int i;

  if (myIOgroup->localrank == 0) {
    currfile = myIOgroup->file[fhandle];
    /* Close all the datasets */
    num_obj = H5Fget_obj_count(currfile->fid,H5F_OBJ_DATASET);
    objects = (hid_t *) malloc(num_obj*sizeof(hid_t));
    if ( num_obj > 0 ) {
        H5Fget_obj_ids(currfile->fid,H5F_OBJ_DATASET,num_obj,objects);
        for(i=0;i<num_obj;i++)
            H5Dclose(objects[i]);
    }
    free(objects);

    /* Close all the groups */
    num_obj = H5Fget_obj_count(currfile->fid,H5F_OBJ_GROUP);
    objects = (hid_t *) malloc(num_obj*sizeof(hid_t));
    if ( num_obj > 0 ) {
        H5Fget_obj_ids(currfile->fid,H5F_OBJ_GROUP,num_obj,objects);
        for(i=0;i<num_obj;i++)
            H5Gclose(objects[i]);
    }
    free(objects);
  }
  scorpio_close_file(fhandle, myIOgroup);
  scorpio_IOgroup_cleanup(myIOgroup);
}

void scorpio_write_simple_attr_int(char *attr_name, int *attr_data,
    int fhandle, char *obj_name, iogroup_t *myIOgroup)
{
  scorpio_write_simple_attr(attr_name, attr_data, SCORPIO_INTEGER,
      fhandle, obj_name, myIOgroup);
}

void scorpio_write_simple_attr_double(char *attr_name, double *attr_data,
    int fhandle, char *obj_name, iogroup_t *myIOgroup)
{
  scorpio_write_simple_attr(attr_name, attr_data, SCORPIO_DOUBLE,
      fhandle, obj_name, myIOgroup);
}

void scorpio_write_simple_attr_string(char *attr_name, char *attr_data,
    int fhandle, char *obj_name, iogroup_t *myIOgroup)
{
  scorpio_write_simple_attr(attr_name, attr_data, SCORPIO_STRING,
      fhandle, obj_name, myIOgroup);
}

void scorpio_write_attr_double(char *attr_name, double *attr_data, int ndims, int *adims,
    int fhandle, char *obj_name, iogroup_t *myIOgroup)
{
  scorpio_write_attr(attr_name, attr_data, SCORPIO_DOUBLE, ndims, adims,
      fhandle, obj_name, myIOgroup);
}


void scorpio_write_dataset2_byte(
    int8_t *vector, int ndims, int *globaldims, int *localdims,
    int fhandle, char *dset_name, iogroup_t *myIOgroup)
{
  scorpio_write_dataset2(vector, SCORPIO_BYTE, ndims, globaldims, localdims,
      fhandle, dset_name, myIOgroup);
}

void scorpio_write_dataset2_int(
    int *vector, int ndims, int *globaldims, int *localdims,
    int fhandle, char *dset_name, iogroup_t *myIOgroup)
{
  scorpio_write_dataset2(vector, SCORPIO_INTEGER, ndims, globaldims, localdims,
      fhandle, dset_name, myIOgroup);
}

void scorpio_write_dataset2_double(
    double *vector, int ndims, int *globaldims, int *localdims,
    int fhandle, char *dset_name, iogroup_t *myIOgroup)
{
  scorpio_write_dataset2(vector, SCORPIO_DOUBLE, ndims, globaldims, localdims,
      fhandle, dset_name, myIOgroup);
}

/* Scorpio-style functions for creating and incrementally writing probe datasets.
 * Assumes data is presented on global rank 0.  Does a parallel write with zero
 * sized arrays on ranks != 0.
 * TODO: change api to be generic -- no hardcoded probes group
 * TODO: serial write on rank 0 only?
 * TODO: separate dataset creation from initial write?
 * TODO: change api to specify global rank data is presented on.
 */

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

void scorpio_write_probe_data_double(
    hid_t pid, int ndims, int *dims, double *data, int fhandle, iogroup_t *myIOgroup)
{
  scorpio_write_probe_data(pid, H5T_NATIVE_DOUBLE, ndims, dims, data, fhandle, myIOgroup);
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

hid_t scorpio_create_probe_double(
    hid_t sid, char *name, int ndims, int *dims, double *data, int fhandle, iogroup_t *myIOgroup)
{
  return scorpio_create_probe(sid, name, H5T_NATIVE_DOUBLE, ndims, dims, data, fhandle, myIOgroup);
}
