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
 * danu_mesh_convert.h
 *
 *  Convert utilities from Truchas format (degenerate hexes) into XDMF style.
 *
 */

#ifndef DANU_MESH_CONVERT_H
#define DANU_MESH_CONVERT_H

#if HAVE_CONFIG_H
# include <danu_config.h>
#endif

#include <hdf5.h>

struct DanuH5Handle {
    hid_t file, group;
};

struct DanuH5Handle * danu_h5_create_handle();
void danu_h5_free_handle(struct DanuH5Handle *h);
void danu_h5_open(struct DanuH5Handle *h, const char *name,
        const char *group_name);
void danu_h5_save_int_1d_array(struct DanuH5Handle *h,
        const char *dataset_name, int* array1d, int n);
void danu_h5_close(struct DanuH5Handle *h);

/*
 * Takes element connectivity in Truchas format (degenerated hexes), converts
 * to XDMF 1D mixed array and saves it into a new HDF5 file.
 *
 * Arguments:
 *   h .............. The DanuH5Handle handle
 *   name ........... filename of the HDF5 file to save into
 *   elements ....... 2D array of dimensions (n1, n2) of the Truchas element
 *       connectivity. n1 should be the number of elements and n2 should be 8.
 *   n1, n2 ......... dimensions of the 'elements' array
 *   group_name ..... the name of the HDF5 group
 *   dataset_name ... the name of the HDF5 dataset. The array will be saved
 *       into group_name/dataset_name. Note: currently nested groups (like
 *       group_name1/group_name2/group_name_3/dataset_name) are not
 *       implemented, so there can only be one group at the moment.
 *
 * The function returns the length of the 1D array (that was saved into HDF5).
 *
 */
int danu_hex_save(struct DanuH5Handle *h, int* elements, int n1, int n2,
        const char *dataset_name);

#endif
