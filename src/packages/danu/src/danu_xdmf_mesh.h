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
 * danu_xdmf_mesh.h
 *
 *  Convert utilities from Truchas format (degenerate hexes) into XDMF style.
 *  We create a new HDF5 file with the arrays that XDMF can understand. We have
 *  to convert the connectivity array as well as element blocks (materials).
 *
 */


#ifndef DANU_MESH_CONVERT_H
#define DANU_MESH_CONVERT_H

#if HAVE_CONFIG_H
# include <danu_config.h>
#endif

#include <hdf5.h>

/*
 * This structure is initialized by 'danu_h5_open'. You have to call
 * 'danu_h5_close' when you are done with it. The structure can be allocated on
 * stack, or you can allocate it on heap using 'danu_h5_free_handle' and
 * 'danu_h5_free_handle' (you still have to call 'danu_h5_open'/'close' after
 * that).
 */
struct DanuXDMFMeshHandle {
    hid_t file, group;
};

/*
 * Returns a pointer to a new heap allocated DanuXDMFMeshHandle struct.
 * This is used in the Python wrappers.
 */
struct DanuXDMFMeshHandle * danu_h5_create_handle();

/*
 * Frees the heap allocated DanuXDMFMeshHandle struct.
 */
void danu_h5_free_handle(struct DanuXDMFMeshHandle *h);

/*
 * Opens the HDF5 mesh file
 *
 * Arguments:
 *   name ........... Filename of the HDF5 file to save into.
 *   group_name ..... The name of the HDF5 group. All datasets will be saved
 *                    under the this group. Note: currently nested groups (like
 *                    group_name1/group_name2/group_name_3/dataset_name) are
 *                    not implemented, so there can only be one group at the
 *                    moment.
 *
 */
void danu_h5_open(struct DanuXDMFMeshHandle *h, const char *name,
        const char *group_name);

/*
 * Saves a 1D array into the HDF5 mesh file (you have to call danu_h5_open
 * first).
 *
 * Arguments:
 *   dataset_name ... the name of the HDF5 dataset where array will be saved.
 */
void danu_h5_save_int_1d_array(struct DanuXDMFMeshHandle *h,
        const char *dataset_name, int* array1d, int n);

/*
 * Closes the DanuXDMFMeshHandle handle 'h'.
 */
void danu_h5_close(struct DanuXDMFMeshHandle *h);

/*
 * Takes element connectivity in Truchas format (degenerated hexes), converts
 * to XDMF 1D mixed array and saves it into a new HDF5 file.
 *
 * Arguments:
 *   h .............. The DanuXDMFMeshHandle handle
 *   elements ....... 2D array of dimensions (n1, n2) of the Truchas element
 *       connectivity. n1 should be the number of elements and n2 should be 8.
 *   n1, n2 ......... dimensions of the 'elements' array
 *   dataset_name ... the name of the HDF5 dataset where array will be saved
 *
 * The function returns the length of the 1D array (that was saved into HDF5).
 *
 */
int danu_hex_save(struct DanuXDMFMeshHandle *h, int* elements, int n1, int n2,
        const char *dataset_name);

#endif
