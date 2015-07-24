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
 * hex_convert.c
 *
 *  Convert utilities from Truchas format (degenerate hexes) into XDMF style.
 *
 */

#if HAVE_CONFIG_H
# include <danu_config.h>
#endif

#include <string.h>

#include <hdf5.h>

#include <danu_error.h>
#include <danu_h5_error.h>
#include <danu_h5_object.h>
#include <danu_types.h>
#include <danu_memory.h>
#include <danu_utils.h>
#include <danu_link.h>
#include <danu_file.h>

#include <danu_group.h>

/*
 * Takes element connectivity in Truchas format (degenerated hexes), converts
 * to XDMF 1D mixed array and saves it into a new HDF5 file.
 *
 * Arguments:
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
int danu_hex_save(const char *name, int* elements, int n1, int n2, const char
    *group_name, const char *dataset_name)
{
    int *elements1d;
    int n, i, nc, idx, list[8];
    /* Convert from Truchas vertex ordering to XDMF vertex ordering (which is
     * the same as the VTK ordering). The XDMF cell type numbers are defined in
     * XdmfTopologyType.cpp in git://public.kitware.com/Xdmf2.git. The vertices
     * are enumerated from 0, so we need to subtract 1 from our node ids. */
    n = n1 * 9;
    elements1d = malloc(n*sizeof(int));
    idx = 0;
    for (nc=0; nc < n1; nc++) {
        for (i=0; i < n2; i++) {
            list[i] = elements[nc*n2+i];
        }
        if (list[0] == list[1]) { /* tet element */
            elements1d[idx] = 6; /* Tet XDMF cell type */
            idx += 1;
            for (i=0; i < 4; i++) elements1d[idx+i] = list[1+i]-1;
            idx += 4;
        } else if (list[4] == list[5]) { /* pyramid element */
            elements1d[idx] = 7; /* Pyramid XDMF cell type */
            idx += 1;
            for (i=0; i < 5; i++) elements1d[idx+i] = list[i]-1;
            idx += 5;
        } else if (list[5] == list[6]) { /* wedge element */
            /* Convert from Truchas ordering to Exodus ordering */
            i = list[1]; list[1] = list[3]; list[3] = i;  /* swap 1 and 3 */
            i = list[2]; list[2] = list[4]; list[4] = i;  /* swap 2 and 4 */
            /* Convert from Exodus ordering to VTK / XDMF ordering */
            i = list[1]; list[1] = list[2]; list[2] = i;  /* swap 1 and 2 */
            i = list[4]; list[4] = list[5]; list[5] = i;  /* swap 4 and 5 */
            elements1d[idx] = 8; /* Wedge XDMF cell type */
            idx += 1;
            for (i=0; i < 6; i++) elements1d[idx+i] = list[i]-1;
            idx += 6;
        } else { /* hex element */
            elements1d[idx] = 9; /* Hex XDMF cell type */
            idx += 1;
            for (i=0; i < 8; i++) elements1d[idx+i] = list[i]-1;
            idx += 8;
        }
    }
    n = idx;


    hid_t file, space, dset, group;
    herr_t status;
    hsize_t dims[1] = {n};

    file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    group = H5Gcreate(file, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    space = H5Screate_simple(1, dims, NULL);
    dset = H5Dcreate(group, dataset_name, H5T_STD_I32LE, space, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                elements1d);
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Fclose(file);

    free(elements1d);

    return n;
}
