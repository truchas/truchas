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

#include <string.h>

#include <danu_error.h>
#include <danu_h5_error.h>
#include <danu_h5_object.h>
#include <danu_types.h>
#include <danu_memory.h>
#include <danu_utils.h>
#include <danu_link.h>
#include <danu_file.h>

#include <danu_group.h>

#include <danu_mesh_convert.h>

struct DanuH5Handle * danu_h5_create_handle()
{
    return malloc(sizeof(struct DanuH5Handle));
}
void danu_h5_free_handle(struct DanuH5Handle *h)
{
    free(h);
}

void danu_h5_open(struct DanuH5Handle *h, const char *name,
        const char *group_name)
{
    h->file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    h->group = H5Gcreate(h->file, group_name, H5P_DEFAULT, H5P_DEFAULT,
            H5P_DEFAULT);
}

void danu_h5_save_int_1d_array(struct DanuH5Handle *h,
        const char *dataset_name, int* array1d, int n)
{
    hid_t space, dset;
    herr_t status;
    hsize_t dims[1] = {n};

    space = H5Screate_simple(1, dims, NULL);
    dset = H5Dcreate(h->group, dataset_name, H5T_STD_I32LE, space, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                array1d);
    status = H5Dclose(dset);
    status = H5Sclose(space);
}

void danu_h5_close(struct DanuH5Handle *h)
{
    herr_t status;
    status = H5Gclose(h->group);
    status = H5Fclose(h->file);
}

int danu_hex_save(struct DanuH5Handle *h, int* elements, int n1, int n2,
        const char *dataset_name)
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

    danu_h5_save_int_1d_array(h, dataset_name, elements1d, n);

    free(elements1d);

    return n;
}
