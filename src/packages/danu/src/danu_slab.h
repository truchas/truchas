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
* danu_slab.h
*
*  DANU hyperslab module
*
*/

#ifndef DANU_SLAB_H
#define DANU_SLAB_H

#include <hdf5.h>

#include <danu_types.h>

/* Hyperslab structure */
typedef struct dslab_struct {
    int dim;          /* Hyperslab dimension MUST match dataspace dimension */

    hsize_t *offset;  /* Offset */
    hsize_t *count;   /* Number of blocks */
    hsize_t *stride;  /* Block stride */
    hsize_t *block;   /* Block dimension */
} dslab_t;

dslab_t* danu_slab_alloc(hid_t dataspace);
void     danu_slab_free(dslab_t *slab);

herr_t danu_slab_select(hid_t dataspace, dslab_t *slab);
herr_t danu_slab_contiguous(dslab_t *slab, hsize_t *offset, hsize_t *size);

#endif




