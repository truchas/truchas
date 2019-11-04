/*
 *  Extensions to Chaparral for extracting the view factor matrix
 *
 *  Neil N. Carlson <nnc@lanl.gov>
 *  Los Alamos National Laboratory.
 *
 *  VF_GetMatrix -- Extracts the view factor matrix in CSR format
 *
 *    encl            enclosure number
 *    vf_cnt          number of nonzero view factors in each row
 *    vf_index        column indices of the nonzero view factors
 *    vf_data         nonzero view factors stored by row
 *    write_vf_diag   if non-zero, the diagonal elements of the matrix are
 *                    written to the vf_diag array.  Otherwise, no data is written.
 *    vf_diag         if write_vf_diag is non-zero, the diagonal elements of the
 *                    matrix are returned in this array, including zero values, and
 *                    excluded from the vf_index/vf_data arrays.
 *    write_vf_virt   if non-zero, the view factors in the virtual patch column are
 *                    written to the vf_virt array.  Othersie, no data is written.
 *                    vf_diag array. Otherwise, no data will be written.
 *    vf_virt         if write_vf_virt is non-zero, the view factors in the virtual
 *                    patch column are are returned in this array, including
 *                    implicit zero values, and excluded from the vf_index/vf_data
 *                    arrays.  In addition, no data is returned for the virtual
 *                    patch row.
 *
 *  The row for each patch is returned on the process the patch was passed to
 *  VF_DefineEnclosure, and in the same order as the global_ids array.  The
 *  column indices returned in vf_index refer to the global ids.  To have the
 *  indices refer to the caller's global index for a patch instead, undefine
 *  the VF_MAP_TO_HOST_GID macro below. The values in each row are sorted with
 *  increasing vf_index value.
 *
 *  The virtual surface patch must be assigned the largest global id.  In
 *  addition it seems to be conventional for it to be the last patch on the
 *  last process.  If this is the case, and a non-zero value for write_vf_virt
 *  is passed, the final entry of vf_cnt, vf_diag, and vf_virt (the one
 *  corresponding to the virtual surface) is never referenced and thus those
 *  arrays can be allocated a length 1 less than normal.  So the matrix
 *  returned really looks like the matrix of the subsystem obtained by
 *  dropping the virtual surface patch as an unknown.
 *
 *  VF_GetRowCounts -- Returns the number of nonzero view factors in each row
 *
 *    encl       enclosure number
 *    mode       the counting mode
 *    count      number of nonzero view factors in each row according to mode
 *
 *  When counting the nonzero view factors in a row, the diagonal element
 *  will be ignored if mode matches the mask VF_MATRIX_EXCL_DIAG.  An
 *  element in the virtual surface patch column will be ignored if mode
 *  matches the mask VF_MATRIX_EXCL_VIRT.   Both can be ignored if mode is
 *  set to the or of the two masks.  A zero value of mode means count all
 *  the nonzeros.  When ignoring the virtual view factors, the count for
 *  the virtual row will always be zero unless it is the last patch, which
 *  is conventional, in which case no value will be assigned and count can
 *  be dimensioned to the size of actual patches.
 *
 *  Typically VF_RowCounts will be called to determine the required length
 *  of the vf_index and vf_data arrays, before passing them to VF_GetMatrix.
 *  Consequently it is essential that the mode value corresponds properly to
 *  the presence of the vf_diag and vf_virt arguments to VF_GetMatrix.
 */

#include <stdlib.h>
#include <stdio.h>

#include "vf.h"

#define VF_MAP_TO_HOST_GID  /* otherwise map to host global indices */

void
VF_GetMatrix(int encl, int vf_cnt[], int vf_index[], float vf_data[],
             int write_vf_diag, float vf_diag[], int write_vf_virt, float vf_virt[])
{
  int i, j, k, n, proc, size, mode, cnt, offset, virt, self, err;
  int *orig_proc, *dest_proc, *index, *row_index, *cnt_g, *index_map, host_npatches, npatches_l;
  float *diag_g, *virt_g, *data, *row_data;
  VFrow *row;
  VFenclosure *e=VF_GetEnclosure(encl);
#ifndef VF_NO_MPI
  MPI_Request *request1=NULL, *request2=NULL;
  MPI_Status status;

  if (VFLIB_Size > 1) {

    mode = 0;
    if (write_vf_diag) {
      mode = mode|VF_MATRIX_EXCL_DIAG;
      diag_g = VF_Newf(e->npatches_g);
    } else {
      diag_g = NULL;
      self = -1;
    }
    if (e->partial && write_vf_virt) {
      mode = mode|VF_MATRIX_EXCL_VIRT;
      virt_g = VF_Newf(e->npatches_g);
      virt = e->npatches_g - 1; /* the global index of the virtual patch */
    } else {
      virt_g = NULL;
      virt = -1;
    }

    /* ADJUST THE NUMBER OF LOCAL PATCHES TO EXCLUDE THE VIRTUAL PATCH IF  */
    /* REQUESTED AND POSSIBLE (= IT'S LAST PATCH ON THE LAST HOST PROCESS) */
    host_npatches  = e->host_npatches;
    npatches_l = e->npatches_l;
    if ((mode&VF_MATRIX_EXCL_VIRT)&&(VFLIB_Rank==VFLIB_Size-1)) {
      if (e->host2vflib_map[host_npatches-1] == e->npatches_g-1) {
        host_npatches -= 1;
        npatches_l -= 1;
      }
    }

    /* ROW COUNT FOR EACH GLOBAL INDEX */
    cnt_g = VF_Newi(e->npatches_g);
    VF_GetRowCounts_Aux(mode,cnt_g);

    /* ORIGINATING PROCESS FOR EACH GLOBAL INDEX */
    orig_proc = VF_Newi(e->npatches_g);
    for (proc=0; proc<VFLIB_Size; proc++) {
      for (j=0; j<e->comm_plan.partners[proc].cnt; j++) {
        orig_proc[e->comm_plan.partners[proc].global_index[j]] = proc;
      }
    }

    /* POST ALL THE RECEIVES */
    request1 = (MPI_Request*)VF_Newv(host_npatches*sizeof(MPI_Request));
    request2 = (MPI_Request*)VF_Newv(host_npatches*sizeof(MPI_Request));
    for (j=0, offset=0; j<host_npatches; j++) {
      n = e->host2vflib_map[j]; /* use the global index as the tag */
      cnt = cnt_g[n];
      vf_cnt[j] = cnt;
      if (cnt<=0) continue;
      err = MPI_Irecv(&vf_index[offset],cnt,MPI_INT,   orig_proc[n],n,VFLIB_Comm,&request1[j]);
      if (err) {
        fprintf(stderr,"MPI error in VF_GetMatrix(): MPI_Irecv() in %s, line %d\n",__FILE__,__LINE__);
        VF_Exit(1);
      }
      err = MPI_Irecv(&vf_data[offset], cnt,MPI_FLOAT, orig_proc[n],n,VFLIB_Comm,&request2[j]);
      if (err) {
        fprintf(stderr,"MPI error in VF_GetMatrix(): MPI_Irecv() in %s, line %d\n",__FILE__,__LINE__);
        VF_Exit(1);
      }
      offset += cnt;
    }

#ifdef VF_MAP_TO_HOST_GID
    /* DESTINATION PROCESS FOR EACH GLOBAL INDEX */
    dest_proc = orig_proc;
    for (proc=0; proc<VFLIB_Size; proc++) {
      if (proc==VFLIB_Rank) {
        size = e->host_npatches;
        index = e->host2vflib_map;
      } else {
        VF_GetINTbuffer_ptr(&index);
      }
      VF_BroadcastInt(&size,1,proc);
      VF_BroadcastInt(index,size,proc);
      for (j=0; j<size; j++) {
        dest_proc[index[j]] = proc;
      }
    }

    /* MAPPING OF GLOBAL INDICES TO HOST GLOBAL IDS */
    index_map = VF_Newi(e->npatches_g);
    row = e->row;
    for (j=0; j<e->npatches_l; j++, row++) {
      index_map[row->global_index] = row->host_gid;
    }
    VF_ExchangeInt(index_map);
#else
    /* DESTINATION PROCESS FOR EACH GLOBAL INDEX AND    */
    /* MAPPING OF GLOBAL INDICES TO HOST GLOBAL INDICES */
    dest_proc = orig_proc;
    index_map = VF_Newi(e->npatches_g);
    for (proc=0, offset=0; proc<VFLIB_Size; proc++) {
      if (proc==VFLIB_Rank) {
        size = e->host_npatches;
        index = e->host2vflib_map;
      } else {
        VF_GetINTbuffer_ptr(&index);
      }
      VF_BroadcastInt(&size,1,proc);
      VF_BroadcastInt(index,size,proc);
      for (j=0; j<size; j++) {
        dest_proc[index[j]] = proc;
        index_map[index[j]] = offset+j;
      }
      offset += size;
    }
#endif

    /* POST ALL THE SENDS */
    VF_GetINTbuffer_ptr(&row_index);
    VF_GetSPbuffer0_ptr(&row_data);
    row = e->row;
    for (j=0; j<npatches_l; j++, row++) {
      n = row->global_index;
      cnt = cnt_g[n];
      if (virt_g!=NULL) virt_g[n] = 0.0;
      if (diag_g!=NULL) {
        self = n;
        diag_g[n] = 0.0;
      }
      /* ASSEMBLE THE SPARSE ROW DATA */
      k = 0;
      index = row->array0.index;
      data  = row->array0.data;
      for (i=0; i<row->array0.cnt; i++) {
        if (index[i]==self) {
          diag_g[n] += data[i];
        } else if (index[i]==virt) {
          virt_g[n] += data[i];
        } else {
          row_index[k] = index_map[index[i]];
          row_data[k]  = data[i];
          k++;
        }
      }
      index = row->array1.index;
      data  = row->array1.data;
      for (i=0; i<row->array1.cnt; i++) {
        if (index[i]==virt) {
          virt_g[n] += data[i];
        } else {
          row_index[k] = index_map[index[i]];
          row_data[k]  = data[i];
          k++;
        }
      }
      if (cnt<=0) continue;
      VF_SortSparseArrayAux(row_index,row_data,cnt);
      /* SEND THE ROW DATA */
      err = MPI_Send(row_index,cnt,MPI_INT,   dest_proc[n],n,VFLIB_Comm);
      if (err) {
        fprintf(stderr,"MPI error in VF_GetMatrix(): MPI_Send() in %s, line %d\n",__FILE__,__LINE__);
        VF_Exit(1);
      }
      err = MPI_Send(row_data, cnt,MPI_FLOAT, dest_proc[n],n,VFLIB_Comm);
      if (err) {
        fprintf(stderr,"MPI error in VF_GetMatrix(): MPI_Send() in %s, line %d\n",__FILE__,__LINE__);
        VF_Exit(1);
      }
    }
    VF_Free(index_map);
    VF_Free(dest_proc);

    /* ASSEMBLE THE MATRIX DIAGONAL IF REQUESTED */
    if (write_vf_diag) {
      VF_ExchangeFloat(diag_g);
      for (j=0; j<host_npatches; j++) {
        vf_diag[j] = diag_g[e->host2vflib_map[j]];
      }
      VF_Free(diag_g);
    }

    /* ASSEMBLE THE VIRTUAL MATRIX COLUMN IF REQUESTED */
    if (write_vf_virt) {
      if (e->partial) {
        VF_ExchangeFloat(virt_g);
        for (j=0; j<host_npatches; j++) {
          vf_virt[j] = virt_g[e->host2vflib_map[j]];
        }
        VF_Free(virt_g);
      } else {
        for (j=0; j<e->host_npatches; j++) {
          vf_virt[j] = 0.0;
        }
      }
    }

    /* WAIT FOR THE RECEIVES TO COMPLETE */
    for (j=0; j<host_npatches; j++) {
      n = e->host2vflib_map[j];
      if (cnt_g[n]<=0) continue;
      if (MPI_Wait(&request1[j],&status)) {
        fprintf(stderr,"MPI error in VF_GetMatrix(): MPI_Wait() in %s, line %d\n",__FILE__,__LINE__);
        VF_Exit(1);
      }
      if (MPI_Wait(&request2[j],&status)) {
        fprintf(stderr,"MPI error in VF_GetMatrix(): MPI_Wait() in %s, line %d\n",__FILE__,__LINE__);
        VF_Exit(1);
      }
    }
    VF_Free(request1);
    VF_Free(request2);
    VF_Free(cnt_g);

  } else {
#endif
    host_npatches = e->host_npatches;
    if (!write_vf_diag) self = -1;
    if (e->partial && write_vf_virt) {
      virt = e->npatches_g - 1; /* the global index of the virtual patch */
      /* adjust the number of patches to exclude the virtual patch if possible */
      if (e->host2vflib_map[host_npatches-1]==virt) host_npatches--;
    } else {
      virt = -1;
    }

#ifdef VF_MAP_TO_HOST_GID
    /* MAPPING OF GLOBAL INDICES TO HOST GLOBAL IDS */
    VF_GetINTbuffer_ptr(&index_map);
    row = e->row;
    for (j=0; j<e->npatches_g; j++, row++) {
      index_map[row->global_index] = row->host_gid;
    }
#else
    /* MAPPING OF GLOBAL INDICES TO HOST GLOBAL INDICES */
    VF_GetINTbuffer_ptr(&index_map);
    for (j=0; j<e->npatches_g; j++) {
      index_map[e->host2vflib_map[j]] = j;
    }
#endif

    for (j=0, size=0; j<host_npatches; j++) {
      n = e->host2vflib_map[j];
      row = &(e->row[n]);  /* assume rows are indexed by global index */
      if (write_vf_virt) vf_virt[j] = 0.0;
      if (write_vf_diag) {
        self = n;
        vf_diag[j] = 0.0;
      }
      /* ASSEMBLE AND STORE THE SPARSE ROW DATA */
      offset = size;
      index = row->array0.index;
      data  = row->array0.data;
      for (i=0; i<row->array0.cnt; i++) {
        if (index[i]==self) {
          vf_diag[j] += data[i];
        } else if (index[i]==virt) {
          vf_virt[j] += data[i];
        } else {
          vf_index[size] = index_map[index[i]];
          vf_data[size]  = data[i];
          size++;
        }
      }
      index = row->array1.index;
      data  = row->array1.data;
      for (i=0; i<row->array1.cnt; i++) {
        if (index[i]==virt) {
          vf_virt[j] += data[i];
        } else {
          vf_index[size] = index_map[index[i]];
          vf_data[size]  = data[i];
          size++;
        }
      }
      vf_cnt[j] = size - offset;
      if (vf_cnt[j]>0) VF_SortSparseArrayAux(&vf_index[offset],&vf_data[offset],vf_cnt[j]);
    }
#ifndef VF_NO_MPI
  }
#endif
}

void
VF_GetRowCounts(int encl, int mode, int *count)
{
  int i, *itmp, host_npatches;
  VFenclosure *e=VF_GetEnclosure(encl);

  if (!e->partial) mode = mode&~VF_MATRIX_EXCL_VIRT;

  host_npatches = e->host_npatches;
  if ((mode&VF_MATRIX_EXCL_VIRT)&&(VFLIB_Rank==VFLIB_Size-1)) {
    if (e->host2vflib_map[host_npatches-1] == e->npatches_g-1) host_npatches--;
  }

  VF_GetINTbuffer_ptr(&itmp);
  VF_GetRowCounts_Aux(mode, itmp);
  for (i=0; i<host_npatches; i++) {
    count[i] = itmp[e->host2vflib_map[i]];
  }
}

void
VF_GetRowCounts_Aux(int mode, int *count)
{
  int i, j, n_l, cnt, virt, *index;
  VFrow *row;
  VFenclosure *e=VF_CurrentEnclosure();

  n_l = e->npatches_l;
  if (mode&VF_MATRIX_EXCL_VIRT) {
    virt = e->npatches_g - 1;
    if (VFLIB_Rank == VFLIB_Size-1) {
      count[virt] = 0;
      n_l -= 1;
    }
  }

  row = e->row;
  for (j=0; j<n_l; j++, row++) {
    cnt = row->array0.cnt + row->array1.cnt;
    /* decrement count for reference to virtual patch element */
    if (mode&VF_MATRIX_EXCL_VIRT) {
      index = row->array0.index;
      for (i=0; i<row->array0.cnt; i++) {
        if (index[i] == virt) cnt--;
      }
      index = row->array1.index;
      for (i=0; i<row->array1.cnt; i++) {
        if (index[i] == virt) cnt--;
      }
    }
    /* decrement count for reference to diagonal element */
    if (mode&VF_MATRIX_EXCL_DIAG) {
      index = row->array0.index;
      for (i=0; i<row->array0.cnt; i++) {
        if (index[i] == row->global_index) cnt--;
      }
    }
    count[row->global_index] = cnt;
  }

  VF_ExchangeInt(count);
}
