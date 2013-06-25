/****************************************************************************\
 *                                                                          *
 *   Copyright (c) 1995, 2000, 2005 Sandia Corporation.                           *
 *   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, *
 *   the U.S. Government retains certain rights in this software.           *
 *   For more info, see the README file in the top-level directory.         * 
 *                                                                          *
\****************************************************************************/

/*
@(#)++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
@(#)
@(#)    $RCSfile: VF_ParallelUtils.c,v $
@(#)    $Revision: 1.3 $  $Date: 2006/03/17 04:59:11 $  $Author: mwglass $
@(#)    $Source: /cvsroot/chaparral/chaparral/vflib/VF_ParallelUtils.c,v $
@(#)
@(#)    DESCRIPTION:  Misc. parallel utility functions.
@(#)
@(#)--------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "vf.h"

void VF_Exit(int exit_code)
{
#ifndef VF_NO_MPI
  MPI_Abort(VFLIB_Comm,exit_code);
#endif
}

void VF_Sync(void)
{
#ifndef VF_NO_MPI
  int err;
  
  err = MPI_Barrier(VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Sync()\n");
    VF_Exit(1);
  }
#endif
}

void VF_PrintSyncStart(void)
{
#ifndef VF_NO_MPI
  int flag=1, tag=10000, src, err;
  MPI_Status status;

  if (VFLIB_Rank != 0) {
    src = VFLIB_Rank-1;
    err = MPI_Recv(&flag, 1, MPI_INT, src, 
                   tag, VFLIB_Comm, &status);
    if (err != 0) {
      fprintf(stderr, "VF_PrintSyncStart: ERROR on processor %d\n", VFLIB_Rank);
      VF_Exit(-1);
    }
  }
#endif
}

void VF_PrintSyncEnd(void)
{
#ifndef VF_NO_MPI
  int flag=1, tag=10000, src, dest, err;
  MPI_Status status;

  fflush(stdout); sleep(2);
  if (VFLIB_Rank<VFLIB_Size-1) {
    dest = VFLIB_Rank+1;
  } else {
    dest = 0;
  }
  err = MPI_Send(&flag, 1, MPI_INT, dest, tag, VFLIB_Comm);
  if (err != 0) {
    fprintf(stderr, "VF_PrintSyncEnd: ERROR on processor %d\n", VFLIB_Rank);
    VF_Exit(-1);
  }
  if (VFLIB_Rank==0) {
    src = VFLIB_Size-1;
    err = MPI_Recv(&flag, 1, MPI_INT, src, 
                   tag, VFLIB_Comm, &status);
    if (err != 0) {
      fprintf(stderr, "PrintSyncStop: ERROR on processor %d\n", VFLIB_Rank);
      VF_Exit(-1);
    }
  }
  VF_Sync();
#endif
}

void VF_StagedIO_Start(void)
{
#ifndef VF_NO_MPI
  int flag=1, tag=10000, src, err;
  MPI_Status status;

  if (!VFLIB_StagedIO) return;
  if (VFLIB_Rank>=VFLIB_Ncntrls) {
    src = VFLIB_Rank-VFLIB_Ncntrls;
    err = MPI_Recv(&flag, 1, MPI_INT, src, 
                   tag+src, VFLIB_Comm, &status);
    if (err != 0) {
      fprintf(stderr, "StagedIO_Start: ERROR on processor %d\n", VFLIB_Rank);
      VF_Exit(-1);
    }
  }
#endif
}

void VF_StagedIO_End(void)
{
#ifndef VF_NO_MPI
  int flag=1, tag=10000, dest, err;

  if (!VFLIB_StagedIO) return;
  if (VFLIB_Rank+VFLIB_Ncntrls < VFLIB_Size) {
    dest = VFLIB_Rank+VFLIB_Ncntrls;
    err  = MPI_Send(&flag, 1, MPI_INT, dest, tag+VFLIB_Rank, VFLIB_Comm);
    if (err != 0) {
      fprintf(stderr, "StagedIO_End: ERROR on processor %d\n", VFLIB_Rank);
      VF_Exit(-1);
    }
  }
#endif
}

void VF_GlobalSum_Int(int *value, char *file, int line)
{
#ifndef VF_NO_MPI
  int err, tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Allreduce(&tmp,value,1,MPI_INT,MPI_SUM,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_GlobalSum_Int() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_GlobalSum_Float(float *value, char *file, int line)
{
#ifndef VF_NO_MPI
  int   err;
  float tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Allreduce(&tmp,value,1,MPI_FLOAT,MPI_SUM,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_GlobalSum_Float() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_GlobalSum_Double(double *value, char *file, int line)
{
#ifndef VF_NO_MPI
  int    err;
  double tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Allreduce(&tmp,value,1,MPI_DOUBLE,MPI_SUM,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_GlobalSum_Double() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_GlobalMin_Int(int *value, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  int tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Allreduce(&tmp,value,1,MPI_INT,MPI_MIN,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_GlobalMin_Int() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_GlobalMin_Float(float *value, char *file, int line)
{
#ifndef VF_NO_MPI
  int   err;
  float tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Allreduce(&tmp,value,1,MPI_FLOAT,MPI_MIN,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_GlobalMin_Float() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_GlobalMin_Double(double *value, char *file, int line)
{
#ifndef VF_NO_MPI
  int    err;
  double tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Allreduce(&tmp,value,1,MPI_DOUBLE,MPI_MIN,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_GlobalMin_Double() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_GlobalMax_Int(int *value, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  int tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Allreduce(&tmp,value,1,MPI_INT,MPI_MAX,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_GlobalMax_Int() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_GlobalMax_Float(float *value, char *file, int line)
{
#ifndef VF_NO_MPI
  int   err;
  float tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Allreduce(&tmp,value,1,MPI_FLOAT,MPI_MAX,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_GlobalMax_Float() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_GlobalMax_Double(double *value, char *file, int line)
{
#ifndef VF_NO_MPI
  int    err;
  double tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Allreduce(&tmp,value,1,MPI_DOUBLE,MPI_MAX,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_GlobalMax_DOUBLE() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Sum_Int(int *value, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int err, tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Reduce(&tmp,value,1,MPI_INT,MPI_SUM,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Sum_Int() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Sum_Float(float *value, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int   err;
  float tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Reduce(&tmp,value,1,MPI_FLOAT,MPI_SUM,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Sum_Float() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Sum_Double(double *value, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int    err;
  double tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Reduce(&tmp,value,1,MPI_DOUBLE,MPI_SUM,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Sum_Double() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Min_Int(int *value, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  int tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Reduce(&tmp,value,1,MPI_INT,MPI_MIN,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Min_Int() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Min_Float(float *value, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int   err;
  float tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Reduce(&tmp,value,1,MPI_FLOAT,MPI_MIN,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Min_Float() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Min_Double(double *value, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int    err;
  double tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Reduce(&tmp,value,1,MPI_DOUBLE,MPI_MIN,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Min_Double() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Max_Int(int *value, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  int tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Reduce(&tmp,value,1,MPI_INT,MPI_MAX,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Max_Int() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Max_Float(float *value, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int   err;
  float tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Reduce(&tmp,value,1,MPI_FLOAT,MPI_MAX,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Max_Float() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Max_Double(double *value, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int    err;
  double tmp;
  
  if (VFLIB_Size==1) return;
  tmp = *value;
  err = MPI_Reduce(&tmp,value,1,MPI_DOUBLE,MPI_MAX,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Max_Double() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Broadcast_Char(char *values, int size, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  
  if (VFLIB_Size==1) return;
  err = MPI_Bcast(values,size,MPI_BYTE,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Broadcast_Char() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Broadcast_Int(int *values, int size, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  
  if (VFLIB_Size==1) return;
  err = MPI_Bcast(values,size,MPI_INT,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Broadcast_Int() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Broadcast_Float(float *values, int size, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  
  if (VFLIB_Size==1) return;
  err = MPI_Bcast(values,size,MPI_FLOAT,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Broadcast_Float() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Broadcast_Double(double *values, int size, int root, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  
  if (VFLIB_Size==1) return;
  err = MPI_Bcast(values,size,MPI_DOUBLE,root,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Broadcast_Double() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Allgather_Int(int *sbuffer, int *rbuffer, int size, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  
  err = MPI_Allgather(sbuffer,size,MPI_INT,rbuffer,size,MPI_INT,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Allgather_Int() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Allgather_Float(float *sbuffer, float *rbuffer, int size, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  
  err = MPI_Allgather(sbuffer,size,MPI_FLOAT,
                      rbuffer,size,MPI_FLOAT,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Allgather_Float() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Allgather_Double(double *sbuffer, double *rbuffer, int size, char *file, int line)
{
#ifndef VF_NO_MPI
  int err;
  
  err = MPI_Allgather(sbuffer,size,MPI_DOUBLE,
                      rbuffer,size,MPI_DOUBLE,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Allgather_Double() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Allgatherv_Int(int *snd_buffer, int *rcv_buffer, 
                       int snd_size, int *rcv_size, int *offsets, 
                       char *file, int line)
{
#ifndef VF_NO_MPI
  int err;

  err = MPI_Allgatherv(snd_buffer,snd_size,MPI_INT,
                       rcv_buffer,rcv_size,offsets,
                       MPI_INT,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Allgatherv_Int() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void VF_Allgatherv_Double(double *snd_buffer, double *rcv_buffer, 
                          int snd_size, int *rcv_size, int *offsets, 
                          char *file, int line)
{
#ifndef VF_NO_MPI
  int err;

  err = MPI_Allgatherv(snd_buffer,snd_size,MPI_DOUBLE,
                       rcv_buffer,rcv_size,offsets,
                       MPI_DOUBLE,VFLIB_Comm);
  if (err) {
    fprintf(stderr,"MPI error in VF_Allgatherv_Double() in %s, line %d\n",file,line);
    VF_Exit(1);
  }
#endif
}

void
VF_Exchange_Int(int x[], char *file, int line)
{
#ifndef VF_NO_MPI
  int    i, *index, proc, err;
  float *buffer;
  int *xtmp;
  VFenclosure *e=VF_CurrentEnclosure();

  if (VFLIB_Size==1) return;
  VF_GetSPbuffer2_ptr(&buffer);
  xtmp = (int*)buffer;
  for (proc=0; proc<VFLIB_Size; proc++) {
    if (e->comm_plan.partners[proc].cnt>0) {
      if (proc==VFLIB_Rank) {
        index = e->comm_plan.partners[proc].global_index;
        for (i=0; i<e->comm_plan.partners[proc].cnt; i++) {
          xtmp[i] = x[*index++];
        }
      }
      err = MPI_Bcast(xtmp,e->comm_plan.partners[proc].cnt,
                      MPI_INT,proc,VFLIB_Comm);
      if (err) {
        fprintf(stderr,"MPI error in VF_Exchange_Int() in %s, line %d\n",file,line);
        VF_Exit(1);
      }
      if (proc!=VFLIB_Rank) {
        index = e->comm_plan.partners[proc].global_index;
        for (i=0; i<e->comm_plan.partners[proc].cnt; i++) {
          x[*index++] = xtmp[i];
        }
      }
    }
  }
#endif
}

void VF_Exchange_Float(float x[], char *file, int line)
{
#ifndef VF_NO_MPI
  int    i, *index, proc, err;
  float *xtmp;
  VFenclosure *e=VF_CurrentEnclosure();

  if (VFLIB_Size==1) return;
  VF_GetSPbuffer2_ptr(&xtmp);
  for (proc=0; proc<VFLIB_Size; proc++) {
    if (e->comm_plan.partners[proc].cnt>0) {
      if (proc==VFLIB_Rank) {
        index = e->comm_plan.partners[proc].global_index;
        for (i=0; i<e->comm_plan.partners[proc].cnt; i++) {
          xtmp[i] = x[*index++];
        }
      }
      err = MPI_Bcast(xtmp,e->comm_plan.partners[proc].cnt,
                      MPI_FLOAT,proc,VFLIB_Comm);
      if (err) {
        fprintf(stderr,"MPI error in VF_Exchange_Float() in %s, line %d\n",file,line);
        VF_Exit(1);
      }
      if (proc!=VFLIB_Rank) {
        index = e->comm_plan.partners[proc].global_index;
        for (i=0; i<e->comm_plan.partners[proc].cnt; i++) {
          x[*index++] = xtmp[i];
        }
      }
    }
  }
#endif
}

void VF_Exchange_Double(double x[], char *file, int line)
{
#ifndef VF_NO_MPI
  int    i, *index, proc, err;
  double *xtmp;
  VFenclosure *e=VF_CurrentEnclosure();

  if (VFLIB_Size==1) return;
  VF_GetDPbuffer2_ptr(&xtmp);
  for (proc=0; proc<VFLIB_Size; proc++) {
    if (e->comm_plan.partners[proc].cnt>0) {
      if (proc==VFLIB_Rank) {
        index = e->comm_plan.partners[proc].global_index;
        for (i=0; i<e->comm_plan.partners[proc].cnt; i++) {
          xtmp[i] = x[*index++];
        }
      }
      err = MPI_Bcast(xtmp,e->comm_plan.partners[proc].cnt,
                      MPI_DOUBLE,proc,VFLIB_Comm);
      if (err) {
        fprintf(stderr,"MPI error in VF_Exchange_Double() in %s, line %d\n",file,line);
        VF_Exit(1);
      }
      if (proc!=VFLIB_Rank) {
        index = e->comm_plan.partners[proc].global_index;
        for (i=0; i<e->comm_plan.partners[proc].cnt; i++) {
          x[*index++] = xtmp[i];
        }
      }
    }
  }
#endif
}

void VF_Gather_VFs(float vf[], int root, char *file, int line)
{
#ifndef VF_NO_MPI
  static MPI_Request *requests=NULL;
  int    i, proc, tag, err, *index;
  float  *tmp;
  MPI_Status  status;
  VFenclosure *e=VF_CurrentEnclosure();

  if (VFLIB_Size==1) return;
  if (requests==NULL) {
    requests = (MPI_Request *)VF_Newv(VFLIB_Size*sizeof(MPI_Request));
  }
  VF_GetSPbuffer2_ptr(&tmp);
  /*=======================*/
  /* POST ALL THE RECEIVES */
  /*=======================*/
  if (root==VFLIB_Rank) {
    for (proc=0; proc<VFLIB_Size; proc++) {
      if (proc==root || e->comm_plan.partners[proc].cnt<=0) continue;
      tag = proc+10000;
      err = MPI_Irecv(&tmp[e->comm_plan.partners[proc].offset],
                      e->comm_plan.partners[proc].cnt,MPI_FLOAT,
                      proc,tag,VFLIB_Comm,&requests[proc]);
      if (err) {
        fprintf(stderr,"MPI error in VF_Gather_VFs():MPI_Irecv() in %s, line %d\n",file,line);
        VF_Exit(1);
      }
    }
    /*=====================================*/
    /* WAIT FOR ALL THE RECEIVES TO FINISH */
    /*=====================================*/
    for (proc=0; proc<VFLIB_Size; proc++) {
      if (proc==root || e->comm_plan.partners[proc].cnt<=0) continue;
      if (MPI_Wait(&requests[proc],&status)) {
        fprintf(stderr,"MPI error in VF_Gather_VFs():MPI_Wait() in %s, line %d\n",file,line);
        VF_Exit(1);
      }
    }
    /*==============================*/
    /* REDISTRIBUTE THE VIEWFACTORS */
    /*==============================*/
    for (i=0; i<e->comm_plan.partners[root].cnt; i++) {
      tmp[e->comm_plan.partners[root].offset+i] = vf[i];
    }
    for (proc=0; proc<VFLIB_Size; proc++) {
      if (e->comm_plan.partners[proc].cnt<=0) continue;
      index = e->comm_plan.partners[proc].global_index;
      for (i=0; i<e->comm_plan.partners[proc].cnt; i++) {
        vf[index[i]] = tmp[e->comm_plan.partners[proc].offset+i];
      }
    }
  } else {
    /*====================*/
    /* POST ALL THE SENDS */
    /*====================*/
    if (e->comm_plan.partners[VFLIB_Rank].cnt>0) {
      tag = VFLIB_Rank+10000;
      err = MPI_Send(vf,e->comm_plan.partners[VFLIB_Rank].cnt,
                     MPI_FLOAT,root,tag,VFLIB_Comm);
      if (err) {
        fprintf(stderr,"MPI error in VF_Gather_VFs():MPI_Send() in %s, line %d\n",file,line);
        VF_Exit(1);
      }
    }
  }
#endif
}

void VF_Gather_ExodusData(double data[], double buffer[], int index[], 
                          int ncnt, int nelem_var, char *file, int line)
{
#ifndef VF_NO_MPI
  int i, j, k, tag, err;
  MPI_Status  status;
    
  if (VFLIB_Size==1) return;
  for (j=1; j<VFLIB_Size; j++) {
    if (VFLIB_Rank==0) {
      err = MPI_Recv(&ncnt,1,MPI_INT,j,10000+j,VFLIB_Comm,&status);
      if (err) {
        fprintf(stderr,"MPI error in VF_Gather_ExodusData():MPI_Recv() in %s, line %d\n",file,line);
        VF_Exit(1);
      }
      if (ncnt>0) {
        for (k=0; k<nelem_var; k++) {
          tag = 100000+k*10000+j;
          err = MPI_Recv(buffer,ncnt,MPI_DOUBLE,j,
                   tag,VFLIB_Comm,&status);
          if (err) {
            fprintf(stderr,"MPI error in VF_Gather_ExodusData():MPI_Recv() in %s, line %d\n",file,line);
            VF_Exit(1);
          }
          tag = 200000+k*10000+j;
          err = MPI_Recv(index,ncnt,MPI_INT,j,
                   tag,VFLIB_Comm,&status);
          if (err) {
            fprintf(stderr,"MPI error in VF_Gather_ExodusData():MPI_Recv() in %s, line %d\n",file,line);
            VF_Exit(1);
          }
          for (i=0; i<ncnt; i++) {
            data[index[i]] = buffer[i];
          }
        }
      }
    } else if (VFLIB_Rank==j) {
      err = MPI_Send(&ncnt,1,MPI_INT,0,10000+j,VFLIB_Comm);
      if (err) {
        fprintf(stderr,"MPI error in VF_Gather_ExodusData():MPI_Send() in %s, line %d\n",file,line);
        VF_Exit(1);
      }
      if (ncnt>0) {
        for (k=0; k<nelem_var; k++) {
          tag = 100000+k*10000+j;
          err = MPI_Send(buffer,ncnt,MPI_DOUBLE,0,tag,VFLIB_Comm);
          if (err) {
            fprintf(stderr,"MPI error in VF_Gather_ExodusData():MPI_Send() in %s, line %d\n",file,line);
            VF_Exit(1);
          }
          tag = 200000+k*10000+j;
          err = MPI_Send(index,ncnt,MPI_INT,0,tag,VFLIB_Comm);
          if (err) {
            fprintf(stderr,"MPI error in VF_Gather_ExodusData():MPI_Send() in %s, line %d\n",file,line);
            VF_Exit(1);
          }
        }
      }
    }
    VF_Sync();
  }  
#endif
}

void VF_Gather_RadiosityVectors(double tsurf[], double eps[], double radq[], 
                                double tsurf_g[], double eps_g[], double radq_g[])
{
  int i, j;
  VFenclosure *enclosure=VF_CurrentEnclosure();
#ifndef VF_NO_MPI
  int    proc, size, *index;
  double *tdata, *edata, *rdata;
#endif
    
#ifndef VF_NO_MPI
  if (VFLIB_Size==1) {
    for (i=0; i<enclosure->host_npatches; i++) {
      j          = enclosure->host2vflib_map[i];
      tsurf_g[j] = tsurf[i];
      eps_g[j]   = eps[i];
      radq_g[j]  = radq[i];
    }
  } else {
    for (proc=0; proc<VFLIB_Size; proc++) {
      if (proc==VFLIB_Rank) {
        size  = enclosure->host_npatches;
        index = enclosure->host2vflib_map;
        tdata = tsurf;
        edata = eps;
        rdata = radq;
      } else {
        VF_GetINTbuffer_ptr(&index);
        VF_GetDPbuffer0_ptr(&tdata);
        VF_GetDPbuffer1_ptr(&edata);
        VF_GetDPbuffer2_ptr(&rdata);
      }
      VF_BroadcastInt(&size,1,proc);
      VF_BroadcastInt(index,size,proc);
      VF_BroadcastDouble(tdata,size,proc);
      VF_BroadcastDouble(edata,size,proc);
      VF_BroadcastDouble(rdata,size,proc);
      for (i=0; i<size; i++) {
        j          = index[i];
        tsurf_g[j] = tdata[i];
        eps_g[j]   = edata[i];
        radq_g[j]  = rdata[i];
      }
    }
  }  
#else
  for (i=0; i<enclosure->host_npatches; i++) {
    j          = enclosure->host2vflib_map[i];
    tsurf_g[j] = tsurf[i];
    eps_g[j]   = eps[i];
    radq_g[j]  = radq[i];
  }
#endif
}

void VF_Gather_RadiosityVectorsAux(double tsurf[], double eps[], double radq[], double radj[], 
                                   double tsurf_g[], double eps_g[], double radq_g[], double radj_g[])
{
  int i, j;
  VFenclosure *enclosure=VF_CurrentEnclosure();
#ifndef VF_NO_MPI
  int    proc, size, *index;
  double *tdata, *edata, *qdata, *jdata;
#endif
    
#ifndef VF_NO_MPI
  if (VFLIB_Size==1) {
    for (i=0; i<enclosure->host_npatches; i++) {
      j          = enclosure->host2vflib_map[i];
      tsurf_g[j] = tsurf[i];
      eps_g[j]   = eps[i];
      radq_g[j]  = radq[i];
      radj_g[j]  = radj[i];
    }
  } else {
    for (proc=0; proc<VFLIB_Size; proc++) {
      if (proc==VFLIB_Rank) {
        size  = enclosure->host_npatches;
        index = enclosure->host2vflib_map;
        tdata = tsurf;
        edata = eps;
        qdata = radq;
        jdata = radj;
      } else {
        VF_GetINTbuffer_ptr(&index);
        VF_GetDPbuffer0_ptr(&tdata);
        VF_GetDPbuffer1_ptr(&edata);
        VF_GetDPbuffer2_ptr(&qdata);
        VF_GetDPbuffer3_ptr(&jdata);
      }
      VF_BroadcastInt(&size,1,proc);
      VF_BroadcastInt(index,size,proc);
      VF_BroadcastDouble(tdata,size,proc);
      VF_BroadcastDouble(edata,size,proc);
      VF_BroadcastDouble(qdata,size,proc);
      VF_BroadcastDouble(jdata,size,proc);
      for (i=0; i<size; i++) {
        j          = index[i];
        tsurf_g[j] = tdata[i];
        eps_g[j]   = edata[i];
        radq_g[j]  = qdata[i];
        radj_g[j]  = jdata[i];
      }
    }
  }  
#else
  for (i=0; i<enclosure->host_npatches; i++) {
    j          = enclosure->host2vflib_map[i];
    tsurf_g[j] = tsurf[i];
    eps_g[j]   = eps[i];
    radq_g[j]  = radq[i];
    radj_g[j]  = radj[i];
  }
#endif
}

void VF_Gather_ExodusVals(double J[], double radq[], double J_g[], double radq_g[])
{
  int i, j;
  VFenclosure *enclosure=VF_CurrentEnclosure();
#ifndef VF_NO_MPI
  int    proc, size, *index;
  double *Jdata, *rdata;
#endif
    
#ifndef VF_NO_MPI
  if (VFLIB_Size==1) {
    for (i=0; i<enclosure->host_npatches; i++) {
      j          = enclosure->host2vflib_map[i];
      radq_g[j]  = radq[i];
      J_g[j]     = J[i];
    }
  } else {
    for (proc=0; proc<VFLIB_Size; proc++) {
      if (proc==VFLIB_Rank) {
        size  = enclosure->host_npatches;
        index = enclosure->host2vflib_map;
        Jdata = J;
        rdata = radq;
      } else {
        VF_GetINTbuffer_ptr(&index);
        VF_GetDPbuffer2_ptr(&rdata);
        VF_GetDPbuffer3_ptr(&Jdata);
      }
      VF_BroadcastInt(&size,1,proc);
      VF_BroadcastInt(index,size,proc);
      VF_BroadcastDouble(Jdata,size,proc);
      VF_BroadcastDouble(rdata,size,proc);
      for (i=0; i<size; i++) {
        j          = index[i];
        J_g[j] = Jdata[i];
        radq_g[j]  = rdata[i];
      }
    }
  }  
#else
  for (i=0; i<enclosure->host_npatches; i++) {
    j         = enclosure->host2vflib_map[i];
    J_g[j]    = J[i];
    radq_g[j] = radq[i];
  }
#endif
}
