#ifdef USE_MPI
#include <stdlib.h>
#include "mpi.h"
MPI_Comm my_comm;
static int *displs=NULL;
#endif
static int my_rank;
static int my_size;

void initialize()
{
#ifdef USE_MPI
  int flag;
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(NULL,NULL);
  MPI_Comm_dup  (MPI_COMM_WORLD, &my_comm);
  MPI_Comm_rank (my_comm, &my_rank);
  MPI_Comm_size (my_comm, &my_size);
  displs = (int*)malloc(sizeof(int)*my_size);
#else
  my_rank = 0;
  my_size = 1;
#endif
}

void finalize()
{
#ifdef USE_MPI
  free(displs);
  MPI_Finalize();
#endif
}

int f77_comm_size()
{
  return my_size;
}

int f77_comm_rank()
{
  return my_rank;
}

void MPI_Send_int(int *buf, int *count, int *dest, int *tag)
{
#ifdef USE_MPI
  MPI_Send(buf,*count,MPI_INT,*dest,*tag,my_comm);
#endif
}

void MPI_Send_float(float *buf, int *count, int *dest, int *tag)
{
#ifdef USE_MPI
  MPI_Send(buf,*count,MPI_FLOAT,*dest,*tag,my_comm);
#endif
}

void MPI_Send_double(double *buf, int *count, int *dest, int *tag)
{
#ifdef USE_MPI
  MPI_Send(buf,*count,MPI_DOUBLE,*dest,*tag,my_comm);
#endif
}

void MPI_Recv_int(int *buf, int *count, int *source, int *tag)
{
#ifdef USE_MPI
  MPI_Status status;
  MPI_Recv(buf,*count,MPI_INT,*source,*tag,my_comm,&status);
#endif
}

void MPI_Recv_float(float *buf, int *count, int *source, int *tag)
{
#ifdef USE_MPI
  MPI_Status status;
  MPI_Recv(buf,*count,MPI_FLOAT,*source,*tag,my_comm,&status);
#endif
}

void MPI_Recv_double(double *buf, int *count, int *source, int *tag)
{
#ifdef USE_MPI
  MPI_Status status;
  MPI_Recv(buf,*count,MPI_DOUBLE,*source,*tag,my_comm,&status);
#endif
}

void MPI_Bcast_int_scalar(int *buffer, int *root)
{
#ifdef USE_MPI
  MPI_Bcast(buffer,1,MPI_INT,*root,my_comm);
#endif
}

void MPI_Bcast_int_vector(int *buffer, int *len, int *root)
{
#ifdef USE_MPI
  MPI_Bcast(buffer,*len,MPI_INT,*root,my_comm);
#endif
}

void MPI_Bcast_logical_scalar(int *buffer, int *root)
{
#ifdef USE_MPI
  MPI_Bcast(buffer,1,MPI_INT,*root,my_comm);
#endif
}

void MPI_Bcast_logical_vector(int *buffer, int *len, int *root)
{
#ifdef USE_MPI
  MPI_Bcast(buffer,*len,MPI_INT,*root,my_comm);
#endif
}

void MPI_Bcast_float_scalar(float *buffer, int *root)
{
#ifdef USE_MPI
  MPI_Bcast(buffer,1,MPI_FLOAT,*root,my_comm);
#endif
}

void MPI_Bcast_float_vector(float *buffer, int *len, int *root)
{
#ifdef USE_MPI
  MPI_Bcast(buffer,*len,MPI_FLOAT,*root,my_comm);
#endif
}

void MPI_Bcast_double_scalar(double *buffer, int *root)
{
#ifdef USE_MPI
  MPI_Bcast(buffer,1,MPI_DOUBLE,*root,my_comm);
#endif
}

void MPI_Bcast_double_vector(double *buffer, int *len, int *root)
{
#ifdef USE_MPI
  MPI_Bcast(buffer,*len,MPI_DOUBLE,*root,my_comm);
#endif
}

void MPI_Bcast_char_vector(char *buffer, int *len, int *root)
{
#ifdef USE_MPI
  MPI_Bcast(buffer,*len,MPI_CHAR,*root,my_comm);
#endif
}

void MPI_Gather_int(int *sendbuf, int *recvbuf, int *root)
{
#ifdef USE_MPI
  MPI_Gather(sendbuf,1,MPI_INT, recvbuf,1,MPI_INT, *root,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Gatherv_int(int *sendbuf, int *sendcount, int *recvbuf, int *recvcounts, int *root)
{
  int i;
#ifdef USE_MPI
  if (my_rank==*root) {
    for (displs[0]=0, i=1; i<my_size; i++) {
      displs[i] = displs[i-1] + recvcounts[i-1];
    }
  }
  MPI_Gatherv(sendbuf,*sendcount,MPI_INT, recvbuf,recvcounts,displs,MPI_INT, *root,my_comm);
#else
  for (i=0; i<*sendcount; i++) {
    recvbuf[i] = sendbuf[i];
  }
#endif
}

void MPI_Gather_float(float *sendbuf, float *recvbuf, int *root)
{
#ifdef USE_MPI
  MPI_Gather(sendbuf,1,MPI_FLOAT, recvbuf,1,MPI_FLOAT, *root,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Gatherv_float(float *sendbuf, int *sendcount, float *recvbuf, int *recvcounts, int *root)
{
  int i;
#ifdef USE_MPI
  if (my_rank==*root) {
    for (displs[0]=0, i=1; i<my_size; i++) {
      displs[i] = displs[i-1] + recvcounts[i-1];
    }
  }
  MPI_Gatherv(sendbuf,*sendcount,MPI_FLOAT, recvbuf,recvcounts,displs,MPI_FLOAT, *root,my_comm);
#else
  for (i=0; i<*sendcount; i++) {
    recvbuf[i] = sendbuf[i];
  }
#endif
}

void MPI_Gather_double(double *sendbuf, double *recvbuf, int *root)
{
#ifdef USE_MPI
  MPI_Gather(sendbuf,1,MPI_DOUBLE, recvbuf,1,MPI_DOUBLE, *root,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Gatherv_double(double *sendbuf, int *sendcount, double *recvbuf, int *recvcounts, int *root)
{
  int i;
#ifdef USE_MPI
  if (my_rank==*root) {
    for (displs[0]=0, i=1; i<my_size; i++) {
      displs[i] = displs[i-1] + recvcounts[i-1];
    }
  }
  MPI_Gatherv(sendbuf,*sendcount,MPI_DOUBLE, recvbuf,recvcounts,displs,MPI_DOUBLE, *root,my_comm);
#else
  for (i=0; i<*sendcount; i++) {
    recvbuf[i] = sendbuf[i];
  }
#endif
}

void MPI_Allgather_int(int *sendbuf, int *recvbuf)
{
#ifdef USE_MPI
  MPI_Allgather(sendbuf,1,MPI_INT, recvbuf,1,MPI_INT, my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Allgatherv_int(int *sendbuf, int *sendcount, int *recvbuf, int *recvcounts)
{
  int i;
#ifdef USE_MPI
  for (displs[0]=0, i=1; i<my_size; i++) {
    displs[i] = displs[i-1] + recvcounts[i-1];
  }
  MPI_Allgatherv(sendbuf,*sendcount,MPI_INT, recvbuf,recvcounts,displs,MPI_INT, my_comm);
#else
  for (i=0; i<*sendcount; i++) {
    recvbuf[i] = sendbuf[i];
  }
#endif
}

void MPI_Allgather_float(float *sendbuf, float *recvbuf)
{
#ifdef USE_MPI
  MPI_Allgather(sendbuf,1,MPI_FLOAT, recvbuf,1,MPI_FLOAT, my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Allgatherv_float(float *sendbuf, int *sendcount, float *recvbuf, int *recvcounts)
{
  int i;
#ifdef USE_MPI
  for (displs[0]=0, i=1; i<my_size; i++) {
    displs[i] = displs[i-1] + recvcounts[i-1];
  }
  MPI_Allgatherv(sendbuf,*sendcount,MPI_FLOAT, recvbuf,recvcounts,displs,MPI_FLOAT, my_comm);
#else
  for (i=0; i<*sendcount; i++) {
    recvbuf[i] = sendbuf[i];
  }
#endif
}

void MPI_Scatter_int(int *sendbuf, int *recvbuf, int *root)
{
#ifdef USE_MPI
  MPI_Scatter(sendbuf,1,MPI_INT, recvbuf,1,MPI_INT, *root,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Scatterv_int(int *sendbuf, int *sendcounts, int *recvbuf, int *recvcount, int *root)
{
  int i;
#ifdef USE_MPI
  if (my_rank==*root) {
    for (displs[0]=0, i=1; i<my_size; i++) {
      displs[i] = displs[i-1] + sendcounts[i-1];
    }
  }
  MPI_Scatterv(sendbuf,sendcounts,displs,MPI_INT, recvbuf,*recvcount,MPI_INT, *root,my_comm);
#else
  for (i=0; i<*recvcount; i++) {
    recvbuf[i] = sendbuf[i];
  }
#endif
}

void MPI_Scatter_float(float *sendbuf, float *recvbuf, int *root)
{
#ifdef USE_MPI
  MPI_Scatter(sendbuf,1,MPI_FLOAT, recvbuf,1,MPI_FLOAT, *root,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Scatterv_float(float *sendbuf, int *sendcounts, float *recvbuf, int *recvcount, int *root)
{
  int i;
#ifdef USE_MPI
  if (my_rank==*root) {
    for (displs[0]=0, i=1; i<my_size; i++) {
      displs[i] = displs[i-1] + sendcounts[i-1];
    }
  }
  MPI_Scatterv(sendbuf,sendcounts,displs,MPI_FLOAT, recvbuf,*recvcount,MPI_FLOAT, *root,my_comm);
#else
  for (i=0; i<*recvcount; i++) {
    recvbuf[i] = sendbuf[i];
  }
#endif
}

void MPI_Scatter_double(double *sendbuf, double *recvbuf, int *root)
{
#ifdef USE_MPI
  MPI_Scatter(sendbuf,1,MPI_DOUBLE, recvbuf,1,MPI_DOUBLE, *root,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Scatterv_double(double *sendbuf, int *sendcounts, double *recvbuf, int *recvcount, int *root)
{
  int i;
#ifdef USE_MPI
  if (my_rank==*root) {
    for (displs[0]=0, i=1; i<my_size; i++) {
      displs[i] = displs[i-1] + sendcounts[i-1];
    }
  }
  MPI_Scatterv(sendbuf,sendcounts,displs,MPI_DOUBLE, recvbuf,*recvcount,MPI_DOUBLE, *root,my_comm);
#else
  for (i=0; i<*recvcount; i++) {
    recvbuf[i] = sendbuf[i];
  }
#endif
}

void MPI_Allreduce_sum_int(int *sendbuf, int *recvbuf)
{
#ifdef USE_MPI
  MPI_Allreduce(sendbuf,recvbuf,1,MPI_INT,MPI_SUM,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Allreduce_sum_float(float *sendbuf, float *recvbuf)
{
#ifdef USE_MPI
  MPI_Allreduce(sendbuf,recvbuf,1,MPI_FLOAT,MPI_SUM,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Allreduce_sum_double(double *sendbuf, double *recvbuf)
{
#ifdef USE_MPI
  MPI_Allreduce(sendbuf,recvbuf,1,MPI_DOUBLE,MPI_SUM,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Allreduce_max_int(int *sendbuf, int *recvbuf)
{
#ifdef USE_MPI
  MPI_Allreduce(sendbuf,recvbuf,1,MPI_INT,MPI_MAX,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Allreduce_max_float(float *sendbuf, float *recvbuf)
{
#ifdef USE_MPI
  MPI_Allreduce(sendbuf,recvbuf,1,MPI_FLOAT,MPI_MAX,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

void MPI_Allreduce_max_double(double *sendbuf, double *recvbuf)
{
#ifdef USE_MPI
  MPI_Allreduce(sendbuf,recvbuf,1,MPI_DOUBLE,MPI_MAX,my_comm);
#else
  *recvbuf = *sendbuf;
#endif
}

