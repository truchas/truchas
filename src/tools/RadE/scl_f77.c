#ifdef USE_MPI
#include <stdlib.h>
#include "mpi.h"
MPI_Comm my_comm;
static int *displs=NULL;
#endif
static int my_rank;
static int my_size;

#include <FortranCInterface_names.h>

#define MPI_Bcast_int_scalar      TR_ROUTINE_GLOBAL_(mpi_bcast_int_scalar,MPI_BCAST_INT_SCALAR)
#define MPI_Bcast_int_vector      TR_ROUTINE_GLOBAL_(mpi_bcast_int_vector,MPI_BCAST_INT_VECTOR)
#define MPI_Bcast_logical_scalar  TR_ROUTINE_GLOBAL_(mpi_bcast_logical_scalar,MPI_BCAST_LOGICAL_SCALAR)
#define MPI_Bcast_logical_vector  TR_ROUTINE_GLOBAL_(mpi_bcast_logical_vector,MPI_BCAST_LOGICAL_VECTOR)
#define MPI_Bcast_float_scalar    TR_ROUTINE_GLOBAL_(mpi_bcast_float_scalar,MPI_BCAST_FLOAT_SCALAR)
#define MPI_Bcast_float_vector    TR_ROUTINE_GLOBAL_(mpi_bcast_float_vector,MPI_BCAST_FLOAT_VECTOR)
#define MPI_Bcast_double_scalar   TR_ROUTINE_GLOBAL_(mpi_bcast_double_scalar,MPI_BCAST_DOUBLE_SCALAR)
#define MPI_Bcast_double_vector   TR_ROUTINE_GLOBAL_(mpi_bcast_double_vector,MPI_BCAST_DOUBLE_VECTOR)
#define MPI_Bcast_char_vector     TR_ROUTINE_GLOBAL_(mpi_bcast_char_vector,MPI_BCAST_CHAR_VECTOR)
#define MPI_Gather_int            TR_ROUTINE_GLOBAL_(mpi_gather_int,MPI_GATHER_INT)
#define MPI_Gatherv_int           TR_ROUTINE_GLOBAL_(mpi_gatherv_int,MPI_GATHERV_INT)
#define MPI_Gather_float          TR_ROUTINE_GLOBAL_(mpi_gather_float,MPI_GATHER_FLOAT)
#define MPI_Gatherv_float         TR_ROUTINE_GLOBAL_(mpi_gatherv_float,MPI_GATHERV_FLOAT)
#define MPI_Gather_double         TR_ROUTINE_GLOBAL_(mpi_gather_double,MPI_GATHER_DOUBLE)
#define MPI_Gatherv_double        TR_ROUTINE_GLOBAL_(mpi_gatherv_double,MPI_GATHERV_DOUBLE)
#define MPI_Allgather_int         TR_ROUTINE_GLOBAL_(mpi_allgather_int,MPI_ALLGATHER_INT)
#define MPI_Allgatherv_int        TR_ROUTINE_GLOBAL_(mpi_allgatherv_int,MPI_ALLGATHERV_INT)
#define MPI_Allgather_float       TR_ROUTINE_GLOBAL_(mpi_allgather_float,MPI_ALLGATHER_FLOAT)
#define MPI_Allgatherv_float      TR_ROUTINE_GLOBAL_(mpi_allgatherv_float,MPI_ALLGATHERV_FLOAT)
#define MPI_Scatter_int           TR_ROUTINE_GLOBAL_(mpi_scatter_int,MPI_SCATTER_INT)
#define MPI_Scatterv_int          TR_ROUTINE_GLOBAL_(mpi_scatterv_int,MPI_SCATTERV_INT)
#define MPI_Scatter_float         TR_ROUTINE_GLOBAL_(mpi_scatter_float,MPI_SCATTER_FLOAT)
#define MPI_Scatterv_float        TR_ROUTINE_GLOBAL_(mpi_scatterv_float,MPI_SCATTERV_FLOAT)
#define MPI_Scatter_double        TR_ROUTINE_GLOBAL_(mpi_scatter_double,MPI_SCATTER_DOUBLE)
#define MPI_Scatterv_double       TR_ROUTINE_GLOBAL_(mpi_scatterv_double,MPI_SCATTERV_DOUBLE)
#define MPI_Allreduce_sum_int     TR_ROUTINE_GLOBAL_(mpi_allreduce_sum_int,MPI_ALLREDUCE_SUM_INT)
#define MPI_Allreduce_sum_float   TR_ROUTINE_GLOBAL_(mpi_allreduce_sum_float,MPI_ALLREDUCE_SUM_FLOAT)
#define MPI_Allreduce_sum_double  TR_ROUTINE_GLOBAL_(mpi_allreduce_sum_double,MPI_ALLREDUCE_SUM_DOUBLE)
#define MPI_Allreduce_max_int     TR_ROUTINE_GLOBAL_(mpi_allreduce_max_int,MPI_ALLREDUCE_MAX_INT)
#define MPI_Allreduce_max_float   TR_ROUTINE_GLOBAL_(mpi_allreduce_max_float,MPI_ALLREDUCE_MAX_FLOAT)
#define MPI_Allreduce_max_double  TR_ROUTINE_GLOBAL_(mpi_allreduce_max_double,MPI_ALLREDUCE_MAX_DOUBLE)
#define MPI_Send_int              TR_ROUTINE_GLOBAL_(mpi_send_int,MPI_SEND_INT)
#define MPI_Send_float            TR_ROUTINE_GLOBAL_(mpi_send_float,MPI_SEND_FLOAT)
#define MPI_Send_double           TR_ROUTINE_GLOBAL_(mpi_send_double,MPI_SEND_DOUBLE)
#define MPI_Recv_int              TR_ROUTINE_GLOBAL_(mpi_recv_int,MPI_RECV_INT)
#define MPI_Recv_float            TR_ROUTINE_GLOBAL_(mpi_recv_float,MPI_RECV_FLOAT)
#define MPI_Recv_double           TR_ROUTINE_GLOBAL_(mpi_recv_double,MPI_RECV_DOUBLE)

void TR_ROUTINE_GLOBAL(initialize,INITIALIZE) ()
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

void TR_ROUTINE_GLOBAL_(finalize,FINALIZE) ()
{
#ifdef USE_MPI
  free(displs);
  MPI_Finalize();
#endif
}

int TR_ROUTINE_GLOBAL_(f77_comm_size,F77_COMM_SIZE) ()
{
  return my_size;
}

int TR_ROUTINE_GLOBAL(f77_comm_rank,F77_COMM_RANK) ()
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

void MPI_Bcast_char_vector(char *buffer, int *root, int len)
{
#ifdef USE_MPI
  MPI_Bcast(buffer,len,MPI_CHAR,*root,my_comm);
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

