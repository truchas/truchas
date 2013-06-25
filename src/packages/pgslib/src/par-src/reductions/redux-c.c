/* This file contains the parallel support for PGSLib reductions routines */

/* $Id: redux-c.c,v 1.1.1.1 2000/10/11 22:44:23 ferrell Exp $ */

/* The routines contained in this file are:
   pgslib_global_min_{int,real,double}_c
   pgslib_global_max_{int,real,double}_c
   pgslib_global_sum_{int,real,double}_c
   pgslib_global_all_{log}_c
   pgslib_global_any_{log}_c
*/

#include <stdlib.h>
#include <stdio.h>
#include "pgslib-include-c.h"

#ifdef USE_SGI_SHMEM_LIB
#include "sm.h"
#endif

#include "redux-c.h"

static char sout[2048];
static int  counter = 0;

void pgslib_global_min_int_c(sptr)
     int *sptr;
{
  int rcvbuf;
  int count=1;

#ifdef USE_SGI_SHMEM_LIB
  rcvbuf = sm_int_min_all(*sptr);
#else  
  MPI_Allreduce(sptr, &rcvbuf, count, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
#ifdef DEBUG_REDUX
  sprintf(sout, "counter, *sptr, rcvbuf = %d, %d, %d", counter, *sptr, rcvbuf);
  pgslib_output_c(sout);
  counter++;
#endif  
  *sptr = rcvbuf;
  return;
}

void pgslib_global_min_float_c(sptr)
     float *sptr;
{
  float rcvbuf;
  int count=1;

#ifdef USE_SGI_SHMEM_LIB
  rcvbuf = sm_float_min_all(*sptr);
#else  
  MPI_Allreduce(sptr, &rcvbuf, count, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif
#ifdef DEBUG_REDUX
  sprintf(sout, "counter, *sptr, rcvbuf = %d, %12.4g, %12.4g", counter, *sptr, rcvbuf);
  pgslib_output_c(sout);
  counter++;
#endif  
  *sptr = rcvbuf;
  return;
}

void pgslib_global_min_double_c(sptr)
     double *sptr;
{
  double rcvbuf;
  int count=1;

#ifdef USE_SGI_SHMEM_LIB
  rcvbuf = sm_double_min_all(*sptr);
#else  
  MPI_Allreduce(sptr, &rcvbuf, count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
#ifdef DEBUG_REDUX
  sprintf(sout, "counter, *sptr, rcvbuf = %d, %12.4g, %12.4g", counter, *sptr, rcvbuf);
  pgslib_output_c(sout);
  counter++;
#endif  
  *sptr = rcvbuf;
  return;
}

void pgslib_global_max_int_c(sptr)
     int *sptr;
{
  int rcvbuf;
  int count=1;

#ifdef USE_SGI_SHMEM_LIB
  rcvbuf = sm_int_max_all(*sptr);
#else  
  MPI_Allreduce(sptr, &rcvbuf, count, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
#ifdef DEBUG_REDUX
  sprintf(sout, "counter, *sptr, rcvbuf = %d, %d, %d", counter, *sptr, rcvbuf);
  pgslib_output_c(sout);
  counter++;
#endif  
  *sptr = rcvbuf;
  return;
}

void pgslib_global_max_float_c(sptr)
     float *sptr;
{
  float rcvbuf;
  int count=1;

#ifdef USE_SGI_SHMEM_LIB
  rcvbuf = sm_float_max_all(*sptr);
#else  
  MPI_Allreduce(sptr, &rcvbuf, count, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
#endif
#ifdef DEBUG_REDUX
  sprintf(sout, "counter, *sptr, rcvbuf = %d, %12.4g, %12.4g", counter, *sptr, rcvbuf);
  pgslib_output_c(sout);
  counter++;
#endif  
  *sptr = rcvbuf;
  return;
}

void pgslib_global_max_double_c(sptr)
     double *sptr;
{
  double rcvbuf;
  int count=1;

#ifdef USE_SGI_SHMEM_LIB
  rcvbuf = sm_double_max_all(*sptr);
#else  
  MPI_Allreduce(sptr, &rcvbuf, count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
#ifdef DEBUG_REDUX
  sprintf(sout, "counter, *sptr, rcvbuf = %d, %12.4g, %12.4g", counter, *sptr, rcvbuf);
  pgslib_output_c(sout);
  counter++;
#endif  
  *sptr = rcvbuf;
  return;
}


void pgslib_global_sum_int_c(sptr)
     int *sptr;
{
  int rcvbuf;
  int count=1;

#ifdef USE_SGI_SHMEM_LIB
  rcvbuf = sm_int_sum_all(*sptr);
#else
  MPI_Allreduce(sptr, &rcvbuf, count, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
#ifdef DEBUG_REDUX
  sprintf(sout, "counter, *sptr, temp, rcvbuf = %d, %d, %d, %d", counter, *sptr, temp, rcvbuf);
  pgslib_output_c(sout);
  counter++;
#endif  
  *sptr = rcvbuf;
  return;
}

void pgslib_global_sum_float_c(sptr)
     float *sptr;
{
  float rcvbuf;
  int count=1;

#ifdef USE_SGI_SHMEM_LIB
  rcvbuf = sm_float_sum_all(*sptr);
#else
  MPI_Allreduce(sptr, &rcvbuf, count, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif
#ifdef DEBUG_REDUX
  sprintf(sout, "counter, *sptr, temp, rcvbuf = %d, %12.4g, %12.4g, %12.4g", counter, *sptr, temp, rcvbuf);
  pgslib_output_c(sout);
  counter++;
#endif  
  *sptr = rcvbuf;
  return;
}

void pgslib_global_sum_double_c(sptr)
     double *sptr;
{
  double rcvbuf;
  int count=1;

#ifdef USE_SGI_SHMEM_LIB
  rcvbuf = sm_double_sum_all(*sptr);
#else
  MPI_Allreduce(sptr, &rcvbuf, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#ifdef DEBUG_REDUX
  sprintf(sout, "counter, *sptr, rcvbuf = %d, %12.4g, %12.4g", counter, *sptr, rcvbuf);
  pgslib_output_c(sout);
  counter++;
#endif  
  *sptr = rcvbuf;
  return;
}


void pgslib_global_all_log_c(sptr)
     int *sptr;
{
  int rcvbuf;
  int count=1;

  MPI_Allreduce(sptr, &rcvbuf, count, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  /* Return TRUE or FALSE for Fortran */
  *sptr = rcvbuf ? TRUE : FALSE;
  return;
}


void pgslib_global_any_log_c(sptr)
     int *sptr;
{
  int rcvbuf;
  int count=1;
  char out_string[1024];

  MPI_Allreduce(sptr, &rcvbuf, count, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
  /* Return TRUE or FALSE for Fortran */
  *sptr = rcvbuf ? TRUE : FALSE;

#ifdef DEBUG_REDUX
  sprintf(out_string,"   global_any_log_c, counter = %d, *sptr = %d, rcvbuf = %d\n", counter, *sptr, rcvbuf);
  pgslib_output_c(out_string);
  counter++;
#endif
  return;
}


void pgslib_global_minloc_int_c(minv_ptr, glbl_index_ptr)
     int *minv_ptr;
     int *glbl_index_ptr;
{     
  struct {
    int	value;
    int index;
  } in, out;
  int count = 1;

  in.value = *minv_ptr;
  in.index = *glbl_index_ptr;

  MPI_Allreduce( &in, &out, count, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);

  *minv_ptr       = out.value;
  *glbl_index_ptr = out.index;

  return;
}

void pgslib_global_minloc_float_c(minv_ptr, glbl_index_ptr)
     float *minv_ptr;
     int   *glbl_index_ptr;
{     
  struct {
    float value;
    int   index;
  } in, out;
  int count = 1;

  in.value = *minv_ptr;
  in.index = *glbl_index_ptr;

  MPI_Allreduce( &in, &out, count, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);

  *minv_ptr       = out.value;
  *glbl_index_ptr = out.index;

  return;
}
void pgslib_global_minloc_double_c(minv_ptr, glbl_index_ptr)
     double *minv_ptr;
     int   *glbl_index_ptr;
{     
  struct {
    double value;
    int   index;
  } in, out;
  int count = 1;

  in.value = *minv_ptr;
  in.index = *glbl_index_ptr;

  MPI_Allreduce( &in, &out, count, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

  *minv_ptr       = out.value;
  *glbl_index_ptr = out.index;

  return;
}

void pgslib_global_maxloc_int_c(minv_ptr, glbl_index_ptr)
     int *minv_ptr;
     int *glbl_index_ptr;
{     
  struct {
    int	value;
    int index;
  } in, out;
  int count = 1;

  in.value = *minv_ptr;
  in.index = *glbl_index_ptr;

  MPI_Allreduce( &in, &out, count, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);

  *minv_ptr       = out.value;
  *glbl_index_ptr = out.index;

  return;
}

void pgslib_global_maxloc_float_c(minv_ptr, glbl_index_ptr)
     float *minv_ptr;
     int   *glbl_index_ptr;
{     
  struct {
    float value;
    int   index;
  } in, out;
  int count = 1;

  in.value = *minv_ptr;
  in.index = *glbl_index_ptr;

  MPI_Allreduce( &in, &out, count, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  *minv_ptr       = out.value;
  *glbl_index_ptr = out.index;

  return;
}
void pgslib_global_maxloc_double_c(minv_ptr, glbl_index_ptr)
     double *minv_ptr;
     int   *glbl_index_ptr;
{     
  struct {
    double value;
    int   index;
  } in, out;
  int count = 1;

  in.value = *minv_ptr;
  in.index = *glbl_index_ptr;

  MPI_Allreduce( &in, &out, count, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  *minv_ptr       = out.value;
  *glbl_index_ptr = out.index;

  return;
}
