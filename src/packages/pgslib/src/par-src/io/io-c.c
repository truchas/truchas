/* This file contains the parallel support for PGSLib io suppor routines.*/

/* $Id: io-c.c,v 1.1.1.1 2000/10/11 22:44:23 ferrell Exp $  */


/* The routines contained in this file are:

   pgslib_bcast_{int, real, double, log}_{scalar,vector}_c
   pgslib_bcast_{char                  }_{       vector}_c

   pgslib_dist_{int, real, double, log}_{scalar, vector}_c

*/
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "pgslib-include-c.h"

#include "utility-c.h"
#include "io-c.h"

/* Broadcast Scalars */
void pgslib_bcast_int_scalar_c(scalar)
     int *scalar;
{
  MPI_Bcast(scalar, 1, MPI_INT, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD );
}
void pgslib_bcast_float_scalar_c(scalar)
     float *scalar;
{
  MPI_Bcast(scalar, 1, MPI_FLOAT, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD );
}
void pgslib_bcast_double_scalar_c(scalar)
     double *scalar;
{
  MPI_Bcast(scalar, 1, MPI_DOUBLE, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD );
}

void pgslib_bcast_log_scalar_c(scalar)
     C_LOG_TYPE *scalar; /* what is a FORTRAN LOGICAL in C?*/
{
  MPI_Bcast(scalar, 1*BYTES_PER_LOG, MPI_LOG_TYPE, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD );
}

void pgslib_bcast_char_scalar_c(scalar)
     char *scalar;
{
  MPI_Bcast(scalar, BYTES_PER_CHAR, MPI_BYTE, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD );
}

/* Broadcast Vectors */
void pgslib_bcast_int_vector_c(vector, len)
     int *vector, *len;
{
  MPI_Bcast(vector, *len, MPI_INT, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD );
}
void pgslib_bcast_float_vector_c(vector, len)
     float *vector;
     int *len;
{
  MPI_Bcast(vector, *len, MPI_FLOAT, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD );
}
void pgslib_bcast_double_vector_c(vector, len)
     double *vector;
     int *len;
{
  MPI_Bcast(vector, *len, MPI_DOUBLE, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD );
}

void pgslib_bcast_log_vector_c(vector, len)
     C_LOG_TYPE *vector;
     int        *len;
{
  MPI_Bcast(vector, (*len)*BYTES_PER_LOG, MPI_LOG_TYPE, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD );
}

void pgslib_bcast_char_vector_c(vector, len)
     char *vector;
     int *len;
{ int bytes_to_send;
  char outstring[256];
  int c;

  bytes_to_send = BYTES_PER_CHAR * (*len);
  MPI_Bcast(vector, bytes_to_send, MPI_BYTE, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD );
}

/********************************************************************/
/* Distribute Routines */

/* Distribute Scalars */
void pgslib_dist_int_scalar_c(scalar_out, scalarv_in)
     int *scalar_out, *scalarv_in;
{
  char errstring[256];

#ifdef DEBUG_IO
  sprintf(errstring, "before call, scalarv_in[1,2] = %d, %d, scalar_out = %d\n", *scalarv_in, *(scalarv_in+1), *scalar_out);
  pgslib_output_c(errstring);
#endif

  MPI_Scatter(scalarv_in, 1, MPI_INT, scalar_out, 1, MPI_INT,
	      PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
#ifdef DEBUG_IO
  sprintf(errstring, "after  call, scalarv_in[1,2] = %d, %d, scalar_out = %d\n", *scalarv_in, *(scalarv_in+1), *scalar_out);
  pgslib_output_c(errstring);
#endif
}
void pgslib_dist_float_scalar_c(scalar_out, scalarv_in)
     float *scalar_out, *scalarv_in;
{
  MPI_Scatter(scalarv_in, 1, MPI_FLOAT, scalar_out, 1, MPI_FLOAT,
	      PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}
void pgslib_dist_double_scalar_c(scalar_out, scalarv_in)
     double *scalar_out, *scalarv_in;
{
  MPI_Scatter(scalarv_in, 1, MPI_DOUBLE, scalar_out, 1, MPI_DOUBLE,
	      PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}
	      
void pgslib_dist_log_scalar_c(scalar_out, scalarv_in)
     C_LOG_TYPE *scalar_out, *scalarv_in; /* What is FORTRAN LOGICAL in C?*/
{
  MPI_Scatter(scalarv_in, 1*BYTES_PER_LOG, MPI_LOG_TYPE,
	      scalar_out, 1*BYTES_PER_LOG, MPI_LOG_TYPE,
	      PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}
	      
void pgslib_dist_char_scalar_c(scalar_out, scalarv_in)
     char *scalar_out, *scalarv_in;
{
  MPI_Scatter(scalarv_in, BYTES_PER_CHAR, MPI_BYTE, scalar_out, BYTES_PER_CHAR, MPI_BYTE,
	      PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}
	      
/* Distribute Vectors */
void pgslib_dist_int_vector_c(vector_out, out_len, vector_in, lengths)
     int *vector_out, *vector_in;
     int *out_len, *lengths;
{
  static int *displs;
  static int init_dist_int_vector = 0;
  static int thisPE, nPE;
  int i;
  
  
  if (!init_dist_int_vector) {
    MPI_Comm_size( MPI_COMM_WORLD, &nPE);
    displs = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(displs, "malloc failed in pgslib_dist_int_vector");
    init_dist_int_vector = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisPE);
  }

  if (thisPE == PGSLib_IO_ROOT_PE) {
    displs[0]=0;
    for(i=1;i<nPE;i++) {
      displs[i] = displs[i-1] + lengths[i-1];
    }
  }

  MPI_Scatterv(vector_in, lengths, displs, MPI_INT, 
               vector_out, *out_len, MPI_INT, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}

void pgslib_dist_float_vector_c(vector_out, out_len, vector_in, lengths)
     float *vector_out, *vector_in;
     int *out_len, *lengths;
{
  static int *displs;
  static int init_dist_float_vector = 0;
  static int thisPE, nPE;
  int i;
  
  if (!init_dist_float_vector) {
    MPI_Comm_size( MPI_COMM_WORLD, &nPE);
    displs = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(displs, "malloc failed in pgslib_dist_float_vector");
    init_dist_float_vector = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisPE);
  }

  if (thisPE == PGSLib_IO_ROOT_PE) {
    displs[0]=0;
    for(i=1;i<nPE;i++) {
      displs[i] = displs[i-1] + lengths[i-1];
    }
  }

  MPI_Scatterv(vector_in, lengths, displs, MPI_FLOAT, 
               vector_out, *out_len, MPI_FLOAT, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}

void pgslib_dist_double_vector_c(vector_out, out_len, vector_in, lengths)
     double *vector_out, *vector_in;
     int *out_len, *lengths;
{
  static int *displs;
  static int init_dist_double_vector = 0;
  static int thisPE, nPE;
  int i;
  
  if (!init_dist_double_vector) {
    MPI_Comm_size( MPI_COMM_WORLD, &nPE);
    displs = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(displs, "malloc failed in pgslib_dist_double_vector");
    init_dist_double_vector = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisPE);
  }

  if (thisPE == PGSLib_IO_ROOT_PE) {
    displs[0]=0;
    for(i=1;i<nPE;i++) {
      displs[i] = displs[i-1] + lengths[i-1];
    }
  }

  MPI_Scatterv(vector_in, lengths, displs, MPI_DOUBLE, 
               vector_out, *out_len, MPI_DOUBLE, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}

void pgslib_dist_log_vector_c(vector_out, out_len, vector_in, lengths)
     C_LOG_TYPE *vector_out, *vector_in; /* What is a F90 logical in C?*/
     int        *out_len, *lengths;
{
  static int *displs, *local_lengths;
  static int init_dist_int_vector = 0;
  static int thisPE, nPE, local_out_len;
  int i;
  
  
  if (!init_dist_int_vector) {
    MPI_Comm_size( MPI_COMM_WORLD, &nPE);
    displs = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(displs, "malloc failed in pgslib_dist_log_vector");
    local_lengths = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(local_lengths, "malloc failed in pgslib_dist_log_vector");
    init_dist_int_vector = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisPE);
  }

  if (thisPE == PGSLib_IO_ROOT_PE) {
    displs[0]=0;
    local_lengths[0] = BYTES_PER_LOG * lengths[0];
    for(i=1;i<nPE;i++) {
      local_lengths[i] = BYTES_PER_LOG * lengths[i];
      displs[i] = displs[i-1] + local_lengths[i-1];
    }
  }

  local_out_len = BYTES_PER_LOG * (*out_len);
  MPI_Scatterv(vector_in, local_lengths, displs, MPI_LOG_TYPE,
               vector_out, local_out_len,       MPI_LOG_TYPE,
	       PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}

/**********************************************************************/
/********** Collate Routines ******************************************/

/* Collate scalars */

void pgslib_collate_int_scalar_c(scalarv_out, scalar_in)
     int *scalarv_out, *scalar_in;
{ 
  MPI_Gather(scalar_in, 1, MPI_INT, scalarv_out, 1, MPI_INT, 
	     PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}


void pgslib_collate_float_scalar_c(scalarv_out, scalar_in)
     float *scalarv_out, *scalar_in;
{ 
  MPI_Gather(scalar_in, 1, MPI_FLOAT, scalarv_out, 1, MPI_FLOAT, 
	     PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}


void pgslib_collate_double_scalar_c(scalarv_out, scalar_in)
     double *scalarv_out, *scalar_in;
{ 
  MPI_Gather(scalar_in, 1, MPI_DOUBLE, scalarv_out, 1, MPI_DOUBLE, 
	     PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}


void pgslib_collate_log_scalar_c(scalarv_out, scalar_in)
     C_LOG_TYPE *scalarv_out, *scalar_in;
{ 
  MPI_Gather(scalar_in, 1*BYTES_PER_LOG, MPI_LOG_TYPE,
	     scalarv_out, 1*BYTES_PER_LOG, MPI_LOG_TYPE,
	     PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}


void pgslib_collate_char_scalar_c(scalarv_out, scalar_in)
     char *scalarv_out, *scalar_in;
{ int bytes_to_send;

  bytes_to_send = BYTES_PER_CHAR;
  MPI_Gather(scalar_in, bytes_to_send, MPI_BYTE, scalarv_out, bytes_to_send, MPI_BYTE, 
	     PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}

/* Collate vectors */
void pgslib_collate_int_vector_c(vector_out, lengths, vector_in, in_length)
     int *vector_out, *vector_in;
     int *lengths, *in_length;
{
  static int *displs;
  static int init_collate_int_vector = 0;
  static int thisPE, nPE;
  int i;
  
  
  if (!init_collate_int_vector) {
    MPI_Comm_size( MPI_COMM_WORLD, &nPE);
    displs = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(displs, "malloc failed in pgslib_collate_int_vector");
    init_collate_int_vector = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisPE);
  }

  if (thisPE == PGSLib_IO_ROOT_PE) {
    displs[0]=0;
    for(i=1;i<nPE;i++) {
      displs[i] = displs[i-1] + lengths[i-1];
    }
  }

  MPI_Gatherv(vector_in, *in_length, MPI_INT, 
               vector_out, lengths, displs, MPI_INT, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}

void pgslib_collate_float_vector_c(vector_out, lengths, vector_in, in_length)
     float *vector_out, *vector_in;
     int *lengths, *in_length;
{
  static int *displs;
  static int init_collate_float_vector = 0;
  static int thisPE, nPE;
  int i;
  
  
  if (!init_collate_float_vector) {
    MPI_Comm_size( MPI_COMM_WORLD, &nPE);
    displs = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(displs, "malloc failed in pgslib_collate_float_vector");
    init_collate_float_vector = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisPE);
  }

  if (thisPE == PGSLib_IO_ROOT_PE) {
    displs[0]=0;
    for(i=1;i<nPE;i++) {
      displs[i] = displs[i-1] + lengths[i-1];
    }
  }

  MPI_Gatherv(vector_in, *in_length, MPI_FLOAT, 
               vector_out, lengths, displs, MPI_FLOAT, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}

void pgslib_collate_double_vector_c(vector_out, lengths, vector_in, in_length)
     double *vector_out, *vector_in;
     int *lengths, *in_length;
{
  static int *displs;
  static int init_collate_double_vector = 0;
  static int thisPE, nPE;
  int i;
  
  
  if (!init_collate_double_vector) {
    MPI_Comm_size( MPI_COMM_WORLD, &nPE);
    displs = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(displs, "malloc failed in pgslib_collate_double_vector");
    init_collate_double_vector = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisPE);
  }

  if (thisPE == PGSLib_IO_ROOT_PE) {
    displs[0]=0;
    for(i=1;i<nPE;i++) {
      displs[i] = displs[i-1] + lengths[i-1];
    }
  }

  MPI_Gatherv(vector_in, *in_length, MPI_DOUBLE, 
               vector_out, lengths, displs, MPI_DOUBLE, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}

void pgslib_collate_log_vector_c(vector_out, lengths, vector_in, in_length)
     C_LOG_TYPE *vector_out, *vector_in;
     int *lengths, *in_length;
{
  static int *displs, *local_lengths;
  static int init_collate_log_vector = 0;
  static int thisPE, nPE, local_in_length;
  int i;
  
  
  if (!init_collate_log_vector) {
    MPI_Comm_size( MPI_COMM_WORLD, &nPE);
    displs = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(displs, "malloc failed in pgslib_collate_log_vector");
    local_lengths = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(local_lengths, "malloc failed in pgslib_collate_log_vector");
    init_collate_log_vector = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisPE);
  }

  if (thisPE == PGSLib_IO_ROOT_PE) {
    displs[0]=0;
    local_lengths[0] = BYTES_PER_LOG * lengths[0];
    for(i=1;i<nPE;i++) {
      local_lengths[i] = BYTES_PER_LOG * lengths[i];
      displs[i] = displs[i-1] + local_lengths[i-1];
    }
  }

  local_in_length = BYTES_PER_LOG * (*in_length);
  MPI_Gatherv(vector_in,  local_in_length, MPI_LOG_TYPE, 
              vector_out, local_lengths, displs, MPI_LOG_TYPE,
	      PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}

void pgslib_collate_char_vector_c(vector_out, lengths, vector_in, in_lengths)
     char *vector_out, *vector_in;
     int  *lengths,    *in_lengths;
{ int bytes_to_send;
  static int *displs;
  static int *bytes_to_receive;
  static int init_collate_char_vector = 0;
  static int thisPE, nPE;
  int i;
  char estr[512];
  
  
  if (!init_collate_char_vector) {
    MPI_Comm_size( MPI_COMM_WORLD, &nPE);
    displs = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(displs, "malloc failed in pgslib_collate_char_vector");
    bytes_to_receive = (int *) malloc(sizeof(int)* nPE);
    pgslib_check_malloc_c(bytes_to_receive, "malloc failed in pgslib_collate_char_vector");
    init_collate_char_vector = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisPE);
  }

  if (thisPE == PGSLib_IO_ROOT_PE) {
    displs[0]=0;
    bytes_to_receive[0] = BYTES_PER_CHAR * lengths[0];
    for(i=1;i<nPE;i++) {
      bytes_to_receive[i] = BYTES_PER_CHAR * lengths[i];
      displs[i] = displs[i-1] + bytes_to_receive[i-1];
    }
  }
    
  bytes_to_send = BYTES_PER_CHAR * (*in_lengths);
#ifdef DEBUG_IO
  sprintf(estr, "bytes_to_send = %d", bytes_to_send);
  pgslib_output_c(estr);
  pgslib_flush_output_c();
#endif
    
  MPI_Gatherv(vector_in, bytes_to_send, MPI_BYTE, 
               vector_out, bytes_to_receive, displs, MPI_BYTE, PGSLib_IO_ROOT_PE, MPI_COMM_WORLD);
}

