/*
   This is the generic procedure which gets included in the type specific files.
   void pgslib_constrained_send_rcv(SendS, RcvR)

This routines takes two inputs, a COMM_SEND and a COMM_RCV structure.
  COMM_SEND contains the data for the messages to be send.
  COMM_RCV  holds the received data.  Depending on the particular
            need, COMM_RCV may be pre-allocated, or parts of it
	    may be allocated on the fly.

A barrier at the start of the routine insures that all PEs are ready for this
communication step.

*/

/* $Id: constrained-send-rcv.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <limits.h>

#include "pgslib-types.h"
#include "tags-c.h"

/*#define USE_SEND_RECV_BARRIER*/

static  int nproc = -1, myrank = -1;
static  float tick_tock;

void PGSLIB_ROUTINE_NAME(pgslib_cnstd_send_rcv_)(SendS, RcvR)
COMM_SEND *SendS;
COMM_RCV  *RcvR;
{ 
  int N_Send, N_Rcvd, N_Rcv;
  int s, r, l, mpi_err;
  int TagBits;

  SEND_BUFFER *this_S_Buff;
  MPI_Request *send_request;
  MPI_Status  *send_status;

  RCV_BUFFER *this_R_Buff;
  MPI_Request *rcv_request;
  MPI_Status  *rcv_status;
  char estr[512];
  struct tms t_buff, start_buff, stop_buff;
  clock_t  start, stop;

  if (nproc < 0)
    { MPI_Comm_size( MPI_COMM_WORLD, &nproc );
      MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
      tick_tock = sysconf(_SC_CLK_TCK);
      }

    
  TagBits = constrainedSendRcv_TagBits();

#ifdef USE_SEND_RECV_BARRIER

#ifdef USE_TIMERS_2
  start = times( &start_buff );
#endif

#ifdef DEBUG_SEND_RCV
  pgslib_output_c("Calling initial barrier in CSR");
#endif
  mpi_err = MPI_Barrier( MPI_COMM_WORLD );
#ifdef DEBUG_SEND_RCV
  pgslib_output_c("Passed initial barrier in CSR");
  pgslib_flush_output_c();
#endif

#ifdef USE_TIMERS_2
  stop  = times( &stop_buff );
  barrier_time += (float)(( stop_buff.tms_utime +  stop_buff.tms_stime) - 
		   (start_buff.tms_utime + start_buff.tms_stime)) /  tick_tock;
#endif

#endif	

#ifdef USE_TIMERS_2
  start = times( &start_buff );
#endif

  N_Send = SendS->N_Send;
  N_Rcv  = RcvR->N_Rcv_Buffers;

/* First post the receives, then the sends.  That should keep the network
   as drained as possible */

  rcv_request = RcvR->rcv_request;
  rcv_status  = RcvR->rcv_status;
    
  for(r=0; r<N_Rcv; r++) {
    this_R_Buff = RcvR->R_Buffers[r];
#ifdef DEBUG_SEND_RCV
    sprintf(estr,
	    "CSR, Posting receive %d: n_rcv=%d, PE=%d, TAG=%d, Iteration Tag=%d",
	    r,
	    this_R_Buff->n_max_rcv_words,
	    this_R_Buff->src_PE.PE,
	    this_R_Buff->rcv_TAG,
	    TagBits);
    pgslib_output_c(estr);	
    pgslib_flush_output_c();	
#endif
    MPI_Irecv(&((this_R_Buff->PGSLIB_TYPE_NAME(rcv_data_))[0]),
	      this_R_Buff->n_max_rcv_words,
	      PGSLIB_MPI_DATA_TYPE,
	      this_R_Buff->src_PE.PE,
	      (this_R_Buff->rcv_TAG || TagBits),
	      MPI_COMM_WORLD,
	      &(rcv_request[r]) );
#ifdef DEBUG_SEND_RCV
    sprintf(estr, "posted receive %d", r);
    pgslib_output_c(estr);
  pgslib_flush_output_c();
#endif
  }
  
/* After posting the receives initiate the sends */

  send_request = SendS->send_request;
  send_status  = SendS->send_status;
    
  for(s=0; s<N_Send; s++) {
    this_S_Buff = SendS->S_Buffers[s];
#ifdef DEBUG_SEND_RCV
    sprintf(estr,
	    "CSR, Sending message %d: n_send=%d, PE=%d, TAG=%d, iteration Tag = %d, first word=%d",
	    s,
	    this_S_Buff->n_send_words,
	    this_S_Buff->dest_PE.PE,
	    this_S_Buff->send_TAG,
	    TagBits,
	    (int)((this_S_Buff->PGSLIB_TYPE_NAME(send_data_))[0]));
    pgslib_output_c(estr);
  pgslib_flush_output_c();
#endif
    MPI_Isend(&((this_S_Buff->PGSLIB_TYPE_NAME(send_data_))[0]),
	      this_S_Buff->n_send_words,
	      PGSLIB_MPI_DATA_TYPE,
	      this_S_Buff->dest_PE.PE,
	      (this_S_Buff->send_TAG || TagBits),
	      MPI_COMM_WORLD,
	      &(send_request[s]) );
#ifdef DEBUG_SEND_RCV
    sprintf(estr, "sent message %d",s);
    pgslib_output_c(estr);
  pgslib_flush_output_c();
#endif
  }
  
/* All the sends are posted, wait for receives to complete */
  if (N_Rcv > 0) {
    MPI_Waitall(N_Rcv, &rcv_request[0], &rcv_status[0]);
  }

/* Finally, finish up the send requests.*/
  if (N_Send > 0) {
    MPI_Waitall(N_Send, &send_request[0], &send_status[0]);

  }

#ifdef USE_TIMERS_2
  stop  = times( &stop_buff );
  sr_time += (float)(( stop_buff.tms_utime +  stop_buff.tms_stime) - 
		   (start_buff.tms_utime + start_buff.tms_stime)) /  tick_tock;
#endif

#ifdef USE_SEND_RECV_BARRIER

#ifdef USE_TIMERS_2
  start = times( &start_buff );
#endif

#ifdef DEBUG_SEND_RCV
  pgslib_output_c("Calling final barrier in CSR");
#endif

  mpi_err = MPI_Barrier( MPI_COMM_WORLD );

#ifdef USE_TIMERS_2
  stop  = times( &stop_buff );
  barrier_time += (float)(( stop_buff.tms_utime +  stop_buff.tms_stime) - 
		   (start_buff.tms_utime + start_buff.tms_stime)) / tick_tock;
#endif

#ifdef DEBUG_SEND_RCV
  pgslib_output_c("Passed final barrier in CSR");
  pgslib_flush_output_c();
#endif

#else
  barrier_time += 0;
#endif  
return;
} /* end of pgslib_constrained_send_receive */

