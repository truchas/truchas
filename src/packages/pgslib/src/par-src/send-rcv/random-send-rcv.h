/* Copyright Robert C. Ferrell, CPCA Ltd.  1995 */

/* $Id: random-send-rcv.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

/*
   This is the generic routine, which is included in the type specific files.
   void pgslib_random_send_rcv(SendS, RcvR)

This routines takes two inputs, a COMM_SEND and a COMM_RCV structure.
  COMM_SEND contains the data for the messages to be send.
  COMM_RCV  holds the received data.  Depending on the particular
            need, COMM_RCV may be pre-allocated, or parts of it
	    may be allocated on the fly.

Since this is a random send, we don't know how many messages we will
receive.  We handle this by doing a global count of the number of
messages to be sent, then keeping track of the number of messages
received.  Once we've received (globally) as many messages were sent, 
were done.

*/

void PGSLIB_ROUTINE_NAME(pgslib_random_send_rcv_)(SendS, RcvR)
COMM_SEND *SendS;
COMM_RCV  *RcvR;
{ 
  int nproc, myrank;
  int N_Send, N_Send_Total, N_Rcvd, N_Rcvd_Total;
  int s, l, Rcvs_Between_Tests;
  SEND_BUFFER *this_S_Buff;
  MPI_Request *send_request;
  MPI_Status  *send_status;
  char estr[512];

  MPI_Comm_size( MPI_COMM_WORLD, &nproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank);

  N_Send = SendS->N_Send;
  if (N_Send <0) N_Send = 0;
  N_Send_Total = 0;
  N_Rcvd = 0;
  N_Rcvd_Total = 0;

/* Global count of the number of messages to be sent. */
  MPI_Allreduce(&N_Send, &N_Send_Total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#ifdef DEBUG_SEND_RCV
  sprintf(estr, "counted %d messages to send", N_Send_Total);
  pgslib_output_c(estr);
#endif

/* Do the sends, and do some attempts some receives */

  /* Only need to use send_? if we have finite number of messages to send. */
  if (N_Send > 0) {
    send_request = (MPI_Request *)malloc(N_Send*sizeof(MPI_Request));
    pgslib_check_malloc_c(send_request, "malloc of send_request failed in random_send_rcv.");

    send_status = (MPI_Status *)malloc(N_Send*sizeof(MPI_Status));
    pgslib_check_malloc_c(send_status, "malloc of send_status failed in random_send_rcv.");
  }

  for(s=0; s<N_Send; s++) {
    this_S_Buff = (SendS->S_Buffers)[s];

#ifdef DEBUG_SEND_RCV
    sprintf(estr,"Sending message %d: n_send=%d, PE=%d, TAG=%d", s,
	    this_S_Buff->n_send_words, this_S_Buff->dest_PE.PE, this_S_Buff->send_TAG);
    pgslib_output_c(estr);
#endif
    MPI_Isend(&((this_S_Buff->PGSLIB_TYPE_NAME(send_data_))[0]),
	      this_S_Buff->n_send_words,
	      PGSLIB_MPI_DATA_TYPE,
	      this_S_Buff->dest_PE.PE,
	      this_S_Buff->send_TAG,
	      MPI_COMM_WORLD,
	      &(send_request[s]) );
#ifdef DEBUG_SEND_RCV
    sprintf(estr, "sent message %d, attempting receive", s);
    pgslib_output_c(estr);
#endif
    PGSLIB_ROUTINE_NAME(pgslib_attempt_receive_)(RcvR, myrank);
  }
  
/* All the sends are complete, so now count how many receives have
   been made? */
  N_Rcvd = RcvR->N_Rcvd;
  MPI_Allreduce(&N_Rcvd,&N_Rcvd_Total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#ifdef DEBUG_SEND_RCV
  sprintf(estr,"Received %d messages, %d total", N_Rcvd, N_Rcvd_Total);
  pgslib_output_c(estr);
#endif

  Rcvs_Between_Tests = 2;
  while(N_Rcvd_Total < N_Send_Total) {
    for(l=0;l<Rcvs_Between_Tests;l++)
      PGSLIB_ROUTINE_NAME(pgslib_attempt_receive_)(RcvR, myrank);
    N_Rcvd = RcvR->N_Rcvd;
    MPI_Allreduce(&N_Rcvd,&N_Rcvd_Total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }

#ifdef DEBUG_SEND_RCV
  sprintf(estr,"After wait loop, received %d messages, %d total", N_Rcvd, N_Rcvd_Total);
  pgslib_output_c(estr);
#endif

/* Finally, finish up the send requests.  This is a formality, since
   we could not have finished the receives if the sends weren't complete.*/
  if (N_Send > 0) {
    MPI_Waitall(N_Send, &send_request[0], &send_status[0]);

/* Done, so free the memory. */
    free(send_status);
    free(send_request);
  }

return;
} /* end of pgslib_random_send_receive */

