/* Copyright Robert C. Ferrell, CPCA Ltd.  1995 */

/* $Id: attempt-receive.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

/*
   This is the generic routine, which gets included in the type specific files.

   int pgslib_attempt_receive(RcvR, myrank)
     COMM_RCV *RcvR;
     int       myrank;

This routine attemps to receive a message.  It first probes for a
message. If one has come in, then it receives it, otherwise it returns;

RETURN VALUE:
       0 = no message received
       > 0  message received

SIDE EFFECTS: 

If return value == 0, then there are no side effects.  

Data in RcvR is set by this routine.  If return value > 0, then a
  message was received.  
  N_Rcvd is incremented by one.
  N_Rcv_Buffers may have been increased, depending on flags.
*/


int PGSLIB_ROUTINE_NAME(pgslib_attempt_receive_)(RcvR, myrank)
     COMM_RCV *RcvR;
     int       myrank;
{
  int flag, tag, count, rcvBIndex;
  MPI_Status pstatus, rstatus;
  RCV_BUFFER *rcvB;
  char estr[512];

  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &pstatus);

/* If a message has come in, receive it. */
  if (flag) {   

/* What to do with the tag? */
    tag = pstatus.MPI_TAG;
 
    /* Put this message in a particular place? */
    if(RcvR->USE_TAG) {
      rcvBIndex = tag;
      rcvB = (RcvR->R_Buffers)[rcvBIndex];
      if(tag != rcvB->rcv_TAG)
	pgslib_error_c("Wrong tag in attempt receive.");
    }
    /* Or put it in the next RCV_BUFFER? */
    else {
#ifdef DEBUG_SEND_RCV
      sprintf(estr, "RcvR->N_Rcv_Buffers = %d, RcvR->N_Rcvd = %d", RcvR->N_Rcv_Buffers, RcvR->N_Rcvd);
      pgslib_output_c(estr);
#endif
      if (RcvR->N_Rcvd >= RcvR->N_Rcv_Buffers) 
	pgslib_fatal_error_c("Not enough RCV_BUFFER in pgslib_attempt_receive.");
      rcvBIndex = (RcvR->N_Rcvd);  /*This is next available slot*/
      rcvB = (RcvR->R_Buffers)[rcvBIndex];
      rcvB->rcv_TAG = tag;
    }

/* Do we know where this message should be coming from?
   (This usually only makes sense if we USE_TAG, but not always.) */
    if(RcvR->KNOWN_SRC){
      if(rcvB->src_PE.PE != pstatus.MPI_SOURCE) {
	pgslib_error_c("message from wrong src PE in pgslib_attemp_receive.");
      }
    }
    else {
      rcvB->src_PE.PE = pstatus.MPI_SOURCE;
#ifdef DEBUG_SEND_RCV
      sprintf(estr, "pstatus.MPI_SOURCE = %d", pstatus.MPI_SOURCE);
      pgslib_output_c(estr);
#endif
    }


/* Is the receive data buffer ready for receiving? */
    MPI_Get_count(&pstatus, PGSLIB_MPI_DATA_TYPE, &count);

    /* Do we allow resizing of the rcv data buffer? */
    if (RcvR->ALLOW_RCV_BUFF_RESIZE) {

      /* If so, have we allocated it already, in which case we have to free the buffer first. */
      if(RcvR->RCV_BUFF_PREALLOCATED) 
	free(rcvB->PGSLIB_TYPE_NAME(rcv_data_));
      
      rcvB->PGSLIB_TYPE_NAME(rcv_data_) = (PGSLIB_DATA_TYPE *)malloc((RcvR->Element_Size)*count);
      pgslib_check_malloc_c(rcvB->PGSLIB_TYPE_NAME(rcv_data_), "rcv_data malloc failed in pgslib_attempt_receive.");
    }	

    /* If we don't allow resizing, is the current buffer large enough? */
    else {
      if(count > rcvB->n_max_rcv_words)
	pgslib_fatal_error_c("Rcv Buffer too small in pgslib_attemp_receive.");
    }
	
/* Receive the message */

#ifdef DEBUG_SEND_RCV
    sprintf(estr, "count=%d, src_PE.PE=%d, rcv_TAG=%d", count, rcvB->src_PE.PE, rcvB->rcv_TAG);
    pgslib_output_c(estr);
#endif
    MPI_Recv(rcvB->PGSLIB_TYPE_NAME(rcv_data_), 
	     count, 
	     PGSLIB_MPI_DATA_TYPE, 
	     rcvB->src_PE.PE,
	     rcvB->rcv_TAG,
	     MPI_COMM_WORLD,
	     &rstatus);
    MPI_Get_count(&rstatus, PGSLIB_MPI_DATA_TYPE, &count);
    rcvB->n_rcvd_words = count;
    (RcvR->N_Rcvd)++;
#ifdef DEBUG_SEND_RCV
    sprintf(estr,"received %d words, first word = %d",count,(int)(rcvB->PGSLIB_TYPE_NAME(rcv_data_)[0]));
    pgslib_output_c(estr);
#endif
  }
  return(flag);
} /* end of pgslib_attempt_receive */
