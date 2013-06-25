/* Support routines for gather scatter */
/* This file contains the  setup routines. */

/* $Id: gs-setup-c.c,v 1.1.1.1 2000/10/11 22:44:25 ferrell Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "pgslib-include-c.h"
#define PGSLIB_DATA_TYPE int
#include "pgslib-types.h"
#include "gs-c.h"

void pgslib_setup_n_duplicate_c(N_Duplicate_ptr, gsTrace_ptrptr)
     int *N_Duplicate_ptr;
     GS_TRACE_STRUCT **gsTrace_ptrptr;
{
  COMM_SEND *SendS;
  COMM_RCV  *RcvR;
  SEND_BUFFER *SBuff;
  RCV_BUFFER *RBuff;
  CmplOwnerPE *CmplOwnPE;
  CmplReceiverPE *CmplRcvrPE;
  GS_TRACE_STRUCT  *gsTrace_ptr;
  int i, cmpl, npe, N_Rcvd, ncmplCommBuffer, Offset;
  int *rcv_data, *send_data, *global_cmplCommBuffer;
  int N_Duplicate;
  char estr[512];

  gsTrace_ptr = *gsTrace_ptrptr;
/* Send to each of N_CmplOwnerPEs the number of cmpls needed. */

/* Need Send and RcvR structures */
  SendS = (COMM_SEND *) malloc(sizeof(COMM_SEND));
  pgslib_check_malloc_c(SendS, "malloc failed for SendS in gs_setup_n_duplicate");

  RcvR = (COMM_RCV *) malloc(sizeof(COMM_RCV));
  pgslib_check_malloc_c(RcvR, "malloc failed for RcvR in gs_setup_n_duplicate");

/* Fill in SendS, which will be passed into pgslib_random_send_rcv */
  SendS->N_Send = gsTrace_ptr->N_CmplOwnerPEs;
  

/* S_Buffers is accessed through double indirection, so that
   we can change the order of the sends if desired. */

#ifdef DEBUG_GATH_SCATT
  sprintf(estr, "SendS->N_Send = %d", SendS->N_Send);
  pgslib_output_c(estr);
#endif

  if (SendS->N_Send > 0) {
    SendS->S_Buffers = (SEND_BUFFER **)malloc(sizeof(SEND_BUFFER *)*SendS->N_Send);
    pgslib_check_malloc_c(SendS->S_Buffers, "malloc failed in gs_setup_n_duplicate");

    SendS->S_Buffer_Data = (SEND_BUFFER *)malloc(sizeof(SEND_BUFFER)*SendS->N_Send);
    pgslib_check_malloc_c((SendS->S_Buffer_Data), "malloc of S_Buffer_Data failed in gs_setup_n_duplicate");

    SendS->S_Order = (int *)malloc(sizeof(int)*SendS->N_Send);
    pgslib_check_malloc_c(SendS->S_Order, "malloc of SendS->S_Order failed in gs_setup_n_duplicate");

    SendS->send_request = (MPI_Request *)malloc(SendS->N_Send*sizeof(MPI_Request));
    pgslib_check_malloc_c(SendS->send_request, "malloc of SendS->send_request failed in gs_setup_n_duplicate.");
    
    SendS->send_status = (MPI_Status *)malloc(SendS->N_Send*sizeof(MPI_Status));
    pgslib_check_malloc_c(SendS->send_status, "malloc of SendS->send_status failed in gs_setup_n_duplicate.");
  }
  else  {
    SendS->S_Buffers     = NULL;
    SendS->S_Buffer_Data = NULL;
    SendS->S_Order       = NULL;
    SendS->send_request  = NULL;
    SendS->send_status   = NULL;
  }

/* For now, just walk through the send buffer in linear order. */
  for(i=0; i< SendS->N_Send; i++)
    (SendS->S_Order)[i] = i;
  
  for(i=0; i< SendS->N_Send; i++)
    SendS->S_Buffers[SendS->S_Order[i]] = &(SendS->S_Buffer_Data[i]);

/* Move data into SendS->  Only one word will be sent, = 
   number of cmpls needed.  The tag sent is the offset
   to be used for the reply.*/

  for(i=0;i<SendS->N_Send;i++) {
    SBuff = &(SendS->S_Buffer_Data[i]);
    CmplOwnPE = gsTrace_ptr->CmplOwnerPEs + i;
    SBuff->dest_PE.PE = CmplOwnPE->PE.PE;
    SBuff->send_TAG = i;
    SBuff->n_send_words = 1;
    SBuff->send_data_int = &(CmplOwnPE->ncmpls);
  }


/* Set up receive buffer for incoming messaes.
   The incoming data will be stored in the order it is received.
   At most one word from each PE will be received, so we
   have a firm upper bound on the size of the buffers we need.*/
  
  npe = pgslib_state.nPE;

  RcvR->R_Buffers = (RCV_BUFFER **)malloc(sizeof(RCV_BUFFER *)*npe);
  pgslib_check_malloc_c(RcvR->R_Buffers, "could not malloc space for R_Buffers * in pgslib_init_comm");

  (RcvR->R_Buffers)[0] = (RCV_BUFFER *)malloc(sizeof(RCV_BUFFER)*npe);
  pgslib_check_malloc_c((RcvR->R_Buffers)[0], "could not malloc space for R_Buffers in pgslib_init_comm");

  /* Since we are using random send-recv, don't use IRecv */
  RcvR -> rcv_request = NULL;
  RcvR -> rcv_status  = NULL;
     
/* Fill in some flags in RcvR */
  RcvR->N_Rcv_Buffers = npe;
  RcvR->N_Rcvd = 0;
  RcvR->USE_TAG = 0;   /* Since messages will be stored in order received */
  RcvR->KNOWN_SRC = 0; /* Since don't know who will be asking us for cmpls */
  RcvR->ALLOW_RCV_BUFF_RESIZE = 0; /* Since we have upper bound on amount of incoming data */
  RcvR->RCV_BUFF_PREALLOCATED = 1; /* But this doesn't matter, because of line above. */
  RcvR->Element_Size = sizeof(int); /* Since that is what is being sent. */

  rcv_data = (int *)malloc(sizeof(int)*npe);
  pgslib_check_malloc_c(rcv_data, "malloc failed for rcv_data in pgslib_init_comm.");

/* Now fill in some fields in R_Buffers */
  for(i=0;i<npe;i++) {
    (RcvR->R_Buffers)[i] = (RcvR->R_Buffers)[0] + i;
    RBuff = (RcvR->R_Buffers)[i];
    RBuff -> n_max_rcv_words = 1;
    RBuff -> rcv_data_int = rcv_data+i;
  }


/* Now do the exchange. */
#ifdef DEBUG_GATH_SCATT
  pgslib_output_c("Calling pgslib_random_send_rcv");
#endif
  pgslib_random_send_rcv_int_c(SendS, RcvR);

/* No barrier is necessary, because random_send doesn't complete until all
   messages have been sent and received. */

/* Unpack the received data, if any */
  N_Rcvd = RcvR->N_Rcvd;
  gsTrace_ptr->N_CmplReceiverPEs = RcvR->N_Rcvd;

  if (N_Rcvd > 0) {
/* Make space in gsTrace_ptr for the CmplReceiverPE info*/
    gsTrace_ptr->CmplReceiverPEs = (CmplReceiverPE *)malloc(sizeof(CmplReceiverPE)*N_Rcvd);
    pgslib_check_malloc_c(gsTrace_ptr->CmplReceiverPEs, "failed to malloc CmplReceiverPEs in pgslib_setup_n_duplicate.");
    
#ifdef DEBUG_GATH_SCATT
    sprintf(estr,"gsTrace_ptr->N_CmplReceiverPEs=%d",gsTrace_ptr->N_CmplReceiverPEs);
    pgslib_output_c(estr);
#endif
    
    ncmplCommBuffer = 0;
    Offset = 0;
    N_Duplicate = 0;
    for(i=0;i<N_Rcvd;i++) {
      RBuff = (RcvR->R_Buffers)[i];
      
#ifdef DEBUG_GATH_SCATT
      sprintf(estr,"(RcvR->R_Buffers)[%d]->  rcv_data[0]=%d, PE=%d, TAG=%d",
	      i, (RcvR->R_Buffers)[i]->rcv_data_int[0], (RcvR->R_Buffers)[i]->src_PE.PE, 
	      (RcvR->R_Buffers)[i]->rcv_TAG );
      pgslib_output_c(estr);
#endif
      
      CmplRcvrPE = (gsTrace_ptr->CmplReceiverPEs) + i;
      CmplRcvrPE->PE.PE = RBuff->src_PE.PE;
      CmplRcvrPE->Tag = RBuff->rcv_TAG;
      CmplRcvrPE->Offset = Offset;
      CmplRcvrPE->ncmpls = RBuff->rcv_data_int[0];
      N_Duplicate += CmplRcvrPE->ncmpls;
      Offset += CmplRcvrPE->ncmpls;
#ifdef DEBUG_GATH_SCATT
      sprintf(estr,"CmplRcvrPE-> PE=%d, Tag=%d, Offset=%d, ncmpls=%d",
	      CmplRcvrPE->PE.PE, CmplRcvrPE->Tag, CmplRcvrPE->Offset, CmplRcvrPE->ncmpls);
      pgslib_output_c(estr);
#endif
    }
    
  }
  else
    N_Duplicate = 0;

  gsTrace_ptr->N_Duplicate = N_Duplicate;
#ifdef DEBUG_GATH_SCATT
  sprintf(estr,"N_Duplicate=%d",gsTrace_ptr->N_Duplicate);
  pgslib_output_c(estr);
#endif
  
  *N_Duplicate_ptr = N_Duplicate;

/*  We are done with this */

  gs_free_send_buffer(SendS);
  gs_free_rcv_buffer(RcvR);
  if (rcv_data != NULL)
    free(rcv_data);	

}

/********************************************************************/
void pgslib_setup_duplicate_buffer_c(gsTrace_ptrptr)
     GS_TRACE_STRUCT **gsTrace_ptrptr;
{
  COMM_SEND *SendS;
  COMM_RCV  *RcvR;
  SEND_BUFFER *SBuff;
  RCV_BUFFER *RBuff;
  CmplOwnerPE *CmplOwnPE;
  CmplReceiverPE *CmplRcvrPE;
  GS_TRACE_STRUCT  *gsTrace_ptr;
  int i, cmpl, npe, N_Rcvd, ncmplCommBuffer, Offset;
  int *rcv_data, *send_data, *global_cmplCommBuffer;
  int N_Duplicate;
  char estr[512];

  gsTrace_ptr = *gsTrace_ptrptr;
/* We are now ready to reply with a tag.  The tag is the index to be used 
   for lookup into gsENdataptr->CmplReceiverPEs.  We are using a simple
   order-received right now, but probably will want to change that. */


/* Need Send and RcvR structures */
  SendS = (COMM_SEND *) malloc(sizeof(COMM_SEND));
  pgslib_check_malloc_c(SendS, "malloc failed for SendS in gs_setup_n_duplicate");

  RcvR = (COMM_RCV *) malloc(sizeof(COMM_RCV));
  pgslib_check_malloc_c(RcvR, "malloc failed for RcvR in gs_setup_n_duplicate");


/* In this case, the number of messages to be sent is N_CmplReceiverPEs, and the number
   of messages to receive is N_CmplOwnerPEs*/
  SendS->N_Send = gsTrace_ptr->N_CmplReceiverPEs;

  if (SendS->N_Send > 0) {
    SendS->S_Buffers = (SEND_BUFFER **)malloc(sizeof(SEND_BUFFER *)*SendS->N_Send);
    pgslib_check_malloc_c(SendS->S_Buffers, "couldn't malloc SendS->S_Buffers in pgslib_exchang_indices_c");
    
    SendS->S_Buffer_Data = (SEND_BUFFER *)malloc(sizeof(SEND_BUFFER)*SendS->N_Send);
    pgslib_check_malloc_c((SendS->S_Buffer_Data), "malloc of S_Buffer_Data failed in gather");

    SendS->S_Order = (int *)malloc(sizeof(int)*SendS->N_Send);
    pgslib_check_malloc_c(SendS->S_Order, "malloc of SendS->S_Order failed in gather");

    send_data = (int *)malloc(sizeof(int)*SendS->N_Send);
    pgslib_check_malloc_c(send_data, "failed to malloc send_data in pgslib_exchang_indices_c");

    SendS->send_request = (MPI_Request *)malloc(SendS->N_Send*sizeof(MPI_Request));
    pgslib_check_malloc_c(SendS->send_request, "malloc of SendS->send_request failed in pgslib_exchang_indices_c.");
    
    SendS->send_status = (MPI_Status *)malloc(SendS->N_Send*sizeof(MPI_Status));
    pgslib_check_malloc_c(SendS->send_status, "malloc of SendS->send_status failed in pgslib_exchang_indices_c.");
  }
  else{
    SendS->S_Buffers     = NULL;
    SendS->S_Buffer_Data = NULL;
    SendS->S_Order       = NULL;
    SendS->send_request  = NULL;
    SendS->send_status   = NULL;
    send_data            = NULL;
  }
  
/* We might want to reorder the sends, but for now, just use linear order */
  for(i=0; i< SendS->N_Send; i++)
    (SendS->S_Order)[i] = i;
  
  for(i=0; i< SendS->N_Send; i++)
    SendS->S_Buffers[SendS->S_Order[i]] = &(SendS->S_Buffer_Data[i]);

  for(i=0;i<SendS->N_Send;i++){
    SBuff = &(SendS->S_Buffer_Data[i]);
    CmplRcvrPE = (gsTrace_ptr->CmplReceiverPEs) + i;
    SBuff->dest_PE.PE = CmplRcvrPE->PE.PE;
    SBuff->send_TAG = CmplRcvrPE->Tag;
    SBuff->n_send_words = 1;
    send_data[i] = i;
    SBuff->send_data_int = send_data + i;
  }

/* For receiving, we will receive exactly N_CmplOwnerPEs messages */
  RcvR->N_Rcv_Buffers = gsTrace_ptr->N_CmplOwnerPEs;
  RcvR->N_Rcvd = 0;

  if (RcvR->N_Rcv_Buffers > 0) {
    RcvR->R_Buffers = (RCV_BUFFER **)malloc(sizeof(RCV_BUFFER *)*RcvR->N_Rcv_Buffers);
    pgslib_check_malloc_c(RcvR->R_Buffers, "failed to malloc Rcv.R_Buffers in pgslib_exchang_indices_c.");
    (RcvR->R_Buffers)[0] = (RCV_BUFFER *)malloc(sizeof(RCV_BUFFER)*RcvR->N_Rcv_Buffers);
    pgslib_check_malloc_c((RcvR->R_Buffers)[0], "failed to malloc Rcv.R_Buffers[0] in pgslib_exchang_indices_c.");

    rcv_data = (int *)malloc(sizeof(int)*RcvR->N_Rcv_Buffers);
    pgslib_check_malloc_c(rcv_data, "malloc failed for rcv_data in pgslib_exchang_indices_c.");

    RcvR->rcv_request = (MPI_Request *)malloc(RcvR->N_Rcv_Buffers*sizeof(MPI_Request));
    pgslib_check_malloc_c(RcvR->rcv_request, "malloc of RcvR->rcv_request failed in pgslib_exchang_indices_c.");
    
    RcvR->rcv_status = (MPI_Status *)malloc(RcvR->N_Rcv_Buffers*sizeof(MPI_Status));
    pgslib_check_malloc_c(RcvR->rcv_status, "malloc of RcvR->rcv_status failed in pgslib_exchang_indices_c.");
  }
  else {
    RcvR->R_Buffers   = NULL;
    RcvR->rcv_request = NULL;
    RcvR->rcv_status  = NULL;
    rcv_data          = NULL;
  }

  RcvR->USE_TAG = 1;    /* message tag is index into RcvR->R_Buffers */
  RcvR->KNOWN_SRC = 1;   /* we know where it is coming from */
  RcvR->ALLOW_RCV_BUFF_RESIZE = 0;  /* only receiving one word */
  RcvR->RCV_BUFF_PREALLOCATED = 1;   /* but this doesn't matter because of line above */
  RcvR->Element_Size = sizeof(int);


/* initialize receive buffers */
  for(i=0;i<RcvR->N_Rcv_Buffers;i++) {
    (RcvR->R_Buffers)[i] = (RcvR->R_Buffers)[0] + i;
    RBuff = RcvR->R_Buffers[i];
    CmplOwnPE = (gsTrace_ptr->CmplOwnerPEs) + i;
    RBuff->src_PE.PE = CmplOwnPE->PE.PE;
    RBuff->n_max_rcv_words = 1;
    RBuff->n_rcvd_words = 0;
    RBuff->rcv_TAG = i;
    RBuff->rcv_data_int = rcv_data + i;
  }

/* Exchange data (which is the tag that CmplOwnerPE uses to identify
   itself to CmplReceiverPE*/
#ifdef DEBUG_GATH_SCATT
  pgslib_output_c("calling pgslib_cnstd_send_rcv");
#endif
  pgslib_cnstd_send_rcv_int_c(SendS, RcvR);

/* Unpack the received data. The tag tells which slot in CmplOwnerPE
   gets which data.  However, on receipt we used the tag to put
   the data in the right slot, so we can just loop linearly here.*/

  for(i=0;i<gsTrace_ptr->N_CmplOwnerPEs;i++) {
    RBuff = RcvR->R_Buffers[i];
    CmplOwnPE = (gsTrace_ptr->CmplOwnerPEs) + i;
    CmplOwnPE->Tag = *(RBuff->rcv_data_int);
  }

/* We are done with SendS->SBuffers and RcvR->RBuffers, so free them up.
   They will be reallocated below.*/

  gs_free_send_buffer(SendS);
  gs_free_rcv_buffer(RcvR);
  if (rcv_data != NULL)
    free(rcv_data);	
  if (send_data != NULL)
    free(send_data);
/* Finally we are ready to exchange lists of needed cmpls */
  
  return;
}

void pgslib_prep_supplement_c(N_SupplementPtr, pes, gsTrace_ptrptr)
     int *N_SupplementPtr;
     int  pes[];
     GS_TRACE_STRUCT **gsTrace_ptrptr;
     /* Given an array pes of size N_Supplement, setup the
	fiedls in the trace.  The fields that get setup in 
	this routine are:
	N_CmplOnerPEs = total number of messages that will be sent/rcvd
	CmplOwnerPEOrder = indirection array of size npe for send/rcv order
	CmplOwnerPEs = array of size N_CmplOnerPEs, info about the messages
	               1. dest pe
		       2. offset of data in data buffer
     */

{	
  int pe, peIndex, offset, ncmpls, npe;
  int N_Supplement, isup, N_CmplOwnerPEs;
  char estr[512];
  GS_TRACE_STRUCT *gsTrace_ptr;

  N_Supplement = *N_SupplementPtr;
  gsTrace_ptr  = *gsTrace_ptrptr;
  npe = pgslib_state.nPE;

  gsTrace_ptr -> N_Supplement = N_Supplement;
  /* Count total number of messages to be sent.  Messages are delimited
     by change in pe number in pes array.
  */
  pe = -1;
  N_CmplOwnerPEs = 0;
  for(isup=0; isup < N_Supplement; isup++) {
    if(pe != pes[isup]) {
      N_CmplOwnerPEs += 1;
      pe = pes[isup];
    }
  }
  gsTrace_ptr -> N_CmplOwnerPEs = N_CmplOwnerPEs;

#ifdef DEBUG_GATH_SCATT
  sprintf(estr, " N_CmplOwnerPEs = %d", N_CmplOwnerPEs);
#endif

  /* Allocate up the PEOrder buffer and the CmplOwnerPEs buffer */
  gsTrace_ptr -> CmplOwnerPEOrder = (int *)malloc(sizeof(int)* npe);
  pgslib_check_malloc_c(gsTrace_ptr -> CmplOwnerPEOrder,
			 "in pgslib_setup_supplement_c, malloc for CmplOwnerPEOrder failed.");

  if (N_CmplOwnerPEs > 0) {
    gsTrace_ptr -> CmplOwnerPEs     = (CmplOwnerPE *)malloc(sizeof(CmplOwnerPE)* N_CmplOwnerPEs);
    pgslib_check_malloc_c(gsTrace_ptr -> CmplOwnerPEs,
			  "in pgslib_setup_supplement_c, malloc for CmplOwnerPEs failed.");
  }
  else	
    gsTrace_ptr -> CmplOwnerPEs = NULL;


  /* Finally, put data into proper message data fields.  This requires another loop
     over the pes buffer.
  */
  offset  = 0;
  peIndex = 0;
  ncmpls  = 0;
  pe      = -1;
  for(isup=0; isup < N_Supplement; isup++) {
    if (pe != pes[isup]) {
      offset += ncmpls;
      ncmpls  = 0;
      pe      = pes[isup];
      gsTrace_ptr -> CmplOwnerPEOrder[pe]         = peIndex;
      gsTrace_ptr -> CmplOwnerPEs[peIndex].PE.PE  = pe;
      gsTrace_ptr -> CmplOwnerPEs[peIndex].Offset = offset;
      peIndex++;
    }
    ncmpls++;
    gsTrace_ptr -> CmplOwnerPEs[peIndex-1].ncmpls = ncmpls;
  }

  return;
}

    
