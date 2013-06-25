/* Support routines for gather scatter */
/* This file contains the setup and utility routines. */

/* $Id: gs-util-c.c,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "pgslib-include-c.h"
#include "pgslib-types.h"
#include "gs-c.h"

void pgslib_gs_init_trace_c(gsTrace_ptrptr)
     GS_TRACE_STRUCT  **gsTrace_ptrptr;
{
  GS_TRACE_STRUCT *gsTrace_ptr;

  /* Make space for gsTrace */
  gsTrace_ptr = (GS_TRACE_STRUCT *)malloc(sizeof(GS_TRACE_STRUCT));
  pgslib_check_malloc_c(gsTrace_ptr, "malloc failed in pgslib_gs_init_trace");
  
  *gsTrace_ptrptr = gsTrace_ptr;

  /* Initialize the fields of gsTrace */
  gsTrace_ptr -> N_CmplOwnerPEs      = -1;
  gsTrace_ptr -> N_CmplReceiverPEs   = -1;
  gsTrace_ptr -> CmplOwnerPEOrder    = NULL;
  gsTrace_ptr -> CmplReceiverPEOrder = NULL;
  gsTrace_ptr -> CmplOwnerPEs        = NULL;
  gsTrace_ptr -> CmplReceiverPEs     = NULL;
  gsTrace_ptr -> N_Supplement        = -1;
  gsTrace_ptr -> Supplements         = NULL;
  gsTrace_ptr -> SupplementsPEs      = NULL;
  gsTrace_ptr -> N_Duplicate         = -1;
  gsTrace_ptr -> Duplicates          = NULL;
  gsTrace_ptr -> GatherSBuf          = NULL;
  gsTrace_ptr -> ScatterSBuf         = NULL;
  gsTrace_ptr -> GatherRBuf          = NULL;
  gsTrace_ptr -> ScatterRBuf         = NULL;
  return;
}

void pgslib_gs_release_trace_c(GS_TRACE_STRUCT **gsTrace_ptrptr)
{
  COMM_SEND *SendS;
  COMM_RCV  *RcvR;
  char estr[512];

  GS_TRACE_STRUCT *gsTrace_ptr;
 
   gsTrace_ptr = *gsTrace_ptrptr;
   /* Deallocate fields used for gather */
#ifdef DEBUG_GATH_SCATT
    sprintf(estr, "Releasing gsTrace_ptr->GatherSBuf");
    pgslib_output_c(estr);
    pgslib_flush_output_c();
#endif

   SendS = gsTrace_ptr ->GatherSBuf;	
   gs_free_send_buffer(SendS);
   gsTrace_ptr -> GatherSBuf = NULL;

#ifdef DEBUG_GATH_SCATT
    sprintf(estr, "Releasing gsTrace_ptr->GatherRBuf");
    pgslib_output_c(estr);
    pgslib_flush_output_c();
#endif

   RcvR = gsTrace_ptr -> GatherRBuf;
   gs_free_rcv_buffer(RcvR);
   gsTrace_ptr -> GatherRBuf = NULL;
     
   /* Deallocate fields used for scatter */

#ifdef DEBUG_GATH_SCATT
    sprintf(estr, "Releasing gsTrace_ptr->ScatterSBuf");
    pgslib_output_c(estr);
    pgslib_flush_output_c();
#endif

   SendS = gsTrace_ptr -> ScatterSBuf;
   gs_free_send_buffer(SendS);
   gsTrace_ptr -> ScatterSBuf = NULL;

#ifdef DEBUG_GATH_SCATT
    sprintf(estr, "Releasing gsTrace_ptr->ScatterRBuf");
    pgslib_output_c(estr);
    pgslib_flush_output_c();
#endif

   RcvR = gsTrace_ptr -> ScatterRBuf;
   gs_free_rcv_buffer(RcvR);
   gsTrace_ptr -> ScatterRBuf = NULL;

   /* Deallocate general comm fields */
   if (gsTrace_ptr -> CmplReceiverPEs != NULL) 
     free(gsTrace_ptr -> CmplReceiverPEs);	
   if (gsTrace_ptr -> CmplReceiverPEOrder != NULL)
     free(gsTrace_ptr -> CmplReceiverPEOrder);
   if (gsTrace_ptr -> CmplOwnerPEs != NULL) 
     free(gsTrace_ptr -> CmplOwnerPEs);
   if (gsTrace_ptr -> CmplOwnerPEOrder != NULL)
     free(gsTrace_ptr -> CmplOwnerPEOrder);

  /* Free the Trace */

#ifdef DEBUG_GATH_SCATT
    sprintf(estr, "Releasing gsTrace_ptrptr");
    pgslib_output_c(estr);
    pgslib_flush_output_c();
#endif

  free(*gsTrace_ptrptr);
  *gsTrace_ptrptr = NULL;

#ifdef DEBUG_GATH_SCATT
    sprintf(estr, "Completed releasing gsTrace_ptr");
    pgslib_output_c(estr);
    pgslib_flush_output_c();
#endif

  return;
}

void gs_free_send_buffer(COMM_SEND *SendS)
{
  char estr[512];
#ifdef DEBUG_GATH_SCATT
    sprintf(estr, "Releasing SendS");
    pgslib_output_c(estr);
    pgslib_flush_output_c();
#endif

   if (SendS != NULL) {
     if ((SendS -> send_status) != NULL)
       free(SendS -> send_status);
     if ((SendS -> send_request) != NULL)
       free(SendS -> send_request);
     if ((SendS -> S_Buffers) != NULL) {
       free(SendS -> S_Buffer_Data);
       free(SendS -> S_Buffers);
     }	
     if ((SendS -> S_Order) != NULL)
       free(SendS -> S_Order);
     free(SendS);
   }
}

void gs_free_rcv_buffer(COMM_RCV *RcvR)
{
  char estr[512];
#ifdef DEBUG_GATH_SCATT
    sprintf(estr, "Releasing RcvR");
    pgslib_output_c(estr);
    pgslib_flush_output_c();
#endif
  
   if (RcvR != NULL) {
     if ((RcvR -> rcv_status) != NULL)
       free(RcvR -> rcv_status);
     if ((RcvR -> rcv_request) != NULL)
       free(RcvR -> rcv_request);
     if ((RcvR -> R_Buffers) != NULL) {
       free((RcvR -> R_Buffers)[0]);
       free(RcvR -> R_Buffers);
       }
     free(RcvR);
   }
}

void pgslib_trace_degree_c(int *scatter_degree, int *gather_degree, GS_TRACE_STRUCT **gsTrace_ptrptr)
{
  *scatter_degree = (*gsTrace_ptrptr)->N_CmplOwnerPEs;
  *gather_degree  = (*gsTrace_ptrptr)->N_CmplReceiverPEs;
}
