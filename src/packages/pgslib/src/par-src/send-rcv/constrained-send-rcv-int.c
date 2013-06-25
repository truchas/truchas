 /*
   void pgslib_constrained_send_rcv_int_c(SendS, RcvR)

This routines takes two inputs, a COMM_SEND and a COMM_RCV structure.
  COMM_SEND contains the data for the messages to be send.
  COMM_RCV  holds the received data.  Depending on the particular
            need, COMM_RCV may be pre-allocated, or parts of it
	    may be allocated on the fly.

A barrier at the start of the routine insures that all PEs are ready for this
communication step.

*/

/* $Id: constrained-send-rcv-int.c,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "pgslib-include-c.h"
#include "timing-c.h"

#define PGSLIB_DATA_TYPE int
#define PGSLIB_ROUTINE_TYPE_POSTFIX int
#define PGSLIB_ROUTINE_NAME(Base_Name) Base_Name ## int_c
#define PGSLIB_TYPE_NAME(Base_Name) Base_Name ## int
#define PGSLIB_MPI_DATA_TYPE MPI_INT

#include "pgslib-types.h"

#include "constrained-send-rcv.h"
