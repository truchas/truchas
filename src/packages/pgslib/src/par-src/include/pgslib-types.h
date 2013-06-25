/* These are the structs needed for the communication buffers */

/* $Id: pgslib-types.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

/* The terminology:
       Requestor:    Requests the communication.  Typically the elements are the requestors.
       Complior :    Complies with the request.  Typically the nodes are the compliors.
   Abbreviations:
       Req = Requestor
       Cmpl = Complior
*/

#ifndef PGSLIB_TYPES_H__
#define PGSLIB_TYPES_H__

#include "hash.h"
#include "send-rcv.h"

typedef int CMPL_NUM;
typedef int CMPL_TYPE;

/* General information struct about state of the system */
typedef struct
{
  int   initialized;
  int   nPE;
  int 	thisPE;
  int   io_pe;
  MPI_Comm PGSLib_Comm; 
} parallel_state;

typedef struct {
  int PE;
} ADDRESS;


typedef struct {
    ADDRESS PE;
    int     Tag;
    int     Offset;
    int     ncmpls;
} CmplOwnerPE;

typedef struct {
    ADDRESS PE;
    int     Tag;
    int     Offset;
    int     ncmpls;
} CmplReceiverPE;

typedef struct {
    int            N_CmplOwnerPEs;
    int            N_CmplReceiverPEs;
    int            *CmplOwnerPEOrder;
    int            *CmplReceiverPEOrder;
    CmplOwnerPE    *CmplOwnerPEs;
    CmplReceiverPE *CmplReceiverPEs;
    int            N_Supplement;
    CMPL_NUM       *Supplements;
    ADDRESS        *SupplementsPEs;
    int            N_Duplicate;
    int            *Duplicates;
/*    HASH_TABLE_ENTRY *Hash_Table; */
    COMM_SEND	   *GatherSBuf, *ScatterSBuf;
    COMM_RCV	   *GatherRBuf, *ScatterRBuf;
} GS_TRACE_STRUCT;  /*Gather Scatter - Trace - Structure */

typedef struct {
  HASH_TABLE_ENTRY *Hash_Table;
} GS_TABLE_STRUCT; /* Hash table container */

#endif
