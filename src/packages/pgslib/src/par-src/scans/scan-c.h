/* C routines to support PREFIX and SUFFIX operations */
/* These routines perform a global scan.
   The routines take only a single item, Src, per processor.  The item
   has two components, (data, segment).

   This routine does an cross-pe scan of the Src.

   The results of the scan are returned in Dest.

   Every PE has exactly one item.  If that PE is "empty" at a higher level of the 
   scan algorithm, code at that level passes the identity, (0, FALSE), to this routine.

   The segment component of each item is a logical, but passed as an integer.  The
   allowed values are: 1 = TRUE
		       0 = FALSE

   The combiner is slightly different for upward (PREFIX) and downward (SUFFIX) scans.

   We are using start bits for segments.  The start bit is a start-when-encountered-bit.  That
   is, for upward scans (PREFIX) the start bit is set at the lowest index of each segment.  For
   downward scans (SUFFIX) the start bit is set at the highest index of each segment.

   Combiner for upward (PREFIX)
        (u,i) * (v,j) = (v + u(!j), i || j)
  
   Combiner for downward (SUFFIX)
        (u,i) * (v,j) = (u + v(!i), i || j)
  
   These combiners can be combined:)  Set (P,S) = (1,0) for PREFIX and (0,1) for SUFFIX scans.
   Then the universal combiner is:
        (u,i) * (v,j) = (u[P*(!j) + S] + v[P + S*(!i)], i || j)


*/

/* $Id: scan-c.h,v 1.1.1.1 2000/10/11 22:44:23 ferrell Exp $ */

#include "tags-c.h"

#ifndef PGSLIB_MPI_DATA_TYPE
#error "PGSLIB_MPI_DATA_TYPE must be defined before scan-c.h is included"
#endif

#ifndef PGSLIB_DATA_TYPE
#error "PGSLIB_DATA_TYPE must be defined before scan-c.h is included"
#endif

#ifndef _SRC_SIGN_
#error "_SRC_SIGN_ must be defined before scan-c.h is included."
#endif

#ifndef _DEST_SIGN_
#error "_DEST_SIGN_ must be defined before scan-c.h is included."
#endif

#ifdef _PREFIX_
#ifdef _SUFFIX_
#error "IF _PREFIX_ is defined, then _SUFFIX_ must not be."
#endif
#endif

#ifdef _SUFFIX_
#ifdef _PREFIX_
#error "IF _SUFFIX_ is defined, then _PREFIX_ must not be."
#endif
#endif

#ifndef _SCAN_ITEM_
typedef struct Scan_Item_t {
    PGSLIB_DATA_TYPE  data;
    int               seg;
  } Scan_Item;
#define _SCAN_ITEM_
#endif


#ifndef _COMBINER_
#define _COMBINER_
static Scan_Item Combiner(Scan_Item Left, Scan_Item Right, int Upward, int Downward)
{
  Scan_Item Result;

  Result.data = Left.data  * (Upward *(!Right.seg) + Downward) +
                Right.data * (Upward + Downward * (!Left.seg));
  Result.seg  = Left.seg || Right.seg;
  return(Result);
}
#endif

void PGSLIB_ROUTINE_NAME ( PGSLIB_DATA_TYPE *Dest_Data, int *Dest_Seg, PGSLIB_DATA_TYPE *Src_Data, int *Src_Seg)
{

  Scan_Item Dest_Item, Send_Item, Rcv_Item, ID_Item;

  int         srcIndex, destIndex, count, tag, base_tag;
  int         log2nPE, nPE, thisPE;
  unsigned    n;
  MPI_Request send_request;
  MPI_Status  s_status, r_status;
  char        estring[1024];

  MPI_Comm_size( MPI_COMM_WORLD, &nPE);
  MPI_Comm_rank( MPI_COMM_WORLD, &thisPE);

  
#ifdef USE_SCAN_BARRIER
#ifdef DEBUG_SCAN
     sprintf(estring,"Calling barrier in scan");
     pgslib_output_c(estring);       
#endif
       
  MPI_Barrier( MPI_COMM_WORLD );
#endif

  n       = nPE>>1;
  log2nPE = 0;
  
  for(n=nPE>>1;n;n=n>>1)
    {	
#ifdef DEBUG_SCAN
      sprintf(estring,"In log loop, n = %d\n", n);
	pgslib_output_c(estring);	
#endif
	log2nPE++;
      }

#ifdef DEBUG_SCAN
  sprintf(estring, "log2nPE = %d \n", log2nPE);
  pgslib_output_c(estring);          
#endif

  /* Frequently we need the ID */
  ID_Item.data = (PGSLIB_DATA_TYPE)0;
  ID_Item.seg  = FALSE;

  /* At start, Rcv_Item is empty and Dest_Item is Src. */
  
  Dest_Item.data = *Src_Data;
  Dest_Item.seg  = *Src_Seg;

  Rcv_Item       = ID_Item;
  
#ifdef DEBUG_SCAN
  sprintf(estring,"Starting scan loop.\n");
  pgslib_output_c(estring);          
#endif

  n = 0;

  srcIndex  = thisPE _SRC_SIGN_  (1<<n); /* (1<<n) = 2**n */
  destIndex = thisPE _DEST_SIGN_ (1<<n);

#ifdef DEBUG_SCAN
  sprintf(estring,"srcIndex = %d, destIndex = %d \n", srcIndex, destIndex);
  pgslib_output_c(estring);       
#endif

/* Loop until all  messages sent and received */

/* The tag is unique for each call, so don't need a barrier. */

  base_tag = scan_TagBits();

  while(  ( (0 <= srcIndex)  && ( srcIndex < nPE) )
	||( (0 <= destIndex) && (destIndex < nPE) ))
    {
#ifdef DEBUG_SCAN
      sprintf (estring, "Top of scan loop, n = %d, log2nPE = %d, \n",n, log2nPE);
      pgslib_output_c(estring);      
      { int largest_tag, tag_flag;
        MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &largest_tag, &tag_flag);
	sprintf(estring,"Largest allowed tag = %d\n", largest_tag);
	pgslib_output_c(estring);
      }
#endif

      /* Sent item is always current Dest */
      Send_Item = Dest_Item;
      /* Some PE's don't receive anything.  In order for the combiner
	 to make sense, give those the ID.  Then post-communicate combiner can 
	 be used on all PEs'. */
      Rcv_Item  = ID_Item;

      tag = base_tag + n;
#ifdef DEBUG_SCAN
     sprintf(estring,"Before ISEND, srcIndex = %d, destIndex = %d, tag = %d\n", srcIndex, destIndex, tag);
     pgslib_output_c(estring);       
#endif

     count = 1;
     if ((0 <= destIndex) && (destIndex < nPE))
       MPI_Isend(&Send_Item,
		 count,
		 PGSLIB_MPI_DATA_TYPE,
		 destIndex,
		 tag,
		 MPI_COMM_WORLD,
		 &send_request);
     if ((0 <= srcIndex) && (srcIndex < nPE))
       MPI_Recv(&Rcv_Item,
		count,
		PGSLIB_MPI_DATA_TYPE,
		srcIndex,
		tag,
		MPI_COMM_WORLD,
		&r_status);
     if ((0 <= destIndex) && (destIndex < nPE))
       MPI_Waitall(1, &send_request, &s_status);
     
     /* Combine Rcv with Dest to make a new Dest */
#ifdef _PREFIX_
     Dest_Item = Combiner(Rcv_Item, Dest_Item, 1, 0);
#endif
#ifdef _SUFFIX_
     Dest_Item = Combiner(Dest_Item, Rcv_Item, 0, 1);     
#endif

     n++;
     srcIndex  = thisPE _SRC_SIGN_  (1<<n); /* (1<<n) = 2**n */
     destIndex = thisPE _DEST_SIGN_ (1<<n);

#ifdef DEBUG_SCAN
     sprintf(estring,"Preparing for next loop, srcIndex = %d, destIndex = %d \n", srcIndex, destIndex);
      pgslib_output_c(estring);       
#endif
   }
       
#ifdef USE_SCAN_BARRIER
  MPI_Barrier( MPI_COMM_WORLD );
#endif
				/* Each PE now has donor, which it must supply
				 to next PE. */
  srcIndex  = thisPE _SRC_SIGN_  1;
  destIndex = thisPE _DEST_SIGN_ 1;

  Send_Item = Dest_Item;
  /* Some PE's don't receive anything.  In order to return a useful result,
     initialize Rcv_Item to ID.  */
  Rcv_Item  = ID_Item;
  
  count     = 1;
  tag       = base_tag + nPE + 1;
  
#ifdef DEBUG_SCAN
  sprintf(estring,"Final shift, srcIndex = %d, destIndex = %d, tag = %d \n", srcIndex, destIndex, tag);
  pgslib_output_c(estring);        
#endif

  if ((0 <= destIndex) && (destIndex < nPE))
    MPI_Isend(&Send_Item,
	      count,
	      PGSLIB_MPI_DATA_TYPE,
	      destIndex,
	      tag,
	      MPI_COMM_WORLD,
	      &send_request);
  if ((0 <= srcIndex) && (srcIndex < nPE))
    MPI_Recv(&Rcv_Item,
	     count,
	     PGSLIB_MPI_DATA_TYPE,
	     srcIndex,
	     tag,
	     MPI_COMM_WORLD,
	     &r_status);
  if ((0 <= destIndex) && (destIndex < nPE))
    MPI_Waitall(1, &send_request, &s_status);
  
  Dest_Item = Rcv_Item;
#ifdef DEBUG_SCAN
  sprintf(estring,"Dest_Item.data = %d\n",(int)(Dest_Item.data));
  pgslib_output_c(estring);         
#endif

  *Dest_Data = Dest_Item.data;
  *Dest_Seg  = Dest_Item.seg;

#ifdef USE_SCAN_BARRIER
  MPI_Barrier( MPI_COMM_WORLD );
#endif

  return;
}
#undef PGSLIB_MPI_DATA_TYPE
#undef PGSLIB_DATA_TYPE
#undef PGSLIB_ROUTINE_NAME
#undef _SRC_SIGN_
#undef _DEST_SIGN_
