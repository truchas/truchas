/* These routines and constants support the tags used in the MPI calls
   in PGSLib.

   Each kind of functionality (constrained-send-rcv, scan, etc...) has it's
   own base tag.  These are contained in the FUNCTION_TAG_BITS.  

   The low ROUTINE_INTERNAL_TAG_BITS are available to be set in each of the functions.

   The highest ITERATION_TAG_BITS bits are for the iteration count.  

   The greatest number of iterations that can be counted is:
   maxIterationTag = 1 << (ITERATION_TAG_BITS-1)

   Note that: sizeof(int) = ROUTINE_INTERNAL_TAG_BITS + FUNCTION_TAG_BITS + ITERATION_TAG_BITS

   A tag is returned according to:

   function_iterationTag = ++function_iterationTag % maxIterationTag;
   tag = (function_iterationTag << (FUNCTION_TAG_BITS + ROUTINE_INTERNAL_TAG_BITS))
       + (function_Tag << ROUTINE_INTERNAL_TAG_BITS);

   or
   tag = ( (function_iterationTag << FUNCTION_TAG_BITS) 
          + function_Tag) << ROUTINE_INTERNAL_TAG_BITS;
       
*/

/* $Id: tags-c.c,v 1.4 2002/03/04 23:49:40 ferrell Exp $ */

#include "pgslib-include-c.h"
#include "tags-c.h"

#define FUNCTION_TAG_BITS  3   /* Allow 8 distinct functions */
#define ITERATION_TAG_BITS 7  /* Allow 128 distinct outstanding (incomplete) calls to this function.
				 Might think to make this larger, but want routines to have at least
				 13 bits for their internal use.  This is already pushing it, since
				 MPI only assures us of 15 tag bits.  Here I'm assuming we get
				 at least 23 tag bits.*/
                               /* ROUTINE_INTERNAL_TAG_BITS are the number of lower order bits available to the user. */
static int ROUTINE_INTERNAL_TAG_BITS = 0;

static int maxIterationTag = 1 << (ITERATION_TAG_BITS - 1);

/* These are the tags for the various functionalities.
   The numbering is arbitrary, but they must be distinct.
   Of course, the number cannot exceed 2**FUNCTION_TAG_BITS -1 */

static int csr_Tag          = 1;
static int scan_Tag         = 2;

/* These are the iteration counts for the various functionalities */

static int csr_iterationTag  = 0;
static int scan_iterationTag = 0;

/* This initializes the tag counters, and must be called when PGSLib is
   initialized.*/
/* Keep track whether this has been called or not. */
static int tags_been_initialized = FALSE;

int initializeTagBits()
{
  int MaxAllowedTag, MaxUsedTag, AllowedTagBits, Tag_Flag;
  int *MaxAllowedTag_ptr;
  char        estring[1024];

  if (! tags_been_initialized) {
    tags_been_initialized = TRUE;

    /* Get the largest allowed tag on this system. */
    /* Have to use MPI_COMM_WORLD, since that has the attributes for the system.*/
    MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &MaxAllowedTag_ptr, &Tag_Flag);
    MaxAllowedTag = *MaxAllowedTag_ptr;

    /* Count the number of bits in the maximum tag. */
    AllowedTagBits = 0;
    MaxUsedTag = 0;
    while (MaxUsedTag <= (MaxAllowedTag>>1)){
      AllowedTagBits++;
      MaxUsedTag = 1<<AllowedTagBits;
    }
  
    /* Now we know how many bits are available to us.  Some of these are reserved
       for function and iteration tags, the rest are available for routines to use
       "internally".*/
    ROUTINE_INTERNAL_TAG_BITS  = AllowedTagBits - ( (int)ITERATION_TAG_BITS + (int)FUNCTION_TAG_BITS );

#ifdef DEBUG_TAGS
    sprintf(estring, "MaxAllowedTag = %d, MaxUsedTag = %d, AllowedTagBits = %d, Tag_Flag = %d \n",
	    MaxAllowedTag, MaxUsedTag, AllowedTagBits, Tag_Flag);
    pgslib_output_c(estring);
    sprintf(estring, "ROUTINE_INTERNAL_TAG_BITS = %d, maxIterationTag = %d \n",
	    ROUTINE_INTERNAL_TAG_BITS, maxIterationTag);
    pgslib_output_c(estring);
#endif      

  }
  return(tags_been_initialized);
}


/* This translates an iteration count into a properly shifted tag */
int tagShifter(int iterator, int function)
{
  return ( ( (iterator << FUNCTION_TAG_BITS) + function) << ROUTINE_INTERNAL_TAG_BITS);

}

/* This returns the iteration tag bits to be used by constrained-send-rcv */
int constrainedSendRcv_TagBits()
{
  /* Increment tag, wrap if necessary */
  csr_iterationTag = ++(csr_iterationTag) % maxIterationTag;
  return( tagShifter(csr_iterationTag, csr_Tag) );
}

   
/* This returns the iteration tag bits to be used by scans */
int scan_TagBits()
{
  /* Increment tag, wrap if necessary */
  scan_iterationTag = ++(scan_iterationTag) % maxIterationTag;
  return( tagShifter(scan_iterationTag, scan_Tag) );
}

   
   
   
   
