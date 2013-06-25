/* C routines to perform a physical circular shift up or down */
/* These routines move data a single PE up or a single PE down.
   The number of data items must be the same on all PE's.
   The routines take 3 arguments, Dest, Src and Count.

*/

/* $Id: shift-c.h,v 1.1.1.1 2000/10/11 22:44:23 ferrell Exp $*/

#ifndef PGSLIB_MPI_DATA_TYPE
#error "PGSLIB_MPI_DATA_TYPE must be defined before scan-c.h is included"
#endif

#ifndef PGSLIB_DATA_TYPE
#error "PGSLIB_DATA_TYPE must be defined before scan-c.h is included"
#endif

#ifndef _DIRECTION_SIGN_
#error "_DIRECTION_SIGN_ must be defined before scan-c.h is included."
#endif

#ifndef SIZE_MULTIPLIER
#error "SIZE_MULTIPLIER must be defined before scan-c.h is included."
#endif


void PGSLIB_ROUTINE_NAME ( PGSLIB_DATA_TYPE *Dest_Data, PGSLIB_DATA_TYPE *Src_Data, int *Count)
{


  int         srcIndex, destIndex, count, tag;
  int         nPE, thisPE;
  MPI_Request send_request;
  MPI_Status  s_status, r_status;
  char        estring[1024];

  MPI_Comm_size( MPI_COMM_WORLD, &nPE);
  MPI_Comm_rank( MPI_COMM_WORLD, &thisPE);

  
#ifdef DEBUG_SHIFT
     sprintf(estring,"Calling barrier in shift");
     pgslib_output_c(estring);       
#endif
       
  MPI_Barrier( MPI_COMM_WORLD );

  srcIndex  = thisPE _DIRECTION_SIGN_ -1;
  destIndex = thisPE _DIRECTION_SIGN_  1;
  /* Circular shift, so loop the pe numbers */
  if (srcIndex == -1  ) srcIndex  = nPE - 1;
  if (srcIndex == nPE ) srcIndex  = 0;
  if (destIndex == -1 ) destIndex = nPE - 1;
  if (destIndex == nPE) destIndex = 0;

#ifdef DEBUG_SCAN
  sprintf(estring,"Physical shift, srcIndex = %d, destIndex = %d \n", srcIndex, destIndex);
  pgslib_output_c(estring);        
#endif

  tag       = nPE + 1;
  
  /* We made sure that src and dest were in range, so don't need to guard
     the sends and receives.*/
  MPI_Isend(Src_Data,
	    *Count *(SIZE_MULTIPLIER),
	    PGSLIB_MPI_DATA_TYPE,
	    destIndex,
	    tag,
	    MPI_COMM_WORLD,
	    &send_request);
  MPI_Recv(Dest_Data,
	   *Count *(SIZE_MULTIPLIER),
	   PGSLIB_MPI_DATA_TYPE,
	   srcIndex,
	   tag,
	   MPI_COMM_WORLD,
	   &r_status);
  MPI_Waitall(1, &send_request, &s_status);
     
#ifdef DEBUG_SCAN
  sprintf(estring,"Shift Received = %d\n",(int)(*Dest_Data))
  pgslib_output_c(estring);         
#endif

  MPI_Barrier( MPI_COMM_WORLD );

  return;
}
#undef PGSLIB_MPI_DATA_TYPE
#undef PGSLIB_DATA_TYPE
#undef PGSLIB_ROUTINE_NAME
#undef _DIRECTION_SIGN_
#undef SIZE_MULTIPLIER
