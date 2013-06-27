! This is included for the all_gather routines.

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif

#ifndef _DEST_DATA_TYPE_
#error "_DEST_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _SRC_DIMENSION_
#error "_SRC_DIMENSION_ must be defined before including this file"
#endif

#ifndef _DST_DIMENSION_
#error "_DST_DIMENSION_ must be defined before including this file"
#endif

#ifndef _BDY_DIMENSION_
#error "_BDY_DIMENSION_ must be defined before including this file"
#endif

#define _DIMENSION_(D) dimension D

  subroutine _ROUTINE_NAME_(DEST, SOURCE, BOUNDARY, RANGE)
    use var_vector_module
    use pgslib_module,only: PGSLib_Size_Of_Dup,   &
                          PGSLib_Size_Of_Sup,   &
                          PGSLib_Dup_Index,     &
                          PGSLib_Buffer_Gather

    ! Arguments
    _DEST_DATA_TYPE_, _DST_DIMENSION_ :: DEST
    _DATA_TYPE_, _SRC_DIMENSION_ :: SOURCE
    _DATA_TYPE_, POINTER, OPTIONAL, _BDY_DIMENSION_:: BOUNDARY

    integer, dimension(:), OPTIONAL :: RANGE

    ! Local variables
    integer :: node, nn, i, lc, uc
    logical :: TEMP_BOUNDARY, NEW_BOUNDARY
    _DATA_TYPE_, POINTER, Dimension(:) :: Supplement_Data
    _DATA_TYPE_, POINTER, Dimension(:) :: Duplicate_Data
!    integer, POINTER, Dimension(:) :: Duplicate_Indx
    integer, POINTER, Dimension(:) :: Node_Ngbrs
    _DATA_TYPE_, POINTER, Dimension(:) :: Dest_Column

    
  ! Do we use local boundary, are use user provided boundary?
  TEMP_BOUNDARY = .NOT. PRESENT(BOUNDARY)

  ! If the BOUNDARY is associated, use it, otherwise get new boundary data.
  IF (TEMP_BOUNDARY) THEN
     NEW_BOUNDARY   = .TRUE.
  ELSE
     NEW_BOUNDARY   = .NOT. ASSOCIATED(BOUNDARY)
  ENDIF

  IF (NEW_BOUNDARY) THEN
     ALLOCATE(Supplement_Data(PGSLib_Size_Of_Sup(NN_All_Ngbr_Trace)))
     ! If we were passed in a boundary, and it wasn't allocated
     ! then point boundary at newly allocated Supplement_Data and save it
     IF (.NOT. TEMP_BOUNDARY) BOUNDARY => Supplement_Data

     ! Gather into comm buffer
     ALLOCATE(Duplicate_Data(PGSLib_Size_Of_Dup(NN_All_Ngbr_Trace)))
     do i = 1, PGSLib_Size_Of_Dup(NN_All_Ngbr_Trace)
        Duplicate_Data(i) = SOURCE(PGSLib_Dup_Index(NN_All_Ngbr_Trace, i))
     end do

     ! The communication, takes data from Duplicate buffer, puts data into Supplement buffer
     call PGSLib_Buffer_Gather(Supplement_Data,  &
                               Duplicate_Data,   &
                               NN_All_Ngbr_Trace)
     ! Clean up done with Duplicate_Data
     DEALLOCATE(Duplicate_Data)

  ELSE
     ! If we were passed in a boundary, and it was allocated already,
     ! then us it
     IF (.NOT. TEMP_BOUNDARY) Supplement_Data => BOUNDARY 

     IF (SIZE(Supplement_Data,1) /= PGSLib_Size_Of_Sup(NN_All_Ngbr_Trace)) THEN
        call TLS_panic ('GATHER: wrong size boundary data')
     END IF
  END IF

  if (PRESENT(RANGE)) then
     lc = RANGE(1)
     uc = RANGE(2)
  else
     lc = 1
     uc = SIZE(Dest,1)
  end if
  
  do node = lc, uc
     ! Grab an array representation of the neighbor index
     ! and of the Destination for neighbor node
     Node_Ngbrs => FLATTEN(Vertex_Ngbr_All(node))
     Dest_Column => FLATTEN(Dest(node))

     do nn = 1, SIZES(Dest(node))
        if (Node_Ngbrs(nn) > 0) then
           Dest_Column(nn) = SOURCE(Node_Ngbrs(nn))
        else
           Dest_Column(nn) = Supplement_Data( -Node_Ngbrs(nn))
        end if
     end do
  end do
  
  ! If BOUNDARY was not passed in, then deallocate space provided
  IF (TEMP_BOUNDARY) DEALLOCATE(Supplement_Data)

 end subroutine _ROUTINE_NAME_

#undef _DEST_DATA_TYPE_
#undef _DATA_TYPE_
#undef _ROUTINE_NAME_
#undef _SRC_DIMENSION_
#undef _DST_DIMENSION_
#undef _BDY_DIMENSION_

