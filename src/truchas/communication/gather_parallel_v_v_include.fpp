! This is included for the gather routines.
! This is the parallel version


#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _OP_ID_
#error "_OP_ID_ must be defined before including this file"
#endif

SUBROUTINE _ROUTINE_NAME_ (Dest, Src, Mesh, TYPE, TRACE, BOUNDARY, MASK, &
                           INITIAL_VALUE, OVERWRITE_MASKED_VALUES)
  !=======================================================================
  ! PURPOSE - 
  !   Gather either cell-centered data or vertex centered _DATA_TYPE_ data
  !   into an _DATA_TYPE_ cell-centered vector of length ShortDim
  !=======================================================================

  ! Global scalars & arrays
  _DATA_TYPE_             , dimension(:,:),                        &
       intent(INOUT)                          :: Dest
  _DATA_TYPE_             , dimension(:,:),                        &
       intent(IN   )                          :: Src
  type(MESH_CONNECTIVITY),  dimension(SIZE(Dest,2)),               &
       intent(IN)                             :: Mesh
  _DATA_TYPE_             , dimension(:,:),                        &         
       POINTER,                               &
       OPTIONAL                               :: BOUNDARY
  logical,                  dimension(SIZE(Dest,1), SIZE(Dest,2)), &
       intent(IN ),                           &
       OPTIONAL                               :: MASK
  type (COMM_TYPE),         intent(IN   )                          :: TYPE
   _DATA_TYPE_, intent(IN   ), OPTIONAL                       :: INITIAL_VALUE
  LOGICAL, intent(IN), OPTIONAL      :: OVERWRITE_MASKED_VALUES

  type (PGSLib_GS_Trace), POINTER :: Trace

  ! Local scalars & arrays
  integer :: f, c
  logical :: PRESENT_MASK, TEMP_BOUNDARY, NEW_BOUNDARY
  _DATA_TYPE_ :: DEFAULT_VALUE

  _DATA_TYPE_, POINTER, Dimension(:,:) :: Supplement_Data
  _DATA_TYPE_, POINTER, Dimension(:,:) :: Duplicate_Data
  logical :: OVERWRITE_MASKED_VALUES_Local

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Initialize relevant quantities
  PRESENT_MASK = PRESENT(MASK)

  ! We may be gathering under a mask.  If so, what should we do with the
  ! values which have mask == .FALSE.?  Either we can overwrite the value 
  ! with some a default value, or we can leave it alone.  The first case is useful
  ! and assumed in many places in T.  OTOH, the second case is useful if
  ! we know that the value already in there is what we want.
  ! The choice is determined by the logical OVERWRITE_MASKED_VALUES_Local
  ! If this is .FALSE., then we do not overwrite masked values.  If this
  ! is .TRUE. we do. 
  ! The default is .TRUE., but may be changed by supplying
  ! OVERWRITE_MASKED_VALUES as an argument.
  ! If we choose to overwrite, then the masked values get set to DEFAULT_VALUE.
  ! This is the identity value for each data type, or it is set to INITIAL_VALUE
  ! if that argument is supplied.

  if (PRESENT(OVERWRITE_MASKED_VALUES)) then
     OVERWRITE_MASKED_VALUES_Local = OVERWRITE_MASKED_VALUES
  else
     OVERWRITE_MASKED_VALUES_Local = .TRUE.

  end if

  if (PRESENT(INITIAL_VALUE) ) then
     DEFAULT_VALUE = INITIAL_VALUE
  else
     DEFAULT_VALUE = _OP_ID_
  end if

  ! If BOUNDARY is not passed in, then we have to use a local
  ! boundaray.  If BOUNDARY is passed in, then if it is NULL
  ! we have to allocate space.  In either case, after we allocate
  ! we have to do a global gather.  If BOUNDARY is passed in
  ! and is associated, we assume it is the right size and
  ! use it for supplement_data.

  ! Do we use local boundary, are use user provided boundary?
  TEMP_BOUNDARY = .NOT. PRESENT(BOUNDARY)

  ! If the BOUNDARY is associated, use it, otherwise get new boundary data.
  IF (TEMP_BOUNDARY) THEN
     NEW_BOUNDARY   = .TRUE.
  ELSE
     NEW_BOUNDARY   = .NOT. ASSOCIATED(BOUNDARY)
  ENDIF

  IF (NEW_BOUNDARY) THEN
     ALLOCATE(Supplement_Data(SIZE(Src,1),PGSLib_Size_Of_Sup(Trace)))
     ! If we were passed in a boundary, and it wasn't allocated
     ! then point boundary at newly allocated Supplement_Data and save it
     IF (.NOT. TEMP_BOUNDARY) BOUNDARY => Supplement_Data

     ! Gather into comm buffer
     ALLOCATE(Duplicate_Data(SIZE(Src,1), PGSLib_Size_Of_Dup(Trace)))
     do c = 1, PGSLib_Size_Of_Dup(Trace)
        do f = 1, SIZE(Src,1)
           Duplicate_Data(f,c) = Src(f,PGSLib_Dup_Index(Trace, c))
        enddo
     enddo

     ! The communication, takes data from Duplicate buffer, puts data into Supplement buffer
     Supplement_Data = PGSLib_gather_buffer(Duplicate_Data, Trace)

     ! Clean up done with Duplicate_Data
     DEALLOCATE(Duplicate_Data)

  ELSE
     ! If we were passed in a boundary, and it was allocated already,
     ! then us it
     IF (.NOT. TEMP_BOUNDARY) Supplement_Data => BOUNDARY 

     IF ( (SIZE(Supplement_Data,1) /= SIZE(Src,1))  .AND. &
          (SIZE(Supplement_Data,1) /= PGSLib_Size_Of_Sup(Trace)) ) THEN
        call TLS_panic ('GATHER: wrong size boundary data')
     END IF
  END IF

  ! Loop over neighbors, gathering into Dest
  SELECT CASE(Type%Type)
  CASE(EE%Type)
     IF (PRESENT_MASK) THEN
        do c = 1, SIZE(Dest,2)
           do f = 1,SIZE(Dest,1)
              if (.NOT. Mask(f,c) ) then
                 if (OVERWRITE_MASKED_VALUES_Local) Dest(f,c) = DEFAULT_VALUE
                 CYCLE
              endif
              if (Mesh(c)%Ngbr_cell(f) > 0) then
                 Dest(f,c) =             Src(Mesh(c)%Ngbr_face(f),  Mesh(c)%Ngbr_cell(f) )
              else
                 Dest(f,c) = Supplement_Data(Mesh(c)%Ngbr_face(f), -Mesh(c)%Ngbr_cell(f) )
              end if
           end do
        end do
     ELSE
        do c = 1, SIZE(Dest,2)
           do f = 1,SIZE(Dest,1)
              if (Mesh(c)%Ngbr_cell(f) > 0) then
                 Dest(f,c) =             Src(Mesh(c)%Ngbr_face(f),  Mesh(c)%Ngbr_cell(f) )
              else
                 Dest(f,c) = Supplement_Data(Mesh(c)%Ngbr_face(f), -Mesh(c)%Ngbr_cell(f) )
              end if
           end do
        end do
     ENDIF

  CASE(EN%Type)
     call TLS_panic ('GATHER_V_V_?: Vector Element<->Node Gathers are not yet supported')

  END SELECT

  ! If BOUNDARY was not passed in, then deallocate space provided
  IF (TEMP_BOUNDARY) DEALLOCATE(Supplement_Data)
  return

END SUBROUTINE _ROUTINE_NAME_
  
#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _OP_ID_

