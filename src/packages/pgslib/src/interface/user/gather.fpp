!!CPP!! This file is included by pgslib_gather_module.F90 and provides the gathers.

!!CPP!! $Id: gather.fpp,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $
  subroutine _ROUTINE_NAME_ (DEST, SOURCE, INDEX, TRACE, MASK)
    implicit none
    _DATA_TYPE_, intent(INOUT),    &
              _DST_DIMENSION_      :: DEST
    _DATA_TYPE_, intent(IN   ),    &
                _SRC_DIMENSION_    :: SOURCE
    integer (PGSLib_Int_Type),     &
                 intent(INOUT),    &
                 TARGET,           &
                 _DST_DIMENSION_   :: INDEX
    logical (PGSLib_Log_Type),     &
                 intent(IN   ),    &
                 OPTIONAL,         &
           dimension(_INDEX_SIZE_)   :: MASK
    type (PGSLib_GS_Trace),        &
                 POINTER,          &
                 OPTIONAL          :: TRACE

    ! Local Variables
    integer (PGSLib_Int_Type), POINTER,          &
                               _DST_DIMENSION_    :: Local_Index
    _DATA_TYPE_, POINTER,          &
                               dimension(:)      :: Sup_Src, Dup_Src
    type (PGSLib_GS_Trace),    POINTER           :: Local_Trace
    integer (PGSLib_Int_Type)                    :: i, j

#ifdef USE_TIMERS_1
    call Enter_Routine(Gather_STATISTICS())
#endif
    ! Index must be returned uncorrupted if we do not have a Trace, since
    ! in that case user expects INDEX values to be maintained.
    if (.NOT. PRESENT(Trace)) then
       ALLOCATE(Local_Index(_INDEX_SIZE_))
       Local_Index =  Index
    else
       Local_Index => Index
    end if

    ! If we got a trace, we use it.  If we did not get a trace, we use a local one
    ! If we got a trace, we may need to initialize it.
    NULLIFY(Local_Trace)
    if (PRESENT(Trace)) then
       if (.NOT. ASSOCIATED(Trace)) Trace => PGSLib_Setup_Trace(Local_Index, SIZE(SOURCE), MASK=MASK)
       Local_Trace => Trace
    else
       Local_Trace => PGSLib_Setup_Trace(Local_INDEX, SIZE(SOURCE), MASK=MASK)
    end if

#ifdef USE_TIMERS_1
    ! Now that we have a trace, we can start accumulating time
    call pgslib_start_timer(Local_Trace%GatherTotal_Timer)
#endif

    ! We need a place to gather the duplicates
    ALLOCATE(Dup_Src(PGSLib_Size_Of_Dup(Local_Trace)))
    Dup_Src = _OP_ID_

    ! Gather what will move between processors.
    Dup_Src = SOURCE(PGSLib_Dup_Index(Local_Trace))

    ! We need some space to receive the supplements
    ALLOCATE(Sup_Src(PGSLib_Size_Of_Sup(Local_Trace)))
    Sup_Src = _OP_ID_

    ! Move buffers between processors
    Sup_Src = PGSLib_Gather_Buffer(Dup_Src, Local_Trace)

    !  Gather, now everything is available
    do j = 1, _J_LOOP_MAX_
       do i = 1, _I_LOOP_MAX_
          if (PRESENT(MASK)) then
             if (.NOT. MASK(_ARRAY_ELEMENT_)) cycle
          end if
          if (Local_Index(_ARRAY_ELEMENT_) > 0) then
             Dest(_ARRAY_ELEMENT_) = SOURCE(Local_Index(_ARRAY_ELEMENT_))
          else
             Dest(_ARRAY_ELEMENT_) = Sup_Src(-Local_Index(_ARRAY_ELEMENT_))
          end if
       end do
    end do

    ! Clean up and go home
    DEALLOCATE(Dup_Src, Sup_Src)

#ifdef USE_TIMERS_1
    call PGSLib_Stop_Timer(Local_Trace%GatherTotal_Timer)
#endif

    ! If we are using a local trace, we must get get rid of it
    if (.NOT. PRESENT(TRACE)) call pgslib_deallocate_trace(Local_Trace)

#ifdef USE_TIMERS_1
    call Exit_Routine(GATHER_STATISTICS())
#endif

    RETURN

  end subroutine _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _INDEX_SIZE_
#undef _ARRAY_ELEMENT_
#undef _I_LOOP_MAX_
#undef _J_LOOP_MAX_
#undef _OP_ID_
