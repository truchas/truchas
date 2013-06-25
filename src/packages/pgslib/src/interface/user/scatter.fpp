!!CPP!! This file is included by pgslib_scatter_module.F90 and provides the scatter ops.

!!CPP!! $Id: scatter.fpp,v 1.2 2001/08/16 20:48:34 lally Exp $

#ifndef _ROUTINE_NAME_
#error "_ROUTINE_NAME_ must be defined before including this file"
#endif

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

#ifndef _DST_DIMENSION_
#error "_DST_DIMENSION_ must be defined before including this file"
#endif

#ifndef _SRC_DIMENSION_
#error "_SRC_DIMENSION_ must be defined before including this file"
#endif

#ifndef _INDEX_SIZE_
#error "_INDEX_SIZE_ must be defined before including this file"
#endif

#ifndef _ARRAY_ELEMENT_
#error "_ARRAY_ELEMENT_ must be defined before including this file"
#endif

#ifndef _I_LOOP_MAX_
#error "_I_LOOP_MAX_ must be defined before including this file"
#endif

#ifndef _J_LOOP_MAX_
#error "_J_LOOP_MAX_ must be defined before including this file"
#endif

#ifndef _OP_
#error "_OP_ must be defined before including this file"
#endif

#ifndef _OP_ID_
#error "_OP_ID_ must be defined before including this file"
#endif

  subroutine _ROUTINE_NAME_ (DEST, SOURCE, INDEX, TRACE, MASK)
    use pgslib_type_module
    use pgslib_globals_module
    use pgslib_gs_module
    use pgslib_timing_module
    implicit none
    _DATA_TYPE_, intent(INOUT),    &
                               _DST_DIMENSION_      :: DEST
    _DATA_TYPE_, intent(IN   ),    &
                               _SRC_DIMENSION_    :: SOURCE
    integer (PGSLib_Int_Type), intent(INOUT),    &
                               TARGET,           &
                               _SRC_DIMENSION_   :: INDEX
    logical (PGSLib_Log_Type), intent(IN   ),    &
                               OPTIONAL,         &
                               dimension(_INDEX_SIZE_)   :: MASK
    type (PGSLib_GS_Trace),    POINTER,          &
                               OPTIONAL          :: TRACE

    ! Local Variables
    integer (PGSLib_Int_Type), POINTER,          &
                               _SRC_DIMENSION_    :: Local_Index
    _DATA_TYPE_, POINTER,          &
                               dimension(:)      :: Sup_Dest, Dup_Dest
    type (PGSLib_GS_Trace),    POINTER           :: Local_Trace  
    integer (PGSLib_Int_Type)                    :: i, j

#ifdef USE_TIMERS_1
    call Enter_Routine(Scatter_STATISTICS())
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
       if (.NOT. ASSOCIATED(Trace)) Trace => PGSLib_Setup_Trace(Local_Index, SIZE(DEST), MASK=MASK)
       Local_Trace => Trace
    else
       Local_Trace => PGSLib_Setup_Trace(Local_INDEX, SIZE(DEST), MASK=MASK)
    end if
    
#ifdef USE_TIMERS_1
    ! Now that we have a trace, we can start accumulating time
    call pgslib_start_timer(Local_Trace%ScatterTotal_Timer)
#endif    

    ! We need a supplemental repository
    ALLOCATE(Sup_Dest(PGSLib_Size_Of_Sup(Local_Trace)))
    Sup_Dest = _OP_ID_

    ! Scatter whatever we can locally
    do j = 1, _J_LOOP_MAX_
       do i = 1, _I_LOOP_MAX_
          if (PRESENT(MASK)) then
             if (.NOT. MASK(_ARRAY_ELEMENT_)) cycle
          end if
          if (Local_Index(_ARRAY_ELEMENT_) > 0) then
             Dest(Local_Index(_ARRAY_ELEMENT_)) = _OP_(Dest(Local_Index(_ARRAY_ELEMENT_)), SOURCE(_ARRAY_ELEMENT_))
          else
             Sup_Dest(-Local_Index(_ARRAY_ELEMENT_)) = _OP_(Sup_Dest(-Local_Index(_ARRAY_ELEMENT_)), SOURCE(_ARRAY_ELEMENT_))
          end if
       end do
    end do

    ! We need some space to get the duplicates
    ALLOCATE(Dup_Dest(PGSLib_Size_Of_Dup(Local_Trace)))
    Dup_Dest = _OP_ID_
    Dup_Dest = PGSLib_Scatter_Buffer(Sup_Dest, Local_Trace)
    ! Finally, sum the extras into the final dest
    do i = 1, SIZE(Dup_Dest)
       Dest(PGSLib_Dup_Index(Local_Trace, i)) = _OP_(Dest(PGSLib_Dup_Index(Local_Trace, i)), Dup_Dest(i))
    end do

    ! Clean up and go home
    DEALLOCATE(Dup_Dest, Sup_Dest)
    
#ifdef USE_TIMERS_1
    call PGSLib_Stop_Timer(Local_Trace%ScatterTotal_Timer)
#endif
    
    ! If we are using a local trace, we must get rid of it
    if (.NOT. PRESENT(TRACE)) call pgslib_deallocate_trace(Local_Trace)
       
#ifdef USE_TIMERS_1
    call Exit_Routine(SCATTER_STATISTICS())
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
