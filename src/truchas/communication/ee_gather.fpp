! This is included by gs_module for the top level gather routines.
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

  SUBROUTINE _ROUTINE_NAME_ (Dest, Src, INITIAL_VALUE, BOUNDARY, &
                                        OVERWRITE_MASKED_VALUES)
    !=======================================================================
    ! PURPOSE - 
    !
    !=======================================================================
    use gs_info_module, only: EE_TRACE, EE
    use mesh_module,    only: Mesh
    use pgslib_module,  only: PGSLib_GS_Trace_Setup_P

    ! Incoming Arguments
    _DATA_TYPE_ , intent(IN   ),                               &
                  _SRC_DIMENSION_                              :: Src
    _DATA_TYPE_ , intent(INOUT),                               &
                  _DST_DIMENSION_                              :: Dest
    _DATA_TYPE_, intent(IN   ), OPTIONAL                       :: INITIAL_VALUE
    _DATA_TYPE_ , _BDY_DIMENSION_,                             &         
                  POINTER,                                     &
                  OPTIONAL                                     :: BOUNDARY
    LOGICAL, intent(IN), OPTIONAL      :: OVERWRITE_MASKED_VALUES



    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialization of Dest is no longer done in here.  Instead it is done
    ! in GATHER.

    if (.not. EE_MASK_Initialized) then
       EE_MASK_Initialized = .true.
       call GS_INIT_EE_MASK()
    end if

    if (.not. PGSLib_GS_Trace_Setup_P(EE_Trace)) call EE_GS_INIT()

    call GATHER (Dest, Src, Mesh, TYPE = EE, TRACE = EE_TRACE, &
                 MASK = El_Nbr_Mask, BOUNDARY = BOUNDARY, INITIAL_VALUE = INITIAL_VALUE, &
                 OVERWRITE_MASKED_VALUES = OVERWRITE_MASKED_VALUES)

  END SUBROUTINE _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _OP_ID_
#undef _SRC_DIMENSION_
#undef _DST_DIMENSION_
#undef _BDY_DIMENSION_
