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

#ifndef _DST_DIMENSION_
#error "_DST_DIMENSION_ must be defined before including this file"
#endif

#define _DIMENSION_(D)   dimension D

  SUBROUTINE _ROUTINE_NAME_ (Dest, Src, INITIAL_VALUE, BOUNDARY)
    !=======================================================================
    ! PURPOSE - 
    !
    !=======================================================================
    use gs_info_module,   only: EN_TRACE, EN
    use mesh_module,      only: Mesh
    use parameter_module, only: ncells, nnodes, nvc
    use pgslib_module,  only: PGSLib_GS_Trace_Setup_P

    ! Incoming Arguments
    _DATA_TYPE_ , intent(IN   ),                               &
                  dimension(nnodes)                            :: Src
    _DATA_TYPE_ , intent(  OUT),                               &
                  _DST_DIMENSION_                              :: Dest
    _DATA_TYPE_, intent(IN   ), OPTIONAL                       :: INITIAL_VALUE
    _DATA_TYPE_ , dimension(:),                                &         
                  POINTER,                                     &
                  OPTIONAL                                     :: BOUNDARY


    ! Local Variables
    _DATA_TYPE_ :: DEFAULT_VALUE
#ifdef _PREFIX_OP_
    integer :: v
    _DATA_TYPE_ , dimension(nvc,ncells) :: TempDest
#endif

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant quantities
    IF (.NOT. PGSLib_GS_Trace_Setup_P(EN_Trace)) Call EN_GS_Init()
    
    if (PRESENT(INITIAL_VALUE)) then
       DEFAULT_VALUE = INITIAL_VALUE
    else
       DEFAULT_VALUE = _OP_ID_
    end if
    

#ifdef _PREFIX_OP_
    TempDest = DEFAULT_VALUE
#endif

#ifdef _PREFIX_OP_
#define _DEST_ TempDest
#else
#define _DEST_ Dest
#endif

    call gather (_DEST_, Src, Mesh, TYPE = EN, TRACE = EN_TRACE, BOUNDARY = BOUNDARY, INITIAL_VALUE = INITIAL_VALUE)

#undef _DEST_

#ifdef _PREFIX_OP_
    Dest     = DEFAULT_VALUE
    do v = 1, nvc
       Dest = _PREFIX_OP_(TempDest(v,:), Dest)
    end do
#endif

    return

  END SUBROUTINE _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _OP_ID_
#undef _DST_DIMENSION_
#ifdef _PREFIX_OP_
#undef _PREFIX_OP_
#endif
