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

#ifndef _OP_NAME_
#error "_OP_NAME_ must be defined before including this file"
#endif

#ifndef _SRC_DIMENSION_
#error "_SRC_DIMENSION_ must be defined before including this file"
#endif

#ifndef _DST_DIMENSION_
#error "_DST_DIMENSION_ must be defined before including this file"
#endif

#define _DIMENSION_(D)   dimension D

  SUBROUTINE _ROUTINE_NAME_ (Dest, Src, Node, BOUNDARY)
    !=======================================================================
    ! PURPOSE - 
    !
    !=======================================================================
    use gs_info_module, only: EN, EN_TRACE
    use mesh_module,    only: Mesh
    use pgslib_module,  only: PGSLib_GS_Trace_Setup_P

    ! Incoming Arguments
    _DATA_TYPE_ , intent(IN   ),                               &
                  _SRC_DIMENSION_                              :: Src
    _DATA_TYPE_ , intent(  OUT),                               &
                  _DST_DIMENSION_                              :: Dest
    _DATA_TYPE_ , dimension(:),                                &         
                  POINTER,                                     &
                  OPTIONAL                                     :: BOUNDARY
    integer, intent(IN), optional :: node

    ! Local Variables
    logical :: present_mask
    logical, dimension(:,:), allocatable :: Node_Mask
    integer :: memstat, v


    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Initialize relevant quantities
    Dest = _OP_ID_

    IF (.NOT. PGSLib_GS_Trace_Setup_P(EN_Trace)) Call EN_GS_Init()

    ! See if node mask is present; if so, then set it
    present_mask = PRESENT(node)
    if (present_mask) then
       allocate (Node_Mask(nvc,ncells), STAT = memstat)
       do v = 1,nvc
          Node_Mask(v,:) = (v == node)
       end do

       call _OP_NAME_ (Dest, Src, Mesh, TYPE = EN, TRACE = EN_TRACE, &
                                MASK = Node_Mask, BOUNDARY = BOUNDARY)

       if (present_mask) deallocate (Node_Mask, STAT = memstat)
    else
       call _OP_NAME_ (Dest, Src, Mesh, TYPE = EN, TRACE = EN_TRACE, &
                                BOUNDARY = BOUNDARY)
    end if

  END SUBROUTINE _ROUTINE_NAME_

#undef _ROUTINE_NAME_
#undef _DATA_TYPE_
#undef _OP_ID_
#undef _OP_NAME_
#undef _DST_DIMENSION_
#undef _SRC_DIMENSION_
