MODULE EN_GATHER_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Supply the EN (element<-node) gather routines
  !
  !=======================================================================
  use gather_module,  only: GATHER
  use gs_util,        only: en_gs_init
  use kind_module
  use parameter_module
  implicit none

  ! Private Module
  private

  Public :: EN_Gather,      &
            EN_MIN_Gather,  &
            EN_MAX_Gather,  &
            EN_OR_Gather


  ! Interface blocks
  INTERFACE EN_GATHER
     !=======================================================================
     ! PURPOSE - 
     !   Gather real vertex data to a real cell-centered
     !   vector of length nvc 
     !   
     !    Input: 
     !          Mesh - Mesh connectivity structure
     !          Src  - Vertex-centered quantity 
     !   Output: 
     !          Dest - Cell-centered vector of length nvc containing
     !                 the value of Src for all vertex neighbors
     !          
     !=======================================================================
     MODULE PROCEDURE EN_GATHER_INT
     MODULE PROCEDURE EN_GATHER_SINGLE
     MODULE PROCEDURE EN_GATHER_DOUBLE
     MODULE PROCEDURE EN_GATHER_LOG
  END INTERFACE

  INTERFACE EN_MIN_GATHER
     !=======================================================================
     ! PURPOSE -
     !   Find the minimum real vertex quantity and store
     !   it in a scalar cell-centered quantity
     !
     !    Input: Mesh - Mesh connectivity structure
     !           Src  - Vertex-centered quantity 
     !   Output: Dest - Cell-centered quantity 
     !
     !=======================================================================
     MODULE PROCEDURE EN_MIN_GATHER_INT
     MODULE PROCEDURE EN_MIN_GATHER_SINGLE
     MODULE PROCEDURE EN_MIN_GATHER_DOUBLE
  END INTERFACE

  INTERFACE EN_MAX_GATHER
     !=======================================================================
     ! PURPOSE -
     !   Find the maximum real vertex quantity and store
     !   it in a scalar cell-centered quantity
     !
     !    Input: Mesh - Mesh connectivity structure
     !           Src  - Vertex-centered quantity 
     !   Output: Dest - Cell-centered quantity 
     !
     !=======================================================================
     MODULE PROCEDURE EN_MAX_GATHER_INT
     MODULE PROCEDURE EN_MAX_GATHER_SINGLE
     MODULE PROCEDURE EN_MAX_GATHER_DOUBLE
  END INTERFACE

  INTERFACE EN_OR_GATHER
     !=======================================================================
     ! PURPOSE - 
     !   Gather the OR of logical ::vertex data and store in
     !   a scalar cell-centered quantity
     !   
     !    Input: 
     !          Mesh - Mesh connectivity structure
     !          Src  - vertex-centered quantity 
     !   Output: 
     !          Dest - cell-centered quantity 
     !          
     !=======================================================================
     MODULE PROCEDURE EN_OR_GATHER_LOG
  END INTERFACE


CONTAINS

  !========== EN_GATHER_INT========================================

#define _ROUTINE_NAME_  EN_GATHER_INT
#define _DATA_TYPE_     integer (int_kind)
#define _OP_ID_         0_int_kind
#define _DST_DIMENSION_ _DIMENSION_((nvc,ncells))

#include "en_gather.fpp"

  !========== EN_GATHER_SINGLE========================================

#define _ROUTINE_NAME_ EN_GATHER_SINGLE
#define _DATA_TYPE_    real (single_kind)
#define _OP_ID_        0.0_single_kind
#define _DST_DIMENSION_ _DIMENSION_((nvc,ncells))

#include "en_gather.fpp"

  !========== EN_GATHER_DOUBLE========================================

#define _ROUTINE_NAME_ EN_GATHER_DOUBLE
#define _DATA_TYPE_    real (double_kind)
#define _OP_ID_        0.0_double_kind
#define _DST_DIMENSION_ _DIMENSION_((nvc,ncells))

#include "en_gather.fpp"

  !========== EN_GATHER_LOG========================================

#define _ROUTINE_NAME_ EN_GATHER_LOG
#define _DATA_TYPE_    logical(log_kind)
#define _OP_ID_        .false.
#define _DST_DIMENSION_ _DIMENSION_((nvc,ncells))

#include "en_gather.fpp"


  !========== EN_MIN_GATHER_INT========================================

#define _ROUTINE_NAME_  EN_MIN_GATHER_INT
#define _DATA_TYPE_     integer (int_kind)
#define _OP_ID_         MINVAL((/0_int_kind/), MASK=.false.)
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MIN

#include "en_gather.fpp"

  !========== EN_MIN_GATHER_SINGLE========================================

#define _ROUTINE_NAME_  EN_MIN_GATHER_SINGLE
#define _DATA_TYPE_     real (single_kind)
#define _OP_ID_         MINVAL((/0.0_single_kind/), MASK=.false.)
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MIN

#include "en_gather.fpp"

  !========== EN_MIN_GATHER_DOUBLE========================================

#define _ROUTINE_NAME_  EN_MIN_GATHER_DOUBLE
#define _DATA_TYPE_     real (double_kind)
#define _OP_ID_         MINVAL((/0.0_double_kind/), MASK=.false.)
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MIN

#include "en_gather.fpp"

  !========== EN_MAX_GATHER_INT========================================

#define _ROUTINE_NAME_  EN_MAX_GATHER_INT
#define _DATA_TYPE_     integer (int_kind)
#define _OP_ID_         MAXVAL((/0_int_kind/), MASK=.false.)
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MAX

#include "en_gather.fpp"

  !========== EN_MAX_GATHER_SINGLE========================================

#define _ROUTINE_NAME_  EN_MAX_GATHER_SINGLE
#define _DATA_TYPE_     real (single_kind)
#define _OP_ID_         MAXVAL((/0.0_single_kind/), MASK=.false.)
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MAX

#include "en_gather.fpp"

  !========== EN_MAX_GATHER_DOUBLE========================================

#define _ROUTINE_NAME_  EN_MAX_GATHER_DOUBLE
#define _DATA_TYPE_     real (double_kind)
#define _OP_ID_         MAXVAL((/0.0_double_kind/), MASK=.false.)
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MAX

#include "en_gather.fpp"

  SUBROUTINE EN_OR_GATHER_LOG (Dest, Src, INITIAL_VALUE, BOUNDARY)
    !=======================================================================
    ! PURPOSE - 
    !
    !=======================================================================
    use kind_module,      only: int_kind, log_kind
    use parameter_module, only: ncells, nnodes, nvc

    implicit none

    ! Arguments
    logical(log_kind), dimension(nnodes),   intent(IN)  :: Src
    logical(log_kind), dimension(ncells),   intent(OUT) :: Dest
    logical(log_kind), dimension(:), pointer, optional  :: BOUNDARY
    logical(log_kind), intent(IN   ), OPTIONAL          :: INITIAL_VALUE


    ! Local Variables
    logical(log_kind) :: DEFAULT_VALUE
    integer(int_kind) :: v
    logical(log_kind), dimension(nvc,ncells) :: TempDest

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant quantities
    if (PRESENT(INITIAL_VALUE)) then
       DEFAULT_VALUE = INITIAL_VALUE
    else
        DEFAULT_VALUE = .FALSE.
    end if
     
    TempDest = DEFAULT_VALUE
    Dest     = DEFAULT_VALUE

    call EN_GATHER (TempDest, Src, BOUNDARY = BOUNDARY)
    
    do v = 1,nvc
       Dest = Dest .or. TempDest(v,:)
    end do

    return

  END SUBROUTINE EN_OR_GATHER_LOG


END MODULE EN_GATHER_MODULE


