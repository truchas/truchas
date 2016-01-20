!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE EN_GATHER_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Supply the EN (element<-node) gather routines
  !
  !=======================================================================
  use kinds, only: r8
  use gather_module,  only: GATHER
  use gs_util,        only: en_gs_init
  use mesh_parameter_module, only: ncells, nnodes, nvc
  use truchas_logging_services
  implicit none
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
#define _DATA_TYPE_     integer
#define _OP_ID_         0
#define _DST_DIMENSION_ _DIMENSION_((nvc,ncells))

#include "en_gather.fpp"

  !========== EN_GATHER_SINGLE========================================

#define _ROUTINE_NAME_ EN_GATHER_SINGLE
#define _DATA_TYPE_    real
#define _OP_ID_        0.0
#define _DST_DIMENSION_ _DIMENSION_((nvc,ncells))

#include "en_gather.fpp"

  !========== EN_GATHER_DOUBLE========================================

#define _ROUTINE_NAME_ EN_GATHER_DOUBLE
#define _DATA_TYPE_    real(r8)
#define _OP_ID_        0.0_r8
#define _DST_DIMENSION_ _DIMENSION_((nvc,ncells))

#include "en_gather.fpp"

  !========== EN_GATHER_LOG========================================

#define _ROUTINE_NAME_ EN_GATHER_LOG
#define _DATA_TYPE_    logical
#define _OP_ID_        .false.
#define _DST_DIMENSION_ _DIMENSION_((nvc,ncells))

#include "en_gather.fpp"


  !========== EN_MIN_GATHER_INT========================================

#define _ROUTINE_NAME_  EN_MIN_GATHER_INT
#define _DATA_TYPE_     integer
#define _OP_ID_         MINVAL([integer::])
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MIN

#include "en_gather.fpp"

  !========== EN_MIN_GATHER_SINGLE========================================

#define _ROUTINE_NAME_  EN_MIN_GATHER_SINGLE
#define _DATA_TYPE_     real
#define _OP_ID_         MINVAL([real::])
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MIN

#include "en_gather.fpp"

  !========== EN_MIN_GATHER_DOUBLE========================================

#define _ROUTINE_NAME_  EN_MIN_GATHER_DOUBLE
#define _DATA_TYPE_     real(r8)
#define _OP_ID_         MINVAL([real(r8)::])
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MIN

#include "en_gather.fpp"

  !========== EN_MAX_GATHER_INT========================================

#define _ROUTINE_NAME_  EN_MAX_GATHER_INT
#define _DATA_TYPE_     integer
#define _OP_ID_         MAXVAL([integer::])
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MAX

#include "en_gather.fpp"

  !========== EN_MAX_GATHER_SINGLE========================================

#define _ROUTINE_NAME_  EN_MAX_GATHER_SINGLE
#define _DATA_TYPE_     real
#define _OP_ID_         MAXVAL([real::])
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MAX

#include "en_gather.fpp"

  !========== EN_MAX_GATHER_DOUBLE========================================

#define _ROUTINE_NAME_  EN_MAX_GATHER_DOUBLE
#define _DATA_TYPE_     real(r8)
#define _OP_ID_         MAXVAL([real(r8)::])
#define _DST_DIMENSION_ _DIMENSION_((ncells))
#define _PREFIX_OP_     MAX

#include "en_gather.fpp"

  SUBROUTINE EN_OR_GATHER_LOG (Dest, Src, INITIAL_VALUE, BOUNDARY)
    !=======================================================================
    ! PURPOSE - 
    !
    !=======================================================================
    use mesh_parameter_module, only: ncells, nnodes, nvc

    ! Arguments
    logical, dimension(nnodes),   intent(IN)  :: Src
    logical, dimension(ncells),   intent(OUT) :: Dest
    logical, dimension(:), pointer, optional  :: BOUNDARY
    logical, intent(IN   ), OPTIONAL          :: INITIAL_VALUE


    ! Local Variables
    logical :: DEFAULT_VALUE
    integer :: v
    logical, dimension(nvc,ncells) :: TempDest

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

  END SUBROUTINE EN_OR_GATHER_LOG

END MODULE EN_GATHER_MODULE


