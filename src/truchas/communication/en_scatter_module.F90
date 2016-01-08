!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE EN_SCATTER_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Supply the EN (element->node) scatter routines
  !
  !=======================================================================
  use gs_util, only: EN_GS_INIT
  use kinds, only: r8
  use parameter_module
  use scatter_module, only: SUM_SCATTER, MIN_SCATTER, MAX_SCATTER, OR_SCATTER
  implicit none
  private

  Public :: EN_SUM_Scatter, &
            EN_MIN_Scatter, &
            EN_MAX_Scatter, &
            EN_OR_Scatter


  ! Interface blocks
  INTERFACE EN_SUM_SCATTER
     !=======================================================================
     ! PURPOSE - 
     !   Scatter  scalar cell-centered data to cell vertices 
     !   with an ADDITION
     !   
     !    Input: 
     !          Mesh - Mesh connectivity structure
     !          Src  - cell quantity
     !   Output: 
     !          Dest - accumulated vertex quantity 
     !          
     !=======================================================================
     MODULE PROCEDURE EN_SUM_SCATTER_SCALAR_INT
     MODULE PROCEDURE EN_SUM_SCATTER_SCALAR_SINGLE
     MODULE PROCEDURE EN_SUM_SCATTER_SCALAR_DOUBLE
     MODULE PROCEDURE EN_SUM_SCATTER_VECTOR_INT
     MODULE PROCEDURE EN_SUM_SCATTER_VECTOR_SINGLE
     MODULE PROCEDURE EN_SUM_SCATTER_VECTOR_DOUBLE
  END INTERFACE

  INTERFACE EN_MIN_SCATTER
     !=======================================================================
     ! PURPOSE - 
     !   Scatter  scalar cell-centered data to cell vertices 
     !   with an MIN
     !   
     !    Input: 
     !          Mesh - Mesh connectivity structure
     !          Src  - cell quantity
     !   Output: 
     !          Dest - accumulated vertex quantity 
     !          
     !=======================================================================
     MODULE PROCEDURE EN_MIN_SCATTER_SCALAR_INT
     MODULE PROCEDURE EN_MIN_SCATTER_SCALAR_SINGLE
     MODULE PROCEDURE EN_MIN_SCATTER_SCALAR_DOUBLE
     MODULE PROCEDURE EN_MIN_SCATTER_VECTOR_INT
     MODULE PROCEDURE EN_MIN_SCATTER_VECTOR_SINGLE
     MODULE PROCEDURE EN_MIN_SCATTER_VECTOR_DOUBLE
  END INTERFACE

  INTERFACE EN_MAX_SCATTER
     !=======================================================================
     ! PURPOSE - 
     !   Scatter  scalar cell-centered data to cell vertices 
     !   with an MAX
     !   
     !    Input: 
     !          Mesh - Mesh connectivity structure
     !          Src  - cell quantity
     !   Output: 
     !          Dest - accumulated vertex quantity 
     !          
     !=======================================================================
     MODULE PROCEDURE EN_MAX_SCATTER_SCALAR_INT
     MODULE PROCEDURE EN_MAX_SCATTER_SCALAR_SINGLE
     MODULE PROCEDURE EN_MAX_SCATTER_SCALAR_DOUBLE
     MODULE PROCEDURE EN_MAX_SCATTER_VECTOR_INT
     MODULE PROCEDURE EN_MAX_SCATTER_VECTOR_SINGLE
     MODULE PROCEDURE EN_MAX_SCATTER_VECTOR_DOUBLE
  END INTERFACE

  INTERFACE EN_OR_SCATTER
     !=======================================================================
     ! PURPOSE - 
     !   Scatter  scalar cell-centered data to cell vertices 
     !   with an OR
     !   
     !    Input: 
     !          Mesh - Mesh connectivity structure
     !          Src  - cell quantity
     !   Output: 
     !          Dest - accumulated vertex quantity 
     !          
     !=======================================================================
     MODULE PROCEDURE EN_OR_SCATTER_SCALAR_LOG
     MODULE PROCEDURE EN_OR_SCATTER_VECTOR_LOG
  END INTERFACE


CONTAINS

  !========== EN_SUM_SCATTER_SCALAR_INT  ==============================
#define _ROUTINE_NAME_  EN_SUM_SCATTER_SCALAR_INT
#define _DATA_TYPE_     integer
#define _OP_ID_         0
#define _OP_NAME_       sum_scatter
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_SUM_SCATTER_SCALAR_SINGLE  ==============================
#define _ROUTINE_NAME_  EN_SUM_SCATTER_SCALAR_SINGLE
#define _DATA_TYPE_     real
#define _OP_ID_         0.0
#define _OP_NAME_       sum_scatter
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_SUM_SCATTER_SCALAR_DOUBLE  ==============================
#define _ROUTINE_NAME_  EN_SUM_SCATTER_SCALAR_DOUBLE
#define _DATA_TYPE_     real(r8)
#define _OP_ID_         0.0_r8
#define _OP_NAME_       sum_scatter
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_SUM_SCATTER_VECTOR_INT  ==============================
#define _ROUTINE_NAME_  EN_SUM_SCATTER_VECTOR_INT
#define _DATA_TYPE_     integer
#define _OP_ID_         0
#define _OP_NAME_       sum_scatter
#define _SRC_DIMENSION_ _DIMENSION_((nvc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_SUM_SCATTER_VECTOR_SINGLE  ==============================
#define _ROUTINE_NAME_  EN_SUM_SCATTER_VECTOR_SINGLE
#define _DATA_TYPE_     real
#define _OP_ID_         0.0
#define _OP_NAME_       sum_scatter
#define _SRC_DIMENSION_ _DIMENSION_((nvc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_SUM_SCATTER_VECTOR_DOUBLE  ==============================
#define _ROUTINE_NAME_  EN_SUM_SCATTER_VECTOR_DOUBLE
#define _DATA_TYPE_     real(r8)
#define _OP_ID_         0.0_r8
#define _OP_NAME_       sum_scatter
#define _SRC_DIMENSION_ _DIMENSION_((nvc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_MIN_SCATTER_SCALAR_INT  ==============================
#define _ROUTINE_NAME_  EN_MIN_SCATTER_SCALAR_INT
#define _DATA_TYPE_     integer
#define _OP_ID_         MINVAL( (/0/), MASK=.false.)
#define _OP_NAME_       min_scatter
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_MIN_SCATTER_SCALAR_SINGLE  ==============================
#define _ROUTINE_NAME_  EN_MIN_SCATTER_SCALAR_SINGLE
#define _DATA_TYPE_     real
#define _OP_ID_         MINVAL( (/0.0/), MASK=.false.)
#define _OP_NAME_       min_scatter
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_MIN_SCATTER_SCALAR_DOUBLE  ==============================
#define _ROUTINE_NAME_  EN_MIN_SCATTER_SCALAR_DOUBLE
#define _DATA_TYPE_     real(r8)
#define _OP_ID_         MINVAL( (/0.0_r8/), MASK=.false.)
#define _OP_NAME_       min_scatter
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_MIN_SCATTER_VECTOR_INT  ==============================
#define _ROUTINE_NAME_  EN_MIN_SCATTER_VECTOR_INT
#define _DATA_TYPE_     integer
#define _OP_ID_         MINVAL( (/0/), MASK=.false.)
#define _OP_NAME_       min_scatter
#define _SRC_DIMENSION_ _DIMENSION_((nvc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"


  !========== EN_MIN_SCATTER_VECTOR_SINGLE  ==============================
#define _ROUTINE_NAME_  EN_MIN_SCATTER_VECTOR_SINGLE
#define _DATA_TYPE_     real
#define _OP_ID_         MINVAL( (/0.0/), MASK=.false.)
#define _OP_NAME_       min_scatter
#define _SRC_DIMENSION_ _DIMENSION_((nvc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_MIN_SCATTER_VECTOR_DOUBLE  ==============================
#define _ROUTINE_NAME_  EN_MIN_SCATTER_VECTOR_DOUBLE
#define _DATA_TYPE_     real(r8)
#define _OP_ID_         MINVAL( (/0.0_r8/), MASK=.false.)
#define _OP_NAME_       min_scatter
#define _SRC_DIMENSION_ _DIMENSION_((nvc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_MAX_SCATTER_SCALAR_INT  ==============================
#define _ROUTINE_NAME_  EN_MAX_SCATTER_SCALAR_INT
#define _DATA_TYPE_     integer
#define _OP_ID_         MAXVAL( (/0/), MASK=.false.)
#define _OP_NAME_       max_scatter
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_MAX_SCATTER_SCALAR_SINGLE  ==============================
#define _ROUTINE_NAME_  EN_MAX_SCATTER_SCALAR_SINGLE
#define _DATA_TYPE_     real
#define _OP_ID_         MAXVAL( (/0.0/), MASK=.false.)
#define _OP_NAME_       max_scatter
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_MAX_SCATTER_SCALAR_DOUBLE  ==============================
#define _ROUTINE_NAME_  EN_MAX_SCATTER_SCALAR_DOUBLE
#define _DATA_TYPE_     real(r8)
#define _OP_ID_         MAXVAL( (/0.0_r8/), MASK=.false.)
#define _OP_NAME_       max_scatter
#define _SRC_DIMENSION_ _DIMENSION_((ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_MAX_SCATTER_VECTOR_INT  ==============================
#define _ROUTINE_NAME_  EN_MAX_SCATTER_VECTOR_INT
#define _DATA_TYPE_     integer
#define _OP_ID_         MAXVAL( (/0/), MASK=.false.)
#define _OP_NAME_       max_scatter
#define _SRC_DIMENSION_ _DIMENSION_((nvc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"


  !========== EN_MAX_SCATTER_VECTOR_SINGLE  ==============================
#define _ROUTINE_NAME_  EN_MAX_SCATTER_VECTOR_SINGLE
#define _DATA_TYPE_     real
#define _OP_ID_         MAXVAL( (/0.0/), MASK=.false.)
#define _OP_NAME_       max_scatter
#define _SRC_DIMENSION_ _DIMENSION_((nvc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  !========== EN_MAX_SCATTER_VECTOR_DOUBLE  ==============================
#define _ROUTINE_NAME_  EN_MAX_SCATTER_VECTOR_DOUBLE
#define _DATA_TYPE_     real(r8)
#define _OP_ID_         MAXVAL( (/0.0_r8/), MASK=.false.)
#define _OP_NAME_       max_scatter
#define _SRC_DIMENSION_ _DIMENSION_((nvc, ncells))
#define _DST_DIMENSION_ _DIMENSION_((nnodes))

#include "en_scatter.fpp"

  SUBROUTINE EN_OR_SCATTER_SCALAR_LOG (Dest, Src, BOUNDARY)
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================
    use gs_info_module,   only: EN, EN_TRACE
    use mesh_module,      only: Mesh
    use parameter_module, only: ncells, nnodes

    ! Arguments
    logical, dimension(ncells), intent(IN)   :: Src
    logical, dimension(nnodes), intent(OUT)  :: Dest
    logical, dimension(:), pointer, optional :: BOUNDARY
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize the vertex quantity
    Dest = .false.

    call OR_SCATTER (Dest, Src, Mesh, TYPE = EN, TRACE = EN_TRACE, BOUNDARY = BOUNDARY)

  END SUBROUTINE EN_OR_SCATTER_SCALAR_LOG

  SUBROUTINE EN_OR_SCATTER_VECTOR_LOG (Dest, Src, BOUNDARY)
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================
    use gs_info_module,   only: EN, EN_TRACE
    use mesh_module,      only: Mesh
    use parameter_module, only: ncells, nnodes, nvc

    ! Arguments
    logical, dimension(nvc, ncells), intent(IN)  :: Src
    logical, dimension(nnodes),      intent(OUT) :: Dest
    logical, dimension(:), pointer, optional     :: BOUNDARY
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize the vertex quantity
    Dest = .false.

    call OR_SCATTER (Dest, Src, Mesh, TYPE = EN, TRACE = EN_TRACE, BOUNDARY = BOUNDARY)

  END SUBROUTINE EN_OR_SCATTER_VECTOR_LOG

END MODULE EN_SCATTER_MODULE


