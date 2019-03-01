!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BC_INFO_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define procedures which provide various information about
  !   boundary condition types and locations.
  !
  ! Contains: BCMatch_Scalar
  !           BCMatch_Array
  !           Boundary_Scalar
  !           Boundary_Array
  !           ExternalBoundary_Scalar
  !           ExternalBoundary_Array
  !           InternalBoundary_Scalar
  !           InternalBoundary_Array
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !            Bryan R. Lally, LANL ESA-EPE (lally@lanl.gov)
  !
  !=======================================================================
  use bc_type_module, only: BC_STRUCTURE
  use legacy_mesh_api, only: ncells
  implicit none

  ! Private Module
  private

  ! Public Variables

  ! Public Subroutines
  public :: BCMatch, Boundary, ExternalBoundary, InternalBoundary

  ! Generic Procedure Interfaces
  interface BCMatch
     module procedure BCMatch_Scalar
     module procedure BCMatch_Array
  end interface

  interface Boundary
     module procedure Boundary_Scalar
     module procedure Boundary_Array
  end interface

  interface ExternalBoundary
     module procedure ExternalBoundary_Scalar
     module procedure ExternalBoundary_Array
  end interface

  interface InternalBoundary
     module procedure InternalBoundary_Scalar
     module procedure InternalBoundary_Array
  end interface

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  FUNCTION BCMatch_Scalar (BC, face, cell)
    !=======================================================================
    ! Purpose(s):
    !   Return BC type identifier; returns a scalar.
    !=======================================================================

    ! Argument List
    type(BC_STRUCTURE),       intent(IN) :: BC
    integer, intent(IN) :: cell
    integer, intent(IN) :: face

    ! Local Variables
    integer :: shift

    ! Function Return
    integer :: BCMatch_Scalar

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    shift = (face - 1) * BC%shift
    BCMatch_Scalar = ISHFT(IAND(BC%Flag(cell),ISHFT(BC%Exist,shift)),-shift)

    return

  END FUNCTION BCMatch_Scalar
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION BCMatch_Array (BC, face)
    !=======================================================================
    ! Purpose(s):
    !   Return BC type identifier; returns an array.
    !=======================================================================

    ! Argument List
    type(BC_STRUCTURE),       intent(IN) :: BC
    integer, intent(IN) :: face

    ! Local Variables
    integer :: shift

    ! Function Return
    integer, dimension(ncells) :: BCMatch_Array

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    shift = (face - 1) * BC%shift
    BCMatch_Array = BC%Exist
    BCMatch_Array = ISHFT(IAND(BC%Flag,ISHFT(BCMatch_Array,shift)),-shift)

    return

  END FUNCTION BCMatch_Array

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION Boundary_Scalar (BC, face, cell)
    !=======================================================================
    ! Purpose(s):
    !   Returns a true if a cell face is on a boundary, else returns a false;
    !   returns a scalar.
    !=======================================================================

    ! Argument List
    type(BC_STRUCTURE),       intent(IN) :: BC
    integer, intent(IN) :: cell
    integer, intent(IN) :: face

    ! Local Variables

    ! Function Return
    logical :: Boundary_Scalar

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    if (IAND(BC%Flag(cell),ISHFT(BC%EXIST,(face-1)*BC%Shift)) /= 0) Then
       Boundary_Scalar = .true.
    else
       Boundary_Scalar = .false.
    end if

    return

  END FUNCTION Boundary_Scalar
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION Boundary_Array (BC, face)
    !=======================================================================
    ! Purpose(s):
    !   Returns a true if a cell face is on a boundary, else returns a false;
    !   returns a mask array.
    !=======================================================================

    ! Argument List
    type(BC_STRUCTURE),       intent(IN) :: BC
    integer, intent(IN) :: face

    ! Local Variables

    ! Function Return
    logical, dimension(ncells) :: Boundary_Array

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    where (IAND(BC%Flag,ISHFT(BC%EXIST,(face-1)*BC%Shift)) /= 0)
       Boundary_Array = .true.
    elsewhere
       Boundary_Array = .false.
    end where

    return

  END FUNCTION Boundary_Array

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION ExternalBoundary_Scalar (BC, face, cell)
    !=======================================================================
    ! Purpose(s):
    !   Returns a true if a cell face is on an external boundary, else
    !   returns a false; returns a scalar.
    !=======================================================================

    ! argument list
    type(BC_STRUCTURE),       intent(IN) :: BC
    integer, intent(IN) :: cell
    integer, intent(IN) :: face
    
    ! function return
    logical :: ExternalBoundary_Scalar

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return

    ExternalBoundary_Scalar = Boundary(BC,face,cell) .and. &
         .not. BTEST(BC%Flag(cell),face*BC%Shift-1)

    return

  END FUNCTION ExternalBoundary_Scalar

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  Function ExternalBoundary_Array (BC, face)
    !=======================================================================
    ! Purpose:
    !   Returns a true if a cell face is on an external boundary, else
    !   returns a false; returns a mask array.
    !=======================================================================

    ! Argument List
    type(BC_STRUCTURE),       intent(IN) :: BC
    integer, intent(IN) :: face

    ! Local Variables
    integer :: cell

    ! Function Return
    logical, dimension(ncells) :: ExternalBoundary_Array

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! function return
    Do cell = 1, ncells
       ExternalBoundary_Array(cell) = ExternalBoundary_Scalar(BC,face,cell)
    End Do

    Return

  END FUNCTION ExternalBoundary_Array

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION InternalBoundary_Scalar (BC, face, cell)
    !=======================================================================
    ! Purpose(s):
    !   Returna a true if a cell face is on an internal boundary, else
    !   returns a false; returns a scalar.
    !=======================================================================

    ! Argument List
    type(BC_STRUCTURE),       intent(IN) :: BC
    integer, intent(IN) :: cell
    integer, intent(IN) :: face

    ! Local Variables

    ! Function Return
    logical :: InternalBoundary_Scalar

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    InternalBoundary_Scalar = Boundary(BC,face,cell) .and. &
         BTEST(BC%Flag(cell),face*BC%Shift-1)

    return

  END FUNCTION InternalBoundary_Scalar

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION InternalBoundary_Array (BC, face)
    !=======================================================================
    ! Purpose(s):
    !   Returns a true if a cell face is on an internal boundary, else
    !   returns a false; returns a mask array.
    !=======================================================================

    ! Argument List
    type(BC_STRUCTURE),       intent(IN) :: BC
    integer, intent(IN) :: face

    ! Local Variables
    integer :: cell

    ! Function Return
    logical, dimension(ncells) :: InternalBoundary_Array

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    Do cell = 1, ncells
       InternalBoundary_Array(cell) = InternalBoundary_Scalar(BC,face,cell)
    End Do

    Return

  END FUNCTION InternalBoundary_Array

END MODULE BC_INFO_MODULE
