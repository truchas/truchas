!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

MODULE MATL_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define the MATL and MATL_SLOT derived types and the procedures
  !   necessary to allocate/deallocate them.
  !
  !   Public Interface(s):
  !
  !     * call GATHER_VOF (m, Vof)
  !         Gather material m volume fractions from Matl; place them in Vof.
  !
  !     * call SLOT_SET (ABC, s)
  !         Initialize slot s of ABC.
  !
  ! Contains: GATHER_VOF
  !           SLOT_SET
  !
  ! Author(s): Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !
  !=======================================================================
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  implicit none
  private

  ! Public Variables
  public :: MATERIAL, MATL
  ! Public Procedures
  public :: GATHER_VOF, SLOT_RESIZE, matl_init, matl_free

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Define MATERIAL Structure
  type MATERIAL
     ! Material Identifier
     integer :: Id = 0

     ! Volume Fraction
     real(r8)   :: Vof = 0.0_r8

  end type MATERIAL

  ! Define MATL_SLOT Structure
  type MATL_SLOT
     type(MATERIAL), dimension(:), pointer :: Cell => NULL()
  end type MATL_SLOT

  integer, parameter :: max_slots = 10
  type(MATL_SLOT), dimension(max_slots), save, Target  :: Matl

  ! The matl_init subroutine sets this value, which is intended for use
  ! only by this module and matl_utilities. This is a temporary hack
  ! until the matl structure can be redesigned.
  integer, protected, public :: ncells = 0

  ! Public variables formerly in parameter_module
  integer, public :: nmat, mat_slot = 0, mat_slot_new = 0

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! From base_types_a_allocate
  subroutine matl_init(n)
    integer, intent(in) :: n
    ncells = n
    mat_slot_new = merge(1, 2, nmat <= 1)
    call slot_resize(matl, mat_slot, mat_slot_new)
  end subroutine

  ! From base_types_a_deallocate
  subroutine matl_free
    mat_slot_new = 0
    call slot_resize(matl, mat_slot, mat_slot_new)
  end subroutine

  SUBROUTINE GATHER_VOF (m, Vof)
    !=======================================================================
    ! Purpose(s):
    !   Gather material m volume fractions from the Matl derived
    !   type and place them into array Vof.
    !=======================================================================

    ! Argument List
    integer,  intent(IN)  :: m
    real(r8), intent(OUT) :: Vof(:)

    ! Local Variables
    integer :: s

    ASSERT(size(vof) == size(matl(1)%cell))

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize
    Vof = 0.0_r8

    ! Loop over slots, gathering material m volume fractions
    do s = 1,mat_slot
       where (Matl(s)%Cell%Id == m) Vof = Matl(s)%Cell%Vof
    end do

  END SUBROUTINE GATHER_VOF

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  subroutine slot_resize (abc, slot, slot_new)

    type(matl_slot), intent(inout) :: abc(:)
    integer,         intent(inout) :: slot
    integer,         intent(in)    :: slot_new

    integer :: s
    character(128) :: message

    !! Check for invalid values for SLOT_NEW.
    if (slot_new < 0 .or. slot_new >= size(abc)) then
       write (message, 10) slot_new, size(abc)
10     format ('SLOT_RESIZE: unacceptable value of ',i2,' for slot_new; must be >= 0 and <= ',i2)
       call TLS_panic (message)
    end if

    if (slot_new < slot) then
      do s = slot_new+1, slot
        deallocate(abc(s)%cell)
      end do
    else if (slot_new > slot) then
      do s = slot+1, slot_new
        allocate(abc(s)%cell(ncells))
        call slot_set (abc, s)
      end do
    end if

    slot = slot_new

  end subroutine slot_resize

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SLOT_SET (ABC, s)
    !=======================================================================
    ! Purpose(s):
    !   Initialize the ABC structure for slot s.
    !=======================================================================

    ! Argument List
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: ABC
    integer, intent(IN) :: s

    ! Local Variables
    integer :: Id_def  = 0
    real(r8)   :: Vof_def = 0.0_r8

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Assign defaults
    ABC(s)%Cell%Id      = Id_def
    ABC(s)%Cell%Vof     = Vof_def

  END SUBROUTINE SLOT_SET
  
END MODULE MATL_MODULE
