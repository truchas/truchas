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
  !     * call SLOT_DECREASE (ABC, slot, slot_new)
  !         Decrease the ABC size from slot to slot_new.
  !
  !     * call SLOT_INCREASE (ABC, slot, slot_new)
  !         Increase the ABC size from slot to slot_new.
  !
  !     * call SLOT_COMPRESS (ABC, slot)
  !         Compress the slot structure of ABC(slot) to remove zeroes
  !
  !     * call SLOT_SET (ABC, s)
  !         Initialize slot s of ABC.
  !
  ! Contains: GATHER_VOF
  !           SLOT_DECREASE
  !           SLOT_INCREASE
  !           SLOT_COMPRESS
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
  public :: MATERIAL, MATL_SLOT, MATL
  ! Public Procedures
  public :: GATHER_VOF, &
            SLOT_DECREASE, SLOT_INCREASE, SLOT_COMPRESS,      &
            SLOT_SET, SLOT_RESIZE, &
            matl_init, matl_free

  interface ASSIGNMENT(=)
     module procedure MATL_ASSIGN
  end interface

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
    call slot_increase(matl, mat_slot, mat_slot_new)
  end subroutine

  ! From base_types_a_deallocate
  subroutine matl_free
    mat_slot_new = 0
    call slot_decrease(matl, mat_slot, mat_slot_new)
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

  SUBROUTINE SLOT_DECREASE (ABC, slot, slot_new)
    !=======================================================================
    ! Purpose(s):
    !   Decrease the ABC structure size from slot to slot_new.
    !=======================================================================

    ! Argument List
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: ABC
    integer, intent(INOUT) :: slot, slot_new

    ! Local Variables
    integer :: s
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Check to make sure slot_new < slot and >= 0; if not, terminate
    if (slot_new >= slot .or. slot_new < 0) then
       write (message, 10) slot_new, slot
10     format ('SLOT_DECREASE: unacceptable value of ',i2,' for slot_new; must be < ',i2,' and >= 0')
       call TLS_panic (message)
    end if

    ! Deallocate a mesh size array for each appropriate slot
    do s = slot,slot_new+1,-1

       if (ASSOCIATED(ABC(s)%Cell)) DEALLOCATE (ABC(s)%Cell)

    end do

    slot = slot_new

  END SUBROUTINE SLOT_DECREASE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SLOT_INCREASE (ABC, slot, slot_new)
    !=======================================================================
    ! Purpose(s):
    !   Increase the ABC structure size from slot to slot_new.
    !=======================================================================

    ! Argument List
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: ABC
    integer, intent(INOUT) :: slot, slot_new

    ! Local Variables
    integer :: memstat, s
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Check to make slot < slot_new <= max_slots; if not, terminate
    if (slot >= slot_new .or. slot_new > max_slots) then
       write (message, 10) slot_new, slot, max_slots
10     format ('SLOT_INCREASE: unacceptable value of ',i2,' for slot_new; must be > ',i2,' and <= ',i2)
       call TLS_panic (message)
    end if

    ! Bump the starting slot number by 1
    slot = slot + 1

    ! Allocate a mesh size array for each appropriate slot s
    do s = slot,slot_new

       ALLOCATE (ABC(s)%Cell(ncells), STAT = memstat)
       if (memstat /= 0) call TLS_panic ('SLOT_INCREASE: memory allocation error')

       ! Initialize the mesh size array
       call SLOT_SET (ABC, s)

    end do

    slot = slot_new

  END SUBROUTINE SLOT_INCREASE
  
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

  SUBROUTINE SLOT_COMPRESS (ABC, slot)
    !=======================================================================
    ! Purpose(s):
    !
    !   This routine compresses the slot structure to remove unused holes.
    !
    !=======================================================================

    ! Argument List
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: ABC
    integer, intent(INOUT) :: slot

    ! Local Variables
    integer :: s_old, s_new
    logical, dimension(ncells) :: Move_slot
    integer, dimension(ncells) :: Num_mat 

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Zero the starting number of materials in each cell.
    Num_mat = 0

    ! Examine each slot to determine if it contains a material that should
    ! have its storage compressed.
    do s_old = 1,slot

       Move_slot = ABC(s_old)%Cell%Id /= 0

       ! If this slot has a valid material, increment the number of
       ! materials in this cell.
       where (Move_slot) Num_mat = Num_mat + 1

       ! If it does have a material, put its properties in the slot
       ! corresponding to the current number of materials in this cell.
       do s_new = 1,slot
          if (s_new == s_old) cycle
          where (Move_slot .and. s_new == Num_mat)
             ABC(s_new)%Cell%Id       = ABC(s_old)%Cell%Id
             ABC(s_new)%Cell%Vof      = ABC(s_old)%Cell%Vof
             ABC(s_old)%Cell%Id       = 0
             ABC(s_old)%Cell%Vof      = 0.0_r8
          end where
       end do

    end do

  END SUBROUTINE SLOT_COMPRESS

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
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE MATL_ASSIGN (Matl_Dest, Matl_Src)
    !====================================================================
    ! Purposes(s):
    !   Assign matl data structure across =.  This involves
    !   copy of appropriate Cell pointers.
    !
    !   WARNING: This routine assumes that matl data structures have
    !   max_slots slots, of which mat_slot are used.
    !
    !====================================================================

    ! Arguments
    type(MATL_SLOT), dimension(max_slots), intent(IN)  :: Matl_Src
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: Matl_Dest

    ! Local variables
    integer :: s

    do s = 1, mat_slot
       Matl_Dest(s)%Cell = Matl_Src(s)%Cell
    end do

  END SUBROUTINE MATL_ASSIGN

END MODULE MATL_MODULE
