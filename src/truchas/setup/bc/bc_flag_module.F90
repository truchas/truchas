!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BC_FLAG_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define procedures which set boundary condition flag bits.
  !
  ! Contains: ASSIGN_BC_BITS
  !           SET_DIRICHLET
  !           SET_FREE_SLIP
  !           SET_DIRICHLET_VEL
  !           SET_INTERNAL_BC
  !           SET_NEUMANN
  !           SET_NEUMANN_VEL
  !           SET_VELOCITY_BC
  !           SET_NO_VEL_BC
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !            Bryan R. Lally, LANL ESA-EPE (lally@lanl.gov)
  !
  !=======================================================================
  use legacy_mesh_api, only: ncells, nfc
  implicit none

  ! Private Module
  private

  ! Public Variables

  ! Public Subroutines
  public :: ASSIGN_BC_BITS, SET_DIRICHLET, SET_FREE_SLIP,    &
            SET_DIRICHLET_VEL, SET_INTERNAL_BC, SET_NEUMANN, &
            SET_NEUMANN_VEL, SET_NO_VEL_BC !, SET_VELOCITY_BC

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE ASSIGN_BC_BITS (Prs_bit, Vel_bit, Conc_bit)
    !=======================================================================
    ! Purpose(s):
    !   Assign bit positions in the Flag portion of the BC structure
    !   (BC%Flag) for pressure, velocity, and concentration BC flags.
    !=======================================================================

    ! Argument List
    integer, dimension(nfc) :: Conc_bit
    integer, dimension(nfc) :: Prs_bit
    integer, dimension(nfc) :: Vel_bit

    ! Local Variables
    integer :: f

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! The following convention holds for the bits in BC%Flag:
    !
    !      Bit     BC Variable      Face                   
    !      ---     -----------      ----
    !      6-11     Pressure         1-6 (1 bit per face)
    !     12-23     Velocity         1-6 (2 bits per face)
    !     24-29   Concentration      1-6 (1 bit per face)
    !
    !      BC Variable   Bit Value      BC Type
    !      -----------   ---------      -------
    !        Pressure        0          Neumann   (Default)
    !        Pressure        1         Dirichlet
    !        Velocity       00           None     (Internal Face Default)
    !        Velocity       01         Free-Slip  (Mesh Boundary Face Default)
    !        Velocity       10         Dirichlet
    !        Velocity       11          Neumann
    !      Concentration     0          Neumann   (Default)
    !      Concentration     1         Dirichlet

    ! Assign Bit Locations
    do f = 1,nfc
       Conc_bit(f) = 23 + f
       Prs_bit(f)  = f + nfc - 1
       Vel_bit(f)  = 2*(f + nfc - 1)
    end do

  END SUBROUTINE ASSIGN_BC_BITS
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_DIRICHLET (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the Dirichlet bit in integer Flag where Mask is true.
    !   One bit (in bit_position) is set to one. This routine
    !   should only be called for temperature and pressure BC.
    !=======================================================================

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer,                    intent(IN)    :: bit_position
    integer, dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) Flag = IBSET(Flag, bit_position)

  END SUBROUTINE SET_DIRICHLET

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_FREE_SLIP (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the Free-Slip bit in integer Flag where Mask is true.
    !   Sets bits to '01'.
    !=======================================================================

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer,                    intent(IN)    :: bit_position
    integer, dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) 
       Flag = IBCLR(Flag, bit_position)
       Flag = IBSET(Flag, bit_position + 1)
    end where

  END SUBROUTINE SET_FREE_SLIP

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_DIRICHLET_VEL (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the velocity Dirichlet bit in integer Flag where Mask is true.
    !   Sets bits to '10'.
    !=======================================================================

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer,                    intent(IN)    :: bit_position
    integer, dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) 
       Flag = IBSET(Flag, bit_position)
       Flag = IBCLR(Flag, bit_position + 1)
    end where

  END SUBROUTINE SET_DIRICHLET_VEL

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_INTERNAL_BC (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the internal BC bit in integer Flag where Mask is true.
    !   One bit (in bit_position) is set to one.
    !=======================================================================

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer,                    intent(IN)    :: bit_position
    integer, dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) Flag = IBSET(Flag, bit_position)

  END SUBROUTINE SET_INTERNAL_BC

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_NEUMANN (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the Neumann bit in integer Flag where Mask is true.
    !   One bit (in bit_position) is set to zero. This routine
    !   should only be called for temperature and pressure BC.
    !=======================================================================

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer,                    intent(IN)    :: bit_position
    integer, dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) Flag = IBCLR(Flag, bit_position)

  END SUBROUTINE SET_NEUMANN
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_NEUMANN_VEL (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the velocity Neumann bit in integer Flag where Mask is true.
    !   Sets bits to '11'.
    !=======================================================================

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer,                    intent(IN)    :: bit_position
    integer, dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) 
       Flag = IBSET(Flag, bit_position)
       Flag = IBSET(Flag, bit_position + 1)
    end where

  END SUBROUTINE SET_NEUMANN_VEL

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!! NNC, Jan 2014.  Time-dependent dirichlet velocity.  This routine is unused,
!! and the change below to use time_step_module to get the time introduced
!! a circular module dependency.  Hence I'm just commenting it out instead.

!  SUBROUTINE SET_VELOCITY_BC (Velocity, bit_position, f)
!    !=======================================================================
!    ! Purpose(s):
!    !   Apply velocity boundary conditions to the velocity vector 
!    !   (U,V,W) lying on face f. BCs are applied according to the 
!    !   bit positions in BC%Flag.
!    !=======================================================================
!    use bc_data_module,   only: BC_Vel, BC
!    use bc_kind_module,   only: FREE_SLIP, DIRICHLET_VEL
!    use kinds, only: r8
!    use legacy_mesh_api, only: ncells, ndim, Cell
!    use time_step_module, only: t
!
!    ! Argument List
!
!    integer, intent(IN) :: bit_position, f
!    real(r8), dimension(ndim,ncells), intent(INOUT) :: Velocity
!
!    ! Local Variables
!
!    logical, dimension(ncells) :: Mask
!    integer :: n, j
!    real(r8), dimension(ncells) :: Mag
!
!    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!    ! FREE_SLIP velocity BC
!    Mask = FREE_SLIP (BC%Flag, bit_position)
!    ! Compute the projection of the velocity and the normal
!    if (ANY(Mask)) then
!       Mag = 0.0_r8
!       do n = 1,ndim
!          where (Mask) Mag = Mag + Velocity(n,:)*Cell%Face_Normal(n,f)
!       end do
!       do n = 1,ndim
!          where (Mask) Velocity(n,:) = Velocity(n,:) - Mag*Cell%Face_Normal(n,f)
!       end do
!    end if
!
!    ! DIRICHLET velocity BC
!    Mask = DIRICHLET_VEL (BC%Flag, bit_position)
!    !! NNC, Jan 2014.  Time-dependent dirichlet velocity.
!    !ORIG: if (ANY(Mask)) then
!    !ORIG:    do n = 1,ndim
!    !ORIG:       where (Mask) Velocity(n,:) = BC_Vel(n,f,:)
!    !ORIG:    end do
!    !ORIG: end if
!    do j = 1, ncells
!       if (Mask(j)) Velocity(:,j) = bndry_vel%get(f,j,t)
!    end do
!
!  END SUBROUTINE SET_VELOCITY_BC
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_NO_VEL_BC (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the velocity BC bit to '00', indicating that no velocity BC is
    !   specified.
    !=======================================================================

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer,                    intent(IN)    :: bit_position
    integer, dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) 
       Flag = IBCLR(Flag, bit_position)
       Flag = IBCLR(Flag, bit_position + 1)
    end where

  END SUBROUTINE SET_NO_VEL_BC

END MODULE BC_FLAG_MODULE
