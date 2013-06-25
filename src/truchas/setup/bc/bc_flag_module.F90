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
  implicit none

  ! Private Module
  private

  ! Public Variables

  ! Public Subroutines
  public :: ASSIGN_BC_BITS, SET_DIRICHLET, SET_FREE_SLIP,    &
            SET_DIRICHLET_VEL, SET_INTERNAL_BC, SET_NEUMANN, &
            SET_NEUMANN_VEL, SET_VELOCITY_BC, SET_NO_VEL_BC

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE ASSIGN_BC_BITS (Prs_bit, Vel_bit, Conc_bit)
    !=======================================================================
    ! Purpose(s):
    !   Assign bit positions in the Flag portion of the BC structure
    !   (BC%Flag) for pressure, velocity, and concentration BC flags.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    integer(KIND = int_kind), dimension(nfc) :: Conc_bit
    integer(KIND = int_kind), dimension(nfc) :: Prs_bit
    integer(KIND = int_kind), dimension(nfc) :: Vel_bit

    ! Local Variables
    integer(KIND = int_kind) :: f

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

    return

  END SUBROUTINE ASSIGN_BC_BITS
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_DIRICHLET (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the Dirichlet bit in integer Flag where Mask is true.
    !   One bit (in bit_position) is set to one. This routine
    !   should only be called for temperature and pressure BC.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer(KIND = int_kind),                    intent(IN)    :: bit_position
    integer(KIND = int_kind), dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) Flag = IBSET(Flag, bit_position)

    return

  END SUBROUTINE SET_DIRICHLET

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_FREE_SLIP (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the Free-Slip bit in integer Flag where Mask is true.
    !   Sets bits to '01'.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer(KIND = int_kind),                    intent(IN)    :: bit_position
    integer(KIND = int_kind), dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) 
       Flag = IBCLR(Flag, bit_position)
       Flag = IBSET(Flag, bit_position + 1)
    end where

    return

  END SUBROUTINE SET_FREE_SLIP

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_DIRICHLET_VEL (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the velocity Dirichlet bit in integer Flag where Mask is true.
    !   Sets bits to '10'.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer(KIND = int_kind),                    intent(IN)    :: bit_position
    integer(KIND = int_kind), dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) 
       Flag = IBSET(Flag, bit_position)
       Flag = IBCLR(Flag, bit_position + 1)
    end where

    return

  END SUBROUTINE SET_DIRICHLET_VEL

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_INTERNAL_BC (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the internal BC bit in integer Flag where Mask is true.
    !   One bit (in bit_position) is set to one.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer(KIND = int_kind),                    intent(IN)    :: bit_position
    integer(KIND = int_kind), dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) Flag = IBSET(Flag, bit_position)

    return

  END SUBROUTINE SET_INTERNAL_BC

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_NEUMANN (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the Neumann bit in integer Flag where Mask is true.
    !   One bit (in bit_position) is set to zero. This routine
    !   should only be called for temperature and pressure BC.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer(KIND = int_kind),                    intent(IN)    :: bit_position
    integer(KIND = int_kind), dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) Flag = IBCLR(Flag, bit_position)

    return

  END SUBROUTINE SET_NEUMANN
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_NEUMANN_VEL (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the velocity Neumann bit in integer Flag where Mask is true.
    !   Sets bits to '11'.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer(KIND = int_kind),                    intent(IN)    :: bit_position
    integer(KIND = int_kind), dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) 
       Flag = IBSET(Flag, bit_position)
       Flag = IBSET(Flag, bit_position + 1)
    end where

    return

  END SUBROUTINE SET_NEUMANN_VEL

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_VELOCITY_BC (Velocity, bit_position, f)
    !=======================================================================
    ! Purpose(s):
    !   Apply velocity boundary conditions to the velocity vector 
    !   (U,V,W) lying on face f. BCs are applied according to the 
    !   bit positions in BC%Flag.
    !=======================================================================
    use bc_data_module,   only: BC_Vel, BC
    use bc_kind_module,   only: FREE_SLIP, DIRICHLET_VEL
    use constants_module, only: zero
    use kind_module,      only: int_kind, log_kind, real_kind
    use mesh_module,      only: Cell
    use parameter_module, only: ncells, ndim

    implicit none

    ! Argument List

    integer(int_kind),                       intent(IN)    :: bit_position, f
    real(real_kind), dimension(ndim,ncells), intent(INOUT) :: Velocity

    ! Local Variables

    logical(log_kind), dimension(ncells) :: Mask
    integer(int_kind)                    :: n
    real(real_kind),   dimension(ncells) :: Mag

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! FREE_SLIP velocity BC
    Mask = FREE_SLIP (BC%Flag, bit_position)
    ! Compute the projection of the velocity and the normal
    if (ANY(Mask)) then
       Mag = zero
       do n = 1,ndim
          where (Mask) Mag = Mag + Velocity(n,:)*Cell%Face_Normal(n,f)
       end do
       do n = 1,ndim
          where (Mask) Velocity(n,:) = Velocity(n,:) - Mag*Cell%Face_Normal(n,f)
       end do
    end if

    ! DIRICHLET velocity BC
    Mask = DIRICHLET_VEL (BC%Flag, bit_position)
    if (ANY(Mask)) then
       do n = 1,ndim
          where (Mask) Velocity(n,:) = BC_Vel(n,f,:)
       end do
    end if

    return

  END SUBROUTINE SET_VELOCITY_BC
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SET_NO_VEL_BC (Mask, Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Set the velocity BC bit to '00', indicating that no velocity BC is
    !   specified.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    logical, dimension(ncells), intent(IN) :: Mask

    integer(KIND = int_kind),                    intent(IN)    :: bit_position
    integer(KIND = int_kind), dimension(ncells), intent(INOUT) :: Flag

    ! Local Variables

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set Flag
    where (Mask) 
       Flag = IBCLR(Flag, bit_position)
       Flag = IBCLR(Flag, bit_position + 1)
    end where

    return

  END SUBROUTINE SET_NO_VEL_BC

END MODULE BC_FLAG_MODULE
