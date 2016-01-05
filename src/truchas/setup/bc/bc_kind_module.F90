!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BC_KIND_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define procedures for the various kinds of boundary conditions.
  !
  ! Contains: DIRICHLET
  !           NEUMANN
  !           FREE_SLIP
  !           DIRICHLET_VEL
  !           NEUMANN_VEL
  !           INTERNAL_BC
  !           IN_FLOW
  !           OUT_FLOW
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
  public :: DIRICHLET, NEUMANN, FREE_SLIP, DIRICHLET_VEL, NEUMANN_VEL, &
            INTERNAL_BC, IN_FLOW, OUT_FLOW

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  FUNCTION DIRICHLET (Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Return true if the Dirichlet bit (in bit_position) in integer
    !   Flag is set; otherwise return false.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    integer, intent(IN) :: bit_position
    integer, intent(IN) :: Flag(:)

    ! Local Variables

    ! Function Return
    logical, dimension(SIZE(Flag)) :: Dirichlet

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    Dirichlet = BTEST(Flag, bit_position)

    return

  END FUNCTION DIRICHLET
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION NEUMANN (Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Return true if the Neumann bit (in bit_position) in integer
    !   Flag is set; otherwise return false.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    integer, intent(IN) :: bit_position
    integer, intent(IN) :: Flag(:)

    ! Local Variables

    ! Function Return
    logical, dimension(SIZE(Flag)) :: Neumann

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    Neumann = .not.(BTEST(Flag, bit_position))

    return

  END FUNCTION NEUMANN
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION FREE_SLIP (Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Return true if the Free-slip bits (in bit_position and
    !   bit_position+1) in integer Flag are set; otherwise return false.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    integer, intent(IN) :: bit_position
    integer, intent(IN) :: Flag(:)

    ! Local Variables

    ! Function Return
    logical, dimension(SIZE(Flag)) :: Free_slip

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    Free_slip = .not.(BTEST(Flag, bit_position)) .and. &
                      BTEST(Flag, bit_position+1)

    return

  END FUNCTION FREE_SLIP
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION DIRICHLET_VEL (Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Return true if the Dirichlet_Vel bits (in bit_position and
    !   bit_position+1) in integer Flag are set; otherwise return false.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    integer, intent(IN) :: bit_position
    integer, intent(IN) :: Flag(:)

    ! Local Variables

    ! Function Return
    logical, dimension(SIZE(Flag)) :: Dirichlet_Vel

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    Dirichlet_Vel =       BTEST(Flag, bit_position) .and. &
                    .not.(BTEST(Flag, bit_position+1))

    return

  END FUNCTION DIRICHLET_VEL
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION NEUMANN_VEL (Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Return true if the Neumann_Vel bits (in bit_position and
    !   bit_position+1) in integer Flag are set; otherwise return false.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    integer, intent(IN) :: bit_position
    integer, intent(IN) :: Flag(:)

    ! Local Variables

    ! Function Return
    logical, dimension(SIZE(Flag)) :: Neumann_Vel

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    Neumann_Vel = BTEST(Flag, bit_position) .and. &
                  BTEST(Flag, bit_position+1)

    return

  END FUNCTION NEUMANN_VEL
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION INTERNAL_BC (Flag, bit_position)
    !=======================================================================
    ! Purpose(s):
    !   Return true if the internal BC bit (in bit_position) in integer
    !   Flag is set; otherwise return false.
    !=======================================================================
    use scalars_module

    implicit none

    ! Argument List
    integer, intent(IN) :: bit_position
    integer, intent(IN) :: Flag(:)

    ! Local Variables

    ! Function Return
    logical, dimension(SIZE(Flag)) :: Internal_BC

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Function Return
    Internal_BC = BTEST(Flag, bit_position)

    return

  END FUNCTION INTERNAL_BC
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION IN_FLOW (f, Fluxing_Velocity)
    !=======================================================================
    ! Purpose(s):
    !   Return true for those faces that have negative values of the
    !   Fluxing_Velocity.
    !=======================================================================
    use bc_data_module, only: BC
    use bc_type_module, only: Vel
    use mesh_module,    only: Mesh
    use scalars_module
    use kinds, only: r8

    ! Argument List

    integer :: f
    real(r8), dimension(nfc,ncells), intent(IN) :: Fluxing_Velocity

    ! Local Variables

    real(r8), dimension(ncells) :: FV
    logical, dimension(ncells) :: Mask

    ! Function Return

    logical, dimension(ncells) :: In_flow

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Mask(:) = Mesh(:)%Ngbr_Face(f) == 0
    Mask = Mask .and. .not. FREE_SLIP (BC%Flag, Vel%Face_bit(f))

    ! set up a local FV face velocity array
    FV(:) = Fluxing_Velocity(f,:)

    ! FV < 0 indicates an inflow face
    In_flow = .false.
    where (Mask .and. FV(:) < 0.0_r8) In_flow(:) = .true.

  END FUNCTION IN_FLOW
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  FUNCTION OUT_FLOW (f, Fluxing_Velocity)
    !=======================================================================
    ! Purpose(s):
    !   Return true for those faces that have a positive values of the
    !   Fluxing_Velocity.
    !=======================================================================
    use bc_data_module, only: BC
    use bc_type_module, only: Vel
    use mesh_module,    only: Mesh
    use scalars_module
    use kinds, only: r8

    ! Argument List

    integer :: f
    real(r8), dimension(nfc,ncells), intent(IN) :: Fluxing_Velocity

    ! Local Variables

    real(r8), dimension(ncells) :: FV
    logical, dimension(ncells) :: Mask

    ! Function Return

    logical, dimension(ncells) :: Out_flow

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Mask(:) = Mesh(:)%Ngbr_Face(f) == 0
    Mask = Mask .and. .not. FREE_SLIP (BC%Flag, Vel%Face_bit(f))

    ! set up a local FV face velocity array
    FV(:) = Fluxing_Velocity(f,:)

    ! FV > 0 indicates an outflow face
    Out_flow = .false.
    where (Mask .and. FV(:) > 0.0_r8) Out_flow(:) = .true.

  END FUNCTION OUT_FLOW

END MODULE BC_KIND_MODULE
