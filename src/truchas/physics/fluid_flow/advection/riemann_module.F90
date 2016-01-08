!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE RIEMANN_MODULE
  !=======================================================================
  ! PURPOSE -
  !   Encapsulate all Riemann solution algorithms
  !=======================================================================
  use kinds, only: r8
  use cutoffs_module
  use parameter_module
  implicit none
  save

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE RIEMANN (Reference_State, Neighbor_State, Interface_State, &
                      method)
    !=======================================================================
    ! PURPOSE -
    !      Solve a Riemann problem based on a solution of the 
    !      inviscid Burger's equation in one dimension
    !=======================================================================
    ! Global scalars & arrays
    character(*), optional, intent(IN) :: method
    real(r8), dimension(ncells), intent(IN) :: Reference_State, &
                                                             Neighbor_State
    real(r8), dimension(ncells), intent(OUT) :: Interface_State

    ! Local scalars & arrays
    logical :: specified_method
    character(80) :: riemann_method = 'van Leer'
    logical, dimension(ncells) :: Mask1, Mask2

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

    ! If the method has been specified, take it
    specified_method = PRESENT(method)
    if (specified_method) riemann_method = method

    select case (riemann_method)

       ! van Leer's solution; see Siam J. Sci. Stat. Comp. 5:1-20 (1984)
       case ('van Leer')

       ! Upwind Home
       Mask1 = Reference_State > 0.0_r8 .and. &
               Reference_State + Neighbor_State > 1.0e-08
       where (Mask1) Interface_State = Reference_State

       ! Upwind Neighbor
       Mask2 = .not.Mask1 .and. &
               Neighbor_State < 0.0_r8 .and. &
               Reference_State + Neighbor_State < -1.0e-08
       where (Mask2) Interface_State = Neighbor_State

       ! Flow away from interface
       where (.not.Mask1 .and. .not.Mask2) &
          Interface_State = 0.5_r8*(Reference_State + Neighbor_State)

       ! Roe's solution; see J. Comput. Phys. 43:573-572 (1981)
       case ('Roe')

          Interface_State = 0.5_r8*(Reference_State + Neighbor_State - &
                                      (Neighbor_State - Reference_State)* &
                            SIGN(1.0_r8, 0.5_r8*(Reference_State + Neighbor_State)))

    end select

  END SUBROUTINE RIEMANN
        
END MODULE RIEMANN_MODULE

