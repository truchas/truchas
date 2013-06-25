MODULE RIEMANN_MODULE
  !=======================================================================
  ! PURPOSE -
  !   Encapsulate all Riemann solution algorithms
  !=======================================================================
  use constants_module
  use cutoffs_module
  use parameter_module
  use kind_module

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
    implicit none

    ! Global scalars & arrays
    character(LEN = *), optional, intent(IN) :: method
    real(KIND = real_kind), dimension(ncells), intent(IN) :: Reference_State, &
                                                             Neighbor_State
    real(KIND = real_kind), dimension(ncells), intent(OUT) :: Interface_State

    ! Local scalars & arrays
    logical(KIND = log_kind) :: specified_method
    character(LEN = 80) :: riemann_method = 'van Leer'
    logical(KIND = log_kind), dimension(ncells) :: Mask1, Mask2

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

    ! If the method has been specified, take it
    specified_method = PRESENT(method)
    if (specified_method) riemann_method = method

    select case (riemann_method)

       ! van Leer's solution; see Siam J. Sci. Stat. Comp. 5:1-20 (1984)
       case ('van Leer')

       ! Upwind Home
       Mask1 = Reference_State > zero .and. &
               Reference_State + Neighbor_State > 1.0e-08
       where (Mask1) Interface_State = Reference_State

       ! Upwind Neighbor
       Mask2 = .not.Mask1 .and. &
               Neighbor_State < zero .and. &
               Reference_State + Neighbor_State < -1.0e-08
       where (Mask2) Interface_State = Neighbor_State

       ! Flow away from interface
       where (.not.Mask1 .and. .not.Mask2) &
          Interface_State = one_half*(Reference_State + Neighbor_State)

       ! Roe's solution; see J. Comput. Phys. 43:573-572 (1981)
       case ('Roe')

          Interface_State = one_half*(Reference_State + Neighbor_State - &
                                      (Neighbor_State - Reference_State)* &
                            SIGN(one, one_half*(Reference_State + Neighbor_State)))

    end select

    return

  END SUBROUTINE RIEMANN
        
END MODULE RIEMANN_MODULE

