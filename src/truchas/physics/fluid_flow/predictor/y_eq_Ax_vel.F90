MODULE Y_EQ_AX_VEL
  !=======================================================================
  ! Purpose:
  !
  !   Define procedures necessary to perform a matrix-vector multiply
  !   (Ax), returning it in vector y
  !
  !   Public Interface:
  !
  !     * call Y_EQ_AX_VELOCITY (X, Y, status)
  !
  !         Returns in y the vector Ax, where A represents a discrete
  !         operator for the viscous terms of the momentum equations
  !
  ! Contains: Y_EQ_AX_VELOCITY
  !
  ! Author(s): Telluride Team
  !            M. A. Christon (christon@lanl.gov)
  !
  !=======================================================================

  Implicit None
  Private

  ! public procedures
  Public :: Y_EQ_AX_VELOCITY

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE Y_EQ_AX_VELOCITY (X_vec, Y, status)
    !=======================================================================
    ! Purpose:
    !
    !=======================================================================

    use bc_operations
    use constants_module,     only: zero
    use fluid_data_module,    only: fluidVof, fluidRho, &
                                    Drag_Coefficient, &
                                    momentum_solidify_implicitness, &
                                    IsPureImmobile
    use flow_phase_change,    only: have_solidifying_flow, solidified_rho
    use parameter_module,     only: ncells, ndim

    use time_step_module,     only: dt

    use kind_module,          only: int_kind, real_kind
    use porous_drag_data,     only: porous_flow, porous_implicitness
    use timing_tree
    use viscous_data_module,  only: viscous_implicitness
 
    use viscous_module,       only: STRESS_GRADIENT

    use UbikSolve_module

    implicit none

    ! Arguments
    type (Ubik_vector_type),                intent(INOUT) :: X_vec
    real (real_kind), dimension(:), target, intent(INOUT) :: Y
    integer (int_kind),                     intent(OUT)   :: status

    ! Local Variables
    real(KIND = real_kind), dimension(ndim,ncells) :: Stress_Grad
    real(KIND = real_kind), dimension(ndim,ncells) :: Vc

    integer(KIND = int_kind)                                :: i, j, n
    real(real_kind), dimension(:), pointer :: X

    ! real(real_kind), dimension(:), pointer :: Stress_Grad_Pointer
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    X => Ubik_values_ptr(X_vec)

    ! Start the timer.
    call start_timer("Timer_Solver_TMP2")

    ! For now just copy X into Vc(ndim,ncells)
    do i = 1, ncells
       do n = 1, ndim
           ! where in X is vel component n for cell i
           j = (i-1)*(ndim)+n
           Vc(n,i) = X(j)
       end do
    end do
    
    ! initialize the Y
    Y = zero

    call STRESS_GRADIENT(Stress_Grad, Vc)

    ! Old Methodology for the mat-vec -- replaced 9/17/07 (MAC)
    ! For now just copy Stress_Grad into Y
    ! do n = 1, ndim
    !    do i = 1, ncells
    !       ! where in Y is vel component n for cell i
    !       j = (i-1)*(ndim)+n
    !       ! Y(J) include fluidVof on porous drag terms for consistency (MAC)
    !       if(FluidVof(i).gt.zero) then
    !          Y(j) = ( (fluidRho(i)*fluidVof(i) &
    !              + dt* porous_implicitness*fluidVof(i)*Drag_Coefficient(n,i))*Vc(n,i)&
    !              - dt*viscous_implicitness*fluidVof(i)*(Stress_Grad(n,i) &
    !              - Stress_Grad_BC(n,i)) )/fluidVof(i)
    !       else
    !          Y(J) = Vc(n,i)
    !       endif
    !    end do
    ! end do
    
    ! New methodology for the mat-vec
    ! Setup the base terms -- intertial and viscous stress
    do i = 1, ncells
       do n = 1, ndim
          ! where in Y is vel component n for cell i
          j = (i-1)*(ndim)+n
          ! The BC component of the stress gradient is not removed here.  See
          ! Viscous StressDescription in the flow modeling docs on source forge.
          Y(j) = fluidRho(i)*Vc(n,i) -  &
                 dt*viscous_implicitness*Stress_Grad(n,i)
       end do
    end do
    ! Handle porous drag if required
    if (porous_flow) then
       do i = 1, ncells
          do n = 1, ndim
             ! where in Y is vel component n for cell i
             j = (i-1)*(ndim)+n
             Y(j) = Y(j) + &
                  dt*porous_implicitness*Drag_Coefficient(n,i)*Vc(n,i)
          end do
       end do
    endif
    ! Handle the phase change terms -- need time-weighting done correctly
    if (have_solidifying_flow()) then
       do i = 1, ncells
          do n = 1, ndim
             ! where in Y is vel component n for cell i
             j = (i-1)*(ndim)+n
             if (fluidVof(i) > zero) then
                Y(j) = Y(j) + &
                momentum_solidify_implicitness*solidified_rho(i)*Vc(n,i)/fluidVof(i)
             endif
          end do
       end do
    endif
    ! Now mask the void cells
    do i = 1, ncells
       do n = 1, ndim
          ! where in Y is vel component n for cell i
          j = (i-1)*(ndim)+n
          if (IsPureImmobile(i) .or. fluidRho(i) == zero) Y(j) = Vc(n,i)
       end do
    end do
    
    call stop_timer("Timer_Solver_TMP2")
    
    status = 0
    return

  END SUBROUTINE Y_EQ_AX_VELOCITY

END MODULE Y_EQ_AX_VEL
