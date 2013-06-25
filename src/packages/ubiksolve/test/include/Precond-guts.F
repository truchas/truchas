    ! -*- Mode: f90 -*-

    ! Args:
    real(Ubik_real_type), intent(in), dimension(:) :: b
    type(Ubik_vector_type), intent(inout) :: x_vec
    integer(Ubik_int_type), intent(inout) :: status

    ! Local:
    integer(Ubik_int_type) :: i
    real(Ubik_real_type), dimension(:), pointer :: x

    ! Initialize return status.
    status = 0
  
    ! Alias for solution vector values.
    x => Ubik_values_ptr(x_vec)

    ! Apply preconditioning as appropriate.
    select case (precond_control%precond)
    case (Ubik_precond_none)

       ! No preconditioning.
       x = b
    case (Ubik_method_Jacobi)

       ! m-step Jacobi preconditioning with relaxation.
       x = 0.0
       do i=1,precond_control%steps
          call Ubik_Jacobi (precond_control%omega, b, coeff, x)
       end do
    case (Ubik_method_SSOR)

       ! m-step SSOR preconditioning.
       x = 0.0
       do i=1,precond_control%steps
          call Ubik_SOR (Ubik_method_SSOR, precond_control%omega, b, coeff, x)
       end do
    case (Ubik_method_IC)

       ! Ap is an incomplete factorization, stored as a lower
       ! triangular matrix and an upper triangular matrix.
       x = b
       call Ubik_SolveLTriang (unit_diagonal=.FALSE., use_transpose=.FALSE., a=precond, x=x)
!!$       if (precond_control%outmode >= 6) then
!!$          call Ubik_Write (precond_control%luout, x, 'Ubik_ApplyPrecond: AFTER FORWARD SWEEP (x):')
!!$          call Ubik_FlushUnit (precond_control%luout, status)
!!$          if (status /= 0) then
!!$             status = -3
!!$             go to 9999
!!$          end if
!!$       end if
       call Ubik_SolveUTriang (unit_diagonal=.FALSE., use_transpose=.TRUE., a=precond, x=x)
    case (Ubik_method_ILU)

       ! Ap is an incomplete factorization, stored as a unit
       ! lower triangular matrix and an upper triangular matrix.
       x = b
       call Ubik_SolveLTriang (unit_diagonal=.TRUE., use_transpose=.FALSE., a=precond, x=x)
!!$       if (precond_control%outmode >= 6) then
!!$          call Ubik_Write (precond_control%luout, x, 'Ubik_ApplyPrecond: AFTER FORWARD SWEEP (x):')
!!$          call Ubik_FlushUnit (precond_control%luout, status)
!!$          if (status /= 0) then
!!$             status = -3
!!$             go to 9999
!!$          end if
!!$       end if
       call Ubik_SolveUTriang (unit_diagonal=.FALSE., use_transpose=.FALSE., a=precond, x=x)
!!$    case (Ubik_method_approxinv)
!!$
!!$       ! Ap contains an approximate inverse.
!!$       x = MatMul(precond, b)
!!$    case (Ubik_method_diag)
!!$
!!$       ! Array precond contains diag(coeff).
!!$       x = b / Ubik_diag(precond)
!!$    case (Ubik_method_invdiag)
!!$
!!$       ! Array precond contains inv(diag(coeff)).
!!$       x = b*Ubik_diag(precond)
    case default

       ! Invalid value for precond_control%precond.
       status = Ubik_invalid_input__
    end select

9999 continue
