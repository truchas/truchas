! -*- Mode: f90 -*-

    type(precond_control_type), intent(in) :: control

    ! Construct preconditioner as appropriate.
    select case (control%precond)
    case (Ubik_method_IC)

       ! Incomplete Cholesky factorization.
       precond = coeff
       call Ubik_IC (control%omega, precond)
    case (Ubik_method_ILU)

       ! Incomplete LU factorization (store unit lower triangular matrix
       ! and upper triangular matrix in ap).
       precond = coeff
       call Ubik_ILU (control%omega, precond)
!!$    case (Ubik_method_approxinv)
!!$
!!$       ! Approximate the inverse of a by inv(diag(a)).
!!$       precond = 0.0
!!$       Ubik_Diag(precond) = 1.0 / Ubik_Diag(coeff)
    end select

!!$    if (control%outmode >= 4) then
!!$       call Ubik_Write (control%luout, precond, 'Ubik_GeneratePreconditioner: PRECONDITIONER:')
!!$    end if

    return
