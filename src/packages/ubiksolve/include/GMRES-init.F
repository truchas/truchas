! -*- Mode: f90 -*-

    ! reset return status
    call set_status_converged (control)

    ! local alias for number of unknowns
    nunk = SIZE(x)

    ! check validity of stopping test - modify if necessary
    call check_stopping_nonstationary (control)

    ! initialize control params
    call initialize (control)

    ! shorthand
    luerr_ = luerr(control)
    luout_ = luout(control)
    norm_type_ = norm_type(control)
    outmode_iterates_ = outmode_iterates(control)
    outmode_full_ = outmode_full(control)

    ! check to see if b is essentially zero
    !
    ! this is dangerous since we're not scaling by anything, but
    ! the default value of smallest(control) is TINY(1.0), or
    ! typically about 1e-38, so the hope is that if the 1-norm of
    ! b (the sum of abs. values of the elements of the rhs) is less
    ! than that, it really is essentially zero, so returning the
    ! trivial solution is reasonable
    call set_bnorm (control, Norm(1, b))
    if (bnorm(control) < smallest(control)) then

       ! b is essentially zero, so declare victory
       x = 0.0

       if (outmode_warnings(control)) then
          call Spew (luerr_, ' ')
          call Spew (luerr_, '*****************************************')
          call Spew (luerr_, TRIM(routine)//': WARNING')
          call Spew (luerr_, '  it appears that the right-hand side (rhs) of the system')
          call Spew (luerr_, '  is essentially zero:')
          call Spew (luerr_, ' ')
          write(out_string,*) '    ||b||_1 =', bnorm(control)
          call Spew (luerr_, out_string)
          call Spew (luerr_, ' ')
          call Spew (luerr_, '  where ||b||_1 is the 1-norm of the rhs, or the sum of')
          call Spew (luerr_, '  the absolute values of the elements of the rhs')
          call Spew (luerr_, ' ')
          call Spew (luerr_, '  if you believe the system should still be valid, try')
          call Spew (luerr_, '  scaling the system before attempting to solve')
          call Spew (luerr_, '*****************************************')
       end if
       go to 999
    else if (need_bnorm(control) .and. norm_type_ /= 1) then

       ! re-compute || b || if using Linf or L2 norm
       call set_bnorm (control, Norm(norm_type_, b))
    end if

    ! flag indicating whether current iterate is used in stopping test
    stopping_uses_x = stopping_axb(control) .or. stopping_x(control)

    ! flag indicating whether updating of iterate is deferred or performed
    ! each iteration
    update_iterate = &
         outmode_iterates(control) .or. &
         residual_compute_only(control) .or. &
         stopping_true(control) .or. &
         stopping_uses_x

    ! additional input checking
    if (subspace(control) <= 0) then
       if (outmode_warnings(control)) then
          call Spew (luerr_, ' ')
          call Spew (luerr_, '*****************************************')
          call Spew (luerr_, TRIM(routine)//': WARNING')
          write(out_string,*) '  Invalid value for subspace: ', subspace(control)
          call Spew (luerr_, out_string)
          call Spew (luerr_, ' ')
          call Spew (luerr_, '  (value must be >0)')
          call Spew (luerr_, 'Continuing with subspace=20')
          call Spew (luerr_, '*****************************************')
       end if
       call set_subspace (control, 20)
    end if

    ! flag indicating whether preconditioning is performed
    precond = PRESENT(ApplyPrecond)

    ! let m be a synonym for subspace (to match other shorthand vars for
    ! components of control_type it should be subspace_ or at least m_;
    ! use of m is historical)
    m = subspace(control)

    ! Allocate local arrays.
    allocate (w(nunk), xt(nunk), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure (control)
       go to 888
    end if
    allocate (c(m+1), s(m+1), r(m+1), rt(m+1), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure (control)
       go to 888
    end if
    allocate (h(m+1,m+1), v(nunk,m+1), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure (control)
       go to 888
    end if

    ! Initialize local arrays.
    w = 0.0
    xt = 0.0
    c = 0.0
    s = 0.0
    r = 0.0
    rt = 0.0
    h = 0.0
    v = 0.0

    ! Initialize vector derived types.
    call set_values_ptr (x_vec, x)
    call set_values_ptr (xt_vec, xt)

    ! Write initial system if desired.
    if (outmode_system(control)) then
       call spew (luout_, x, b, TRIM(routine)//': ON ENTRY (x,b):')
    end if

    ! Initialize xold if needed.
    if (need_xold(control)) then
       call set_xold (control, x)
    end if

    ! Copy initial guess into xt if needed for stopping test.
    if (stopping_uses_x) then
       xt = x
    end if
