program test_alloy_back_diff_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use parameter_list_json
  use alloy_back_diff_type
  
  integer :: inlun, stat, n
  character(128) :: arg
  character(:), allocatable :: infile
  type(parameter_list), pointer :: params, plist
  character(:), allocatable :: errmsg
  type(alloy_back_diff) :: pd !

  call get_command_argument(1, arg)
  infile = trim(arg)

  open(newunit=inlun,file=infile,action='read',access='stream')
  call parameter_list_from_json_stream(inlun, params, errmsg)
  if (.not.associated(params)) then
    write(*,*) errmsg
    stop 1
  end if
  close(inlun)

  block
    use material_factory, only: alloc_material
    use material_class
    class(material), allocatable :: matl
    plist => params%sublist('alloy')
    call alloc_material(matl, 'alloy', plist, stat, errmsg)
    if (stat /= 0) then
      write(*,*) errmsg
      stop 2
    end if
    call matl%add_enthalpy_prop(stat, errmsg)
    if (stat /= 0) then
      write(*,*) errmsg
      stop
    end if
    plist => params%sublist('phase-diagram')
    call pd%init(matl, plist, stat, errmsg)
    if (stat /= 0) then
      write(*,*) errmsg
      stop
    end if
  end block

  call timestep

contains

  subroutine timestep
  
    use state_history_type
    use nka_type
    use block_linear_solver
    use iso_fortran_env

    real(r8) :: T0, T1, Hdot, t, dt, tol, error
    real(r8), allocatable :: u0(:), u(:), udot(:), du(:), jac(:,:)
    real(r8), allocatable :: C0(:), C_sol(:), C_liq(:), C_liq_max(:)
    type(state_history) :: uhist
    type(nka) :: accel
    integer :: iter, max_iter
    type(alloy_back_diff_jac) :: newjac
    real(r8) :: dfTdH, dFTdT
    real(r8), allocatable :: dfpdp(:), dfpdg(:), dfgdp(:), Cdot(:)

    n = pd%num_comp + 3
    allocate(u0(n), u(n), udot(n), du(n), jac(n,n))
    n = pd%num_comp
    allocate(Cdot(n), C_sol(n), C_liq(n), C_liq_max(n))
    allocate(dfpdp(n), dfpdg(n), dfgdp(n))

    call newjac%init(pd%num_comp)

    call params%get('C0', C0)
    call params%get('T0', T0)
    call params%get('T1', T1)
    call params%get('Hdot', Hdot)

    !! Starting state (pure liquid
    u(1:n) = C0
    u(n+1) = 1
    u(n+2) = pd%H_liq%eval([T0])
    u(n+3) = T0

    udot(1:n) = 0
    udot(n+1) = 0
    udot(n+2) = Hdot
    udot(n+3) = Hdot ! assumes unit heat capacity (per unit volume)
    
    t = 0
    call uhist%init(2, t, u, udot) ! maintain 2 solution vectors
    call accel%init(size(u), 4)

    call params%get('dt', dt)
    call params%get('tol', tol)
    call params%get('max-iter', max_iter)

    C_liq_max = 0
    Cdot = 0
    
    do while (u(n+3) > T1)
    
      t = t+dt
      
      u0 = u
      call uhist%interp_state(t, u, order=1) ! linear extrapolation
      
      block ! truchas-like
        call pd%compute_f_jac(C0, Cdot, u(1:n), u(n+1), u(n+2), u(n+3), udot(1:n), udot(n+1), dt, newjac)
        call newjac%lu_factor
        dfTdH = 1/dt
        dfTdH = dfTdH/newjac%dfHdH
        dfTdT = -dFTdH*newjac%dFHdT
      end block

      call accel%restart
      do iter = 1, max_iter
        udot = (u - u0)/dt
        call pd%compute_f(C0, Cdot, u(1:n), u(n+1), u(n+2), u(n+3), udot(1:n), udot(n+1), udot(n+2), &
            udot(n+3), du(1:n), du(n+1), du(n+2))
        du(n+3) = udot(n+2) - Hdot
        
        ! apply preconditioner to du
        block ! truchas-like
          call newjac%lower_solve(du(1:n), du(n+1), du(n+2))
          du(n+3) = du(n+3) - dfTdH * du(n+2)
          du(n+3) = du(n+3) / dfTdT
          call newjac%upper_solve(du(1:n), du(n+1), du(n+2), du(n+3))
        end block
        
        call accel%accel_update(du)
        u = u - du
        
        error = maxval(abs(du))
        !write(error_unit,'(a,i0,a,es8.2,a,*(es8.2,:,1x))') 'iter=', iter, ', max |du|=', error, ', |du|=', abs(du)
        if (error <= tol) then
          write(error_unit,'(a,i0,a,es8.2)') 'converged: ', iter, ' iterations, error=', error
          !write(output_unit,'(*(es12.5,:,1x))') t, u
          if (u(n+1) > 0) then
            C_liq = u(1:n)/u(n+1)
            C_liq_max = max(C_liq, C_liq_max)
          else
            C_liq = C_liq_max
          end if
          if (u(n+1) < 1) then
            C_sol = (C0 - u(1:n)) / (1-u(n+1))
          else
            C_sol = 0
          end if
          write(output_unit,'(*(es12.5,:,1x))') u(n+3), u(n+1), u(1:n), C_liq, C_sol, u(n+2)
          exit
        end if        
      end do
 
      if (iter > max_iter) then
        write(error_unit,*) 'error=', error
        stop 10
      end if
      call uhist%record_state(t, u)
      
    end do

  end subroutine timestep

end program
