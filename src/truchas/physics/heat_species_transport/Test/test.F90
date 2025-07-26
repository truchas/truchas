program test_multicomp_lever_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use parameter_list_json
  use alloy_back_diff_type
  
  integer :: inlun, stat
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

    real(r8) :: T0, T1, tlast, t, dt, tol, error
    real(r8), allocatable :: u0(:), u(:), udot(:), du(:), jac(:,:), C_sol(:), C_liq(:), C_liq_max(:)
    type(state_history) :: uhist
    type(nka) :: accel
    integer :: n, iter, max_iter

    n = pd%num_comp + 3
    allocate(u0(n), u(n), udot(n), du(n), jac(n,n))
    n = pd%num_comp
    allocate(C_sol(n), C_liq(n), C_liq_max(n))

    call params%get('T0', T0)
    call params%get('T1', T1)

    !! Starting state (pure liquid
    u(1:n) = pd%C0
    u(n+1) = 1
    u(n+2) = pd%H_liq%eval([T0])
    u(n+3) = T0

    udot(1:n) = 0
    udot(n+1) = 0
    udot(n+2) = pd%Hdot
    udot(n+3) = pd%Hdot ! assumes unit heat capacity (per unit volume)
    
    t = 0
    call uhist%init(2, t, u, udot) ! maintain 2 solution vectors
    call accel%init(size(u), 4)

    call params%get('dt', dt)
    call params%get('tol', tol)
    call params%get('max-iter', max_iter)

    C_liq_max = 0
    
    do while (u(4) > T1)
    
      tlast = t
      t = t+dt
      
      u0 = u
      call uhist%interp_state(t, u, order=1) ! linear extrapolation
      
      call pd%compute_f1_jac(u, udot, dt, jac)
      call fct(jac)

      call accel%restart
      do iter = 1, max_iter
        udot = (u - u0)/dt
        call pd%compute_f1(u, udot, du)
        
        ! apply preconditioner to du
        call slv(jac, du)
        
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
            C_sol = (pd%C0 - u(1:n)) / (1-u(n+1))
          else
            C_sol = pd%C0
          end if
          write(output_unit,'(*(es12.5,:,1x))') u(n+3), u(n+1), u(1:n), C_liq, C_sol, u(n+2)
          exit
        end if        
      end do
      
      if (iter > max_iter) then
        stop 1
      end if
      call uhist%record_state(t, u)
      
    end do

  end subroutine timestep

end program
