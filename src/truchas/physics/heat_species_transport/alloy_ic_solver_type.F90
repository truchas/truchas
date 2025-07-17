!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module alloy_ic_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use alloy_vector_type
  use alloy_model_type
  use truchas_logging_services
  use parallel_communication
  use parameter_list_type
  implicit none
  private

  type, public :: alloy_ic_solver
    private
    type(unstr_mesh), pointer :: mesh  => null()  ! reference only -- do not own
    type(alloy_model), pointer :: model => null()  ! reference only -- do not own
    type(parameter_list), pointer :: params => null() ! reference only -- do not own
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_udot
  end type alloy_ic_solver

contains

  subroutine init (this, model, params)

    class(alloy_ic_solver), intent(out) :: this
    type(alloy_model), intent(in), target :: model
    type(parameter_list), pointer :: params

    this%model => model
    this%mesh  => model%mesh
    this%params => params

  end subroutine init


  subroutine compute(this, t, temp, u, udot)

    class(alloy_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: temp(:)
    type(alloy_vector), intent(inout) :: u, udot ! data is intent(out)
    target :: u

    real(r8), pointer :: state(:,:)

    call TLS_info('')
    call TLS_info('Computing consistent initial state for HT solver ...')

    call u%setval(0.0_r8)

    !! Set cell temperatures (U%TC)
    ASSERT(size(temp) >= this%mesh%ncell_onP)
    u%tc(:this%mesh%ncell_onP) = temp(:this%mesh%ncell_onP)
    call this%mesh%cell_imap%gather_offp(u%tc)

    !! Solve for consistent face temperatures.
    !! Averaging adjacent cell temps provides a cheap initial guess.
    call average_to_faces(this%mesh, u%tc, u%tf)
    if (allocated(this%model%bc_dir)) then
      call this%model%bc_dir%compute(t)
      u%tf(this%model%bc_dir%index) = this%model%bc_dir%value
    end if
    state(1:this%mesh%ncell,1:1) => u%tc
    call compute_face_temp(this%model, t, state, u, this%params)

    !! Compute the cell enthalpy density (U%HC)
    call this%model%H_of_T%compute_value(state, u%hc)

    call this%model%alloy%compute_g(u%tc, u%lf)

    call compute_udot(this, t, u, udot)

  end subroutine compute


  subroutine compute_udot (this, t, u, udot)

    class(alloy_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    type(alloy_vector), intent(inout) :: u    ! data is intent(in)
    type(alloy_vector), intent(inout) :: udot ! data is intent(out)
    target :: u, udot

    integer :: j, n
    real(r8), pointer :: state(:,:)
    real(r8) :: dt, Tmin, Tmax

    type(alloy_vector), target :: f

    call TLS_info('')
    call TLS_info('Computing consistent initial state derivative for HT solver ...')

    !TODO: The existing prop_mesh_func%compute_value function expects a rank-2
    !      state array. This is a workaround until prop_mesh_func is redesigned.
    state(1:this%mesh%ncell,1:1) => u%tc

    call udot%setval(0.0_r8)
    call f%init(u)
    call this%model%compute_f(t, u, udot, f)  ! only care about f%tc

    !! Compute Hdot
    udot%hc = -f%tc/this%mesh%volume

    !! Forward Euler time step (only cell enthalpy is advanced)
    call this%params%get('dt', dt)
    call udot%update(1.0_r8, u, dt) ! udot = u + dt*udot

    !! Compute the advanced cell temps from the advanced cell enthalpies
    do j = 1, size(u%tc)
      Tmin = u%tc(j) - 1
      Tmax = u%tc(j) + 1
      call this%model%T_of_H%compute(j, udot%hc(j), Tmin, Tmax, udot%tc(j))
    end do
    state(1:this%mesh%ncell,1:1) => udot%tc ! at advanced solution

    !! Set consistent advanced face temps.  Use their initial
    !! conditions as the initial guess for the solution procedure.
    call compute_face_temp(this%model, t+dt, state, udot, this%params)

    call this%model%alloy%compute_g(udot%tc, udot%lf)

    !! Forward Euler approximation to the time derivative at T.
    !f = (f - u) / dt
    call udot%update(-1.0_r8, u)
    call udot%scale(1.0_r8/dt)

  end subroutine compute_udot

  !! Given a cell-based field, this auxiliary subroutine produces a face-based
  !! field whose value on each face is the average of the values on adjacent
  !! cells.

  subroutine average_to_faces(mesh, ucell, uface)

    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in)  :: ucell(:)
    real(r8), intent(out) :: uface(:)

    integer :: j
    integer :: scale(size(uface))

    ASSERT(size(ucell) == mesh%ncell)
    ASSERT(size(uface) == mesh%nface)

    uface = 0.0_r8
    scale = 0.0_r8
    do j = 1, mesh%ncell
      associate (cface => mesh%cface(mesh%xcface(j):mesh%xcface(j+1)-1))
        uface(cface) = uface(cface) + ucell(j)
        scale(cface) = scale(cface) + 1
      end associate
    end do
    call mesh%face_imap%gather_offp(uface)
    call mesh%face_imap%gather_offp(scale)

    where (scale == 0)
      uface = 0.0_r8
    elsewhere
      uface = uface / scale
    end where

  end subroutine average_to_faces

  !! This auxillary procedure solves for the consistent face temperatures given
  !! cell temperatures.

  subroutine compute_face_temp(this, t, state, u, params)

    use mfd_diff_matrix_type
    use pcsr_precon_class
    use pcsr_precon_factory
    use parameter_list_type
    use nka_type
    use parallel_communication, only: global_sum, global_maxval, global_all

    type(alloy_model), intent(inout) :: this
    real(r8), intent(in) :: t, state(:,:) ! lower bound is now 1
    type(alloy_vector), intent(inout) :: u
    type(parameter_list) :: params

    integer :: n, j, stat, num_itr, max_iter, iter
    type(mfd_diff_matrix), target :: matrix
    real(r8) :: atol, rtol, error, r0_err, r_err, dT_max
    real(r8), allocatable :: z(:)
    type(alloy_vector) :: udot, f
    type(nka) :: nka_obj
    procedure(pardp), pointer :: dp
    character(80) :: string
    character(:), allocatable :: errmsg
    class(pcsr_precon), allocatable :: precon

    call TLS_info ('  Computing consistent face temperatures ...')

    !! Initial residual of the face temperature equations.
    call udot%init(u)
    call udot%setval(0.0_r8)
    call f%init(u)
    call this%compute_f(t, u, udot, f)
    associate (r => f%tf(1:this%mesh%nface_onP))
      r0_err = sqrt(global_sum(norm2(r)**2))
    end associate
    if (r0_err == 0.0_r8) return
    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      write(string,'(4x,a,es9.2)') '||Rface(0)||_2 = ', r0_err
      call TLS_info (string)
    end if

    !! Form the Jacobian (approx) of the face temperature system.
    call matrix%init(this%disc)
    call make_matrix(matrix)

    call alloc_pcsr_precon(precon, matrix%a22, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_FACE_TEMP: ' // errmsg)
    call precon%compute

    call params%get('atol-temp', atol, default=0.0_r8)
    call params%get('rtol-temp', rtol, default=0.0_r8)
    call params%get('max-iter', max_iter)

    !! Simple Picard iteration.
    allocate(z(this%mesh%nface_onP))
    call nka_obj%init(size(z), 10)
    dp => pardp
    call nka_obj%set_dot_prod(dp)

    associate (Tface => u%tf(1:this%mesh%nface_onP), &
               r => f%tf(1:this%mesh%nface_onP))
      iter = 0
      do ! until converged

        iter = iter + 1
        if (iter > max_iter) exit

        z = r
        call precon%apply(z)
        call nka_obj%accel_update(z) ! nonlinear krylov acceleration (NKA)

        Tface = Tface - z

        dT_max = global_maxval(abs(z))

        !! Recompute the face temp residual.
        call this%compute_f(t, u, udot, f)
        r_err = sqrt(global_dot_product(r,r))
        if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
          write(string,'(4x,a,i0,2(a,es9.2))') &
              '||Rface(', iter, ')||=', r_err, ', ||ΔTface||_max=', dT_max
          call TLS_info (string)
        end if

        !! Check for convergence.
        if (global_all(abs(z) <= atol + rtol * abs(Tface))) exit
      end do
    end associate

    write(string,'(4x,a,i0,3(a,es9.2))') &
        '||Rface(', iter, ')||=', r_err, ', ||Rface(0)||=', r0_err, &
        ', ||ΔTface||_max=', dT_max
    call TLS_info (string)

    if (iter > max_iter) call TLS_info ('    WARNING: face temperatures not converged')

  contains

    !! This auxillary routine computes an approximate Jacobian matrix for the
    !! time-independent heat equation (diffusion operator). This is a cell and
    !! face system matrix, but we are only interested in the (2,2) block (face-
    !! face).
    !!
    !! Note that this is nearly identical to code in HTSD_precon_type:compute
    !! that computes the preconditioner matrix.  CAN WE AVOID DUPLICATION?

    subroutine make_matrix (matrix)

      type(mfd_diff_matrix), intent(inout) :: matrix

      integer :: n, index
      real(r8) :: D(this%mesh%ncell)
      real(r8), allocatable :: values(:)

      !! Compute the diffusion coefficient.
      call this%conductivity%compute_value(state, D)

      !! Compute the basic diffusion matrix
      call matrix%compute (D)

      !! Dirichlet boundary condition fixups.
      if (allocated(this%bc_dir)) then
        call this%bc_dir%compute(t)
        call matrix%set_dir_faces(this%bc_dir%index)
      end if

      !! External HTC boundary condition contribution.
      if (allocated(this%bc_htc)) then
        call this%bc_htc%compute(t, u%tf)
        associate (index => this%bc_htc%index, &
                   deriv => this%bc_htc%deriv)
          call matrix%incr_face_diag(index, deriv)
        end associate
      end if

      !! Simple radiation boundary condition contribution.
      if (allocated(this%bc_rad)) then
        call this%bc_rad%compute(t, u%tf)
        associate (index => this%bc_rad%index, &
                   deriv => this%bc_rad%deriv)
          call matrix%incr_face_diag(index, deriv)
        end associate
      end if

      !! Experimental evaporation heat flux contribution.
      if (allocated(this%evap_flux)) then
        call this%evap_flux%compute_deriv(t, u%tf)
        associate (index => this%evap_flux%index, &
                   deriv => this%evap_flux%deriv)
          call matrix%incr_face_diag(index, this%mesh%area(index)*deriv)
        end associate
      endif

      !! Internal HTC interface condition contribution.
      if (allocated(this%ic_htc)) then
        call this%ic_htc%compute(t, u%tf)
        associate (index => this%ic_htc%index, &
                   deriv => this%ic_htc%deriv)
          call matrix%incr_interface_flux3(index, deriv) !TODO: rename these methods
        end associate
      end if

      !! Internal gap radiation condition contribution.
      if (allocated(this%ic_rad)) then
        call this%ic_rad%compute(t, u%tf)
        associate (index => this%ic_rad%index, &
                   deriv => this%ic_rad%deriv)
          call matrix%incr_interface_flux3(index, deriv) !TODO: rename these methods
        end associate
      end if

    end subroutine make_matrix

  end subroutine compute_face_temp

  function pardp (a, b) result (dp)
    use parallel_communication, only: global_dot_product
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    dp = global_dot_product(a, b)
  end function pardp

end module alloy_ic_solver_type
