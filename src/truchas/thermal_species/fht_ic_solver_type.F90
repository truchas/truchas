!!
!! FHT_IC_SOLVER_TYPE
!!
!! This module defines the initial-condition solver used by the flow-coupled
!! thermal transport solver. It takes user-facing cell temperature initial
!! conditions and computes the full FHT vector state and initial time
!! derivative needed by the nonadaptive integrator.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! Face temperatures and radiosities are solved to satisfy their algebraic
!! residual equations. Cell temperature derivatives, face temperature
!! derivatives, and radiosity derivatives are computed by finite differencing
!! advanced states with similarly solved algebraic variables.
!!

#include "f90_assert.fpp"

module fht_ic_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use fht_model_type
  use fht_vector_type
  use truchas_logging_services
  use parallel_communication
  use parameter_list_type
  implicit none
  private

  type, public :: fht_ic_solver
    private
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(fht_model), pointer :: model => null() ! unowned reference
    type(parameter_list), pointer :: params => null() ! unowned reference
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_udot
  end type

contains

  subroutine init(this, model, params)
    class(fht_ic_solver), intent(out) :: this
    type(fht_model), intent(in), target :: model
    type(parameter_list), intent(inout), target :: params
    this%model => model
    this%mesh => model%mesh
    this%params => params
  end subroutine

  subroutine compute(this, t, temp, u, udot)

    class(fht_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: temp(:)
    type(fht_vector), intent(inout) :: u, udot

    call TLS_info('')
    call TLS_info('Computing consistent initial state for FHT solver ...')

    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      ASSERT(size(temp) == ncell_onP)
      ASSERT(associated(this%model%void_cell))
      ASSERT(associated(this%model%void_face))

      call u%setval(0.0_r8)
      u%tc(:ncell_onP) = temp
      where (this%model%void_cell(:ncell_onP)) u%tc(:ncell_onP) = this%model%void_temp
      call this%mesh%cell_imap%gather_offp(u%tc)

      !! Solve for consistent face temperatures and radiosities.
      !! Averaging adjacent cell temps provides a cheap initial guess.
      call average_to_faces(this%mesh, u%tc, u%tf, this%model%void_cell)
      if (allocated(this%model%thermal%bc_dir)) then
        call this%model%thermal%bc_dir%compute(t)
        u%tf(this%model%thermal%bc_dir%index) = this%model%thermal%bc_dir%value
      end if
      where (this%model%void_face(:nface_onP)) u%tf(:nface_onP) = this%model%void_temp
      call compute_face_temp(this%model, t, u%tc, u, this%params)
    end associate

    call compute_udot(this, t, u, udot)

  end subroutine compute


  subroutine compute_udot(this, t, u, udot)

    class(fht_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    type(fht_vector), intent(inout) :: u    ! data is intent(in)
    type(fht_vector), intent(inout) :: udot ! data is intent(out)

    integer :: j
    real(r8) :: dt, Hadv, Tmin, Tmax
    real(r8), allocatable :: Hcell(:), Hdot(:)
    type(fht_vector) :: f, uadv

    call TLS_info('')
    call TLS_info('Computing consistent initial state derivative for FHT solver ...')

    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      ASSERT(associated(this%model%void_cell))
      ASSERT(associated(this%model%void_face))

      call u%gather_offp

      call udot%setval(0.0_r8)
      call f%init(u)
      allocate(Hcell(ncell_onP), Hdot(ncell_onP))
      Hdot = 0.0_r8
      call this%model%residual(t, u, Hdot, f)

      !! Compute Hdot.
      Hdot = -f%tc(:ncell_onP) / this%mesh%volume(:ncell_onP)
      where (this%model%void_cell(:ncell_onP)) Hdot = 0.0_r8

      call this%params%get('dt', dt)

      !! Forward Euler time step of the cell enthalpy, followed by inversion
      !! back to cell temperature. FHT vectors do not carry cell enthalpy, so
      !! the enthalpy values are held locally here.
      call this%model%thermal%H_of_T%compute_value(u%tc(:ncell_onP), Hcell)
      call uadv%init(u)
      call uadv%copy(u)
      do j = 1, ncell_onP
        if (this%model%void_cell(j)) then
          uadv%tc(j) = this%model%void_temp
          cycle
        end if
        Hadv = Hcell(j) + dt*Hdot(j)
        Tmin = u%tc(j) - 1.0_r8
        Tmax = u%tc(j) + 1.0_r8
        call this%model%T_of_H%compute(j, Hadv, Tmin, Tmax, uadv%tc(j))
      end do
      call this%mesh%cell_imap%gather_offp(uadv%tc)

      !! Set consistent advanced face temps and radiosities. Use their initial
      !! conditions as the initial guess for the solution procedure.
      call compute_face_temp(this%model, t+dt, uadv%tc, uadv, this%params)

      call uadv%update(-1.0_r8, u)
      call uadv%scale(1.0_r8/dt)

      udot%tc(:ncell_onP) = uadv%tc(:ncell_onP)
      udot%tf(:nface_onP) = uadv%tf(:nface_onP)
      if (this%model%vf_rad%is_active()) then
        ASSERT(allocated(udot%qrad))
        udot%qrad(:) = uadv%qrad(:)
      end if
      call udot%gather_offp
    end associate

  end subroutine compute_udot


  !! Given a cell-based field, this auxiliary subroutine produces a face-based
  !! field whose value on each face is the average of the values on adjacent
  !! cells. SKIP is an optional cell-based mask array identifying cells whose
  !! values are to be ignored. Faces with no adjacent cells have the value 0.

  subroutine average_to_faces(mesh, ucell, uface, skip)

    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in) :: ucell(:)
    real(r8), intent(out) :: uface(:)
    logical, intent(in), optional :: skip(:)

    integer :: j
    integer :: scale(size(uface))

    ASSERT(size(ucell) == mesh%ncell)
    ASSERT(size(uface) == mesh%nface)

    uface = 0.0_r8
    scale = 0
    do j = 1, mesh%ncell
      if (present(skip)) then
        if (skip(j)) cycle
      end if
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


  !! This auxiliary procedure solves for the consistent face temperatures given
  !! cell temperatures. It also computes the consistent radiosities if enclosure
  !! radiation problems are present.

  subroutine compute_face_temp(this, t, temp, u, params)

    use mfd_diff_matrix_type
    use pcsr_precon_class
    use pcsr_precon_factory
    use parameter_list_type
    use nka_type
    use parallel_communication, only: global_dot_product, global_maxval, global_all

    type(fht_model), intent(inout), target :: this
    real(r8), intent(in) :: t, temp(:)
    type(fht_vector), intent(inout) :: u
    type(parameter_list) :: params

    integer :: stat, num_itr, max_iter, iter
    type(mfd_diff_matrix), target :: matrix
    real(r8) :: atol, rtol, error, r0_err, r_err, dT_max
    real(r8), allocatable :: z(:), zero_hdot(:)
    type(fht_vector) :: f
    type(nka) :: nka_obj
    procedure(pardp), pointer :: dp
    character(160) :: string
    character(:), allocatable :: errmsg
    class(pcsr_precon), allocatable :: precon

    call TLS_info('  Computing consistent face temperatures and radiosities ...')
    ASSERT(associated(this%void_cell))
    ASSERT(associated(this%void_face))

    !! Solve the radiosity components given the approximate face temperatures.
    if (this%vf_rad%is_active()) then
      u%qrad(:) = 0.0_r8
      call this%vf_rad%solve_radiosity(t, u%tf, u%qrad, stat, num_itr, error)
      if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
        write(string,'(4x,a,es9.2," (",i0,")")') 'radiosity: |r|/|b|=', error, num_itr
        call TLS_info(string)
      end if
      if (stat /= 0) call TLS_info('      WARNING: radiosities not converged')
    end if

    !! Initial residual of the face temperature equations.
    call f%init(u)
    allocate(zero_hdot(this%mesh%ncell_onP))
    zero_hdot = 0.0_r8
    call this%residual(t, u, zero_hdot, f)
    associate (r => f%tf(1:this%mesh%nface_onP))
      r0_err = sqrt(global_dot_product(r, r))
    end associate
    if (r0_err == 0.0_r8) return
    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      write(string,'(4x,a,es9.2)') '||Rface(0)||_2 = ', r0_err
      call TLS_info(string)
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
    r_err = r0_err
    dT_max = 0.0_r8

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

        !! Solve the radiosity components given the new face temperatures.
        if (this%vf_rad%is_active()) then
          call this%vf_rad%solve_radiosity(t, u%tf, u%qrad, stat, num_itr, error)
          if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
            write(string,'(4x,a,es9.2," (",i0,")")') 'radiosity: |r|/|b|=', error, num_itr
            call TLS_info(string)
          end if
          if (stat /= 0) call TLS_info('    WARNING: radiosities not converged')
        end if

        !! Recompute the face temp residual.
        call this%residual(t, u, zero_hdot, f)
        r_err = sqrt(global_dot_product(r, r))
        if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
          write(string,'(4x,a,i0,2(a,es9.2))') &
              '||Rface(', iter, ')||=', r_err, ', ||Delta Tface||_max=', dT_max
          call TLS_info(string)
        end if

        if (global_all(abs(z) <= atol + rtol * abs(Tface))) exit
      end do
    end associate

    write(string,'(4x,a,i0,3(a,es9.2))') &
        '||Rface(', iter, ')||=', r_err, ', ||Rface(0)||=', r0_err, &
        ', ||Delta Tface||_max=', dT_max
    call TLS_info(string)

    if (iter > max_iter) call TLS_info('    WARNING: face temperatures not converged')

  contains

    !! This auxiliary routine computes an approximate Jacobian matrix for the
    !! time-independent heat equation (diffusion operator). This is a cell and
    !! face system matrix, but we are only interested in the (2,2) block (face-
    !! face). That part is approximate in that it ignores the dependence of
    !! the radiosity on face temperatures (assuming the radiosities have been
    !! solved for). Moreover the Jacobian would be constant except for the T^4
    !! terms present in radiation BC.
    !!
    !! This mirrors the preconditioner matrix construction; keep the two paths
    !! aligned when changing thermal flux or radiation Jacobian contributions.

    subroutine make_matrix(matrix)

      type(mfd_diff_matrix), intent(inout) :: matrix

      integer :: j, n
      real(r8) :: D(this%mesh%ncell)
      integer, allocatable :: void_dir_faces(:)

      call this%thermal%conductivity%compute_value(temp, D)
      where (this%void_cell) D = 0.0_r8

      call matrix%compute(D)

      if (allocated(this%thermal%bc_dir)) then
        call this%thermal%bc_dir%compute(t)
        call matrix%set_dir_faces(this%thermal%bc_dir%index)
      end if

      call this%thermal%add_flux_bc_jacobian(t, u%tf, this%mesh%area, matrix, &
          this%void_face)

      n = count(this%void_face)
      if (n > 0) then
        allocate(void_dir_faces(n))
        n = 0
        do j = 1, this%mesh%nface
          if (this%void_face(j)) then
            n = n + 1
            void_dir_faces(n) = j
          end if
        end do
        call matrix%set_dir_faces(void_dir_faces)
      end if

      if (this%vf_rad%is_active()) then
        call this%vf_rad%add_heat_precon_deriv(t, u%tf, this%mesh%area, matrix)
      end if

    end subroutine make_matrix

  end subroutine compute_face_temp


  function pardp(a, b) result(dp)
    use parallel_communication, only: global_dot_product
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    dp = global_dot_product(a, b)
  end function

end module fht_ic_solver_type
