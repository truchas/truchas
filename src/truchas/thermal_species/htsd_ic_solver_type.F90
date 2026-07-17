!!
!! HTSD_IC_SOLVER_TYPE
!!
!! This module defines the initial-condition solver used by the coupled
!! thermal/species transport solver. It takes user-facing cell temperature and
!! concentration initial conditions and computes the full coupled vector state
!! and initial time derivative needed by the implicit integrator.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The thermal face temperatures and radiosities are solved to satisfy their
!! algebraic residual equations. Species face concentrations are solved to
!! satisfy their algebraic residual equations, and their time derivatives are
!! computed by finite differencing advanced states with similarly solved face
!! values. Cell temperature derivatives are computed by finite differencing
!! advanced cell enthalpies and concentrations and inverting the enthalpy
!! relation cell by cell.
!!

#include "f90_assert.fpp"

module htsd_ic_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use htsd_model_type
  use htsd_vector_type
  use truchas_logging_services
  use parallel_communication
  use parameter_list_type
  use string_utilities, only: i_to_c
  implicit none
  private

  type, public :: htsd_ic_solver
    private
    type(unstr_mesh), pointer :: mesh  => null() ! unowned reference
    type(htsd_model), pointer :: model => null() ! unowned reference
    type(parameter_list), pointer :: params => null() ! unowned reference
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_udot
  end type htsd_ic_solver

contains

  subroutine init (this, model, params)

    class(htsd_ic_solver), intent(out) :: this
    type(htsd_model), intent(in), target :: model
    type(parameter_list), intent(inout), target :: params

    this%model => model
    this%mesh  => model%mesh
    this%params => params

  end subroutine init


  subroutine compute (this, t, temp, conc, u, udot)

    class(htsd_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), optional :: temp(:), conc(:,:)
    type(htsd_vector), intent(inout) :: u, udot

    integer :: n
    type(htsd_vector) :: f
    real(r8), allocatable :: Tface(:), Cface(:)
    real(r8) :: dt

    call TLS_info ('')
    call TLS_info ('Computing consistent initial state for HT/SD solver ...')

    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      call u%setval(0.0_r8)

      !! Set cell temperatures.
      ASSERT(present(temp))
      ASSERT(size(temp) == ncell_onP)
      u%tc(:ncell_onP) = temp
      if (allocated(this%model%void_cell)) &
          where (this%model%void_cell(:ncell_onP)) u%tc(:ncell_onP) = this%model%void_temp

      !! Set cell concentrations.
      ASSERT(present(conc))
      ASSERT(size(conc,dim=1) == ncell_onP)
      ASSERT(size(conc,dim=2) == this%model%num_comp)
      do n = 1, this%model%num_comp
        u%cc(:ncell_onP,n) = conc(:,n)
        if (allocated(this%model%void_cell)) &
            where (this%model%void_cell(:ncell_onP)) u%cc(:ncell_onP,n) = 0.0_r8
      end do

      !! Gather cell temperatures and concentrations used by model parameters.
      call u%gather_offp

      !! Set consistent face temperatures and radiosities.
      !! Approximate face temps by averaging adjacent cell temps; provides
      !! a cheap initial guess for the solution procedure that follows.
      allocate(Tface(this%mesh%nface))
      call average_to_faces (this%mesh, u%tc, Tface, this%model%void_cell)
      !! Overwrite face temperatures with Dirichlet boundary data.
      if (allocated(this%model%thermal%bc_dir)) then
        call this%model%thermal%bc_dir%compute(t)
        Tface(this%model%thermal%bc_dir%index) = this%model%thermal%bc_dir%value
      end if
      !! Overwrite void face components with dummy value.
      if (allocated(this%model%void_face)) where (this%model%void_face) Tface = this%model%void_temp
      u%tf(:nface_onP) = Tface(:nface_onP)
      call u%gather_offp
      deallocate(Tface)
      !! Now solve the algebraic equations for face temps and radiosities.
      call compute_face_temp (this%model, t, u%tc, u%cc, u, this%params)

      !! Initial guess for the algebraic face concentration unknowns.
      allocate(Cface(this%mesh%nface))
      do n = 1, this%model%num_comp
        call average_to_faces (this%mesh, u%cc(:,n), Cface, this%model%void_cell)
        if (allocated(this%model%species(n)%bc_dir)) then
          call this%model%species(n)%bc_dir%compute(t)
          Cface(this%model%species(n)%bc_dir%index) = this%model%species(n)%bc_dir%value
        end if
        if (allocated(this%model%void_face)) where (this%model%void_face) Cface = 0.0_r8
        u%cf(:nface_onP,n) = Cface(:nface_onP)
      end do
      deallocate(Cface)
      call u%gather_offp
      do n = 1, this%model%num_comp
        call compute_face_conc(this%model, n, t, u, this%params)
      end do

      !! Set the cell heat density.
      u%hc(:ncell_onP) = 0.0_r8
      call udot%setval(0.0_r8)
      call f%init(u)
      call this%model%residual(t, u, udot, f)
      !! Extract the heat densities from F.
      u%hc(:ncell_onP) = -f%hc(:ncell_onP)
      if (allocated(this%model%void_cell)) &
          where (this%model%void_cell(:ncell_onP)) u%hc(:ncell_onP) = 0.0_r8
      call u%gather_offp

      call TLS_info ('')
      call TLS_info ('Computing consistent initial state derivative for HT/SD solver ...')

      !! Set the cell heat density time derivative.
      udot%hc(:ncell_onP) = -f%tc(:ncell_onP) / this%mesh%volume(:ncell_onP)
      if (allocated(this%model%void_cell)) &
          where (this%model%void_cell(:ncell_onP)) udot%hc(:ncell_onP) = 0.0_r8

      !! Set the cell concentration time derivative, extracted from F.
      do n = 1, this%model%num_comp
        udot%cc(:ncell_onP,n) = -f%cc(:ncell_onP,n) / this%mesh%volume(:ncell_onP)
        if (allocated(this%model%void_cell)) &
            where (this%model%void_cell(:ncell_onP)) udot%cc(:ncell_onP,n) = 0.0_r8
      end do

      call this%params%get ('dt', dt)
      call compute_cell_temp_deriv(this, u, udot, dt)
      call udot%gather_offp

      !! Compute discrete initial time derivatives of the remaining algebraic
      !! variables (face temps, concs, and radiosities) by finite difference.
      !! The state is advanced by a small time step using the cell derivatives,
      !! and then the associated advanced face values are solved for.

      call f%copy(u)
      call f%update(dt, udot)
      call f%gather_offp

      !! Set consistent advanced face temps and radiosities.  Use their initial
      !! conditions as the initial guess for the solution procedure.
      call compute_face_temp (this%model, t+dt, f%tc, f%cc, f, this%params)

      !! Set consistent advanced face concentrations. Use their initial
      !! conditions as the initial guess for the solution procedure.
      do n = 1, this%model%num_comp
        call compute_face_conc(this%model, n, t+dt, f, this%params)
      end do

      call f%update(-1.0_r8, u)
      call f%scale(1.0_r8/dt)

      !! Finite difference approx of face temp and radiosity derivatives.
      udot%tf(:nface_onP) = f%tf(:nface_onP)
      if (this%model%vf_rad%is_active()) then
        ASSERT(allocated(udot%qrad))
        udot%qrad(:) = f%qrad(:)
      end if

      !! Finite difference approximation of face concentration derivatives.
      udot%cf(:nface_onP,:) = f%cf(:nface_onP,:)
      call udot%gather_offp
    end associate

  end subroutine compute


  subroutine compute_udot (this, t, u, udot)

    class(htsd_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    type(htsd_vector), intent(inout) :: u
    type(htsd_vector), intent(inout) :: udot

    integer :: n
    type(htsd_vector) :: f
    real(r8) :: dt

    call TLS_info ('')
    call TLS_info ('Computing consistent initial state derivative for HT/SD solver ...')

    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      !! Gather cell temperatures and concentrations used by model parameters.
      call u%gather_offp

      call udot%setval(0.0_r8)
      call f%init(u)
      call this%model%residual(t, u, udot, f)

      !! Set the cell heat density time derivative.
      udot%hc(:ncell_onP) = -f%tc(:ncell_onP) / this%mesh%volume(:ncell_onP)
      if (allocated(this%model%void_cell)) &
          where (this%model%void_cell(:ncell_onP)) udot%hc(:ncell_onP) = 0.0_r8

      !! Set the cell concentration time derivative, extracted from F.
      do n = 1, this%model%num_comp
        udot%cc(:ncell_onP,n) = -f%cc(:ncell_onP,n) / this%mesh%volume(:ncell_onP)
        if (allocated(this%model%void_cell)) &
            where (this%model%void_cell(:ncell_onP)) udot%cc(:ncell_onP,n) = 0.0_r8
      end do

      call this%params%get ('dt', dt)
      call compute_cell_temp_deriv(this, u, udot, dt)
      call udot%gather_offp

      !! Compute discrete initial time derivatives of the remaining algebraic
      !! variables (face temps, concs, and radiosities) by finite difference.
      !! The state is advanced by a small time step using the cell derivatives,
      !! and then the associated advanced face values are solved for.

      call f%copy(u)
      call f%update(dt, udot)
      call f%gather_offp

      !! Set consistent advanced face temps and radiosities.  Use their initial
      !! conditions as the initial guess for the solution procedure.
      call compute_face_temp (this%model, t+dt, f%tc, f%cc, f, this%params)

      !! Set consistent advanced face concentrations. Use their initial
      !! conditions as the initial guess for the solution procedure.
      do n = 1, this%model%num_comp
        call compute_face_conc(this%model, n, t+dt, f, this%params)
      end do

      call f%update(-1.0_r8, u)
      call f%scale(1.0_r8/dt)

      !! Finite difference approx of face temp and radiosity derivatives.
      udot%tf(:nface_onP) = f%tf(:nface_onP)
      if (this%model%vf_rad%is_active()) then
        ASSERT(allocated(udot%qrad))
        udot%qrad(:) = f%qrad(:)
      end if

      !! Finite difference approximation of face concentration derivatives.
      udot%cf(:nface_onP,:) = f%cf(:nface_onP,:)
      call udot%gather_offp
    end associate

  end subroutine compute_udot

  subroutine compute_cell_temp_deriv (this, u, udot, dt)

    class(htsd_ic_solver), intent(inout) :: this
    type(htsd_vector), intent(in) :: u
    type(htsd_vector), intent(inout) :: udot
    real(r8), intent(in) :: dt

    integer :: j
    real(r8) :: Hadv, Tadv, Tmin, Tmax
    real(r8) :: Cadv(this%model%num_comp)

    do j = 1, this%mesh%ncell_onP
      if (allocated(this%model%void_cell)) then
        if (this%model%void_cell(j)) then
          udot%tc(j) = 0.0_r8
          cycle
        end if
      end if
      Hadv = u%hc(j) + dt*udot%hc(j)
      Cadv = u%cc(j,:) + dt*udot%cc(j,:)
      Tmin = u%tc(j) - 1.0_r8
      Tmax = u%tc(j) + 1.0_r8
      call this%model%T_of_H%compute(j, Cadv, Hadv, Tmin, Tmax, Tadv)
      udot%tc(j) = (Tadv - u%tc(j)) / dt
    end do

  end subroutine compute_cell_temp_deriv

  !! Given a cell-based field, this auxiliary subroutine produces a face-based
  !! field whose value on each face is the average of the values on adjacent
  !! cells.  SKIP is an optional cell-based mask array identifying cells whose
  !! values are to be ignored.  Faces with no adjacent cells have the value 0.

  subroutine average_to_faces (mesh, ucell, uface, skip)

    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in)  :: ucell(:)
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
  !! cell temperatures (and concentrations -- heat equation parameters may
  !! depend on concentration). It also computes the consistent radiosities if
  !! enclosure radiation problems are present.

  subroutine compute_face_temp (this, t, temp, conc, u, params)

    use mfd_diff_matrix_type
    use pcsr_precon_class
    use pcsr_precon_factory
    use parameter_list_type
    use nka_type
    use parallel_communication, only: global_maxval, global_all

    type(htsd_model), intent(inout), target :: this
    real(r8), intent(in) :: t, temp(:), conc(:,:)
    type(htsd_vector), intent(inout) :: u
    type(parameter_list) :: params

    integer :: stat, num_itr, max_iter, iter
    type(mfd_diff_matrix), target :: matrix
    real(r8) :: atol, rtol, error, r0_err, r_err, dT_max
    real(r8), allocatable :: z(:)
    type(htsd_vector) :: udot, f
    type(nka) :: nka_obj
    procedure(pardp), pointer :: dp
    character(80) :: string
    character(:), allocatable :: errmsg
    class(pcsr_precon), allocatable :: precon

    call TLS_info ('  Computing consistent face temperatures and radiosities ...')

    !! Solve for the radiosity components given the (approx) face temperatures.
    if (this%vf_rad%is_active()) then
      u%qrad(:) = 0.0_r8
      call this%vf_rad%solve_radiosity(t, u%tf, u%qrad, stat, num_itr, error)
      if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
        write(string,'(2x,a,es9.2," (",i0,")")') 'radiosity: |r|/|b|=', error, num_itr
        call TLS_info (string)
      end if
      if (stat /= 0) call TLS_info ('  WARNING: radiosities not converged')
    end if

    !! Initial residual of the face temperature equations.
    call udot%init(u)
    call udot%setval(0.0_r8)
    call f%init(u)
    call this%residual(t, u, udot, f)
    associate (r => f%tf(1:this%mesh%nface_onP))
      r0_err = sqrt(global_dot_product(r, r))
    end associate
    if (r0_err == 0.0_r8) return
    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      write(string,'(2x,a,es9.2)') '||Rface(0)||_2 = ', r0_err
      call TLS_info (string)
    end if

    !! Form the Jacobian (approx) of the face temperature system.
    call matrix%init (this%disc)
    call make_matrix (matrix)

    call alloc_pcsr_precon (precon, matrix%a22, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_FACE_TEMP: ' // errmsg)
    call precon%compute

    call params%get ('atol-temp', atol, default=0.0_r8)
    call params%get ('rtol-temp', rtol, default=0.0_r8)
    call params%get ('max-iter', max_iter)

    !! Simple Picard iteration.
    allocate(z(this%mesh%nface_onP))
    call nka_obj%init (size(z), 10)
    dp => pardp
    call nka_obj%set_dot_prod (dp)

    associate (Tface => u%tf(1:this%mesh%nface_onP), &
               r => f%tf(1:this%mesh%nface_onP))
      iter = 0
      do ! until converged

        iter = iter + 1
        if (iter > max_iter) exit

        z = r
        call precon%apply (z)
        call nka_obj%accel_update (z) ! nonlinear krylov acceleration (NKA)

        Tface = Tface - z

        dT_max = global_maxval(abs(z))

        !! Solve the radiosity components given the new face temperatures.
        if (this%vf_rad%is_active()) then
          call this%vf_rad%solve_radiosity(t, u%tf, u%qrad, stat, num_itr, error)
          if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
            write(string,'(2x,a,es9.2," (",i0,")")') 'radiosity: |r|/|b|=', error, num_itr
            call TLS_info (string)
          end if
          if (stat /= 0) call TLS_info ('  WARNING: radiosities not converged')
        end if

        !! Recompute the face temp residual.
        call this%residual(t, u, udot, f)
        r_err = sqrt(global_dot_product(r, r))
        if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
          write(string,'(2x,a,i0,2(a,es9.2))') &
              '||Rface(', iter, ')||=', r_err, ', ||ΔTface||_max=', dT_max
          call TLS_info (string)
        end if

        !! Check for convergence.
        if (global_all(abs(z) <= atol + rtol * abs(Tface))) exit
      end do
    end associate

    write(string,'(2x,a,i0,3(a,es9.2))') &
        '||Rface(', iter, ')||=', r_err, ', ||Rface(0)||=', r0_err, &
        ', ||ΔTface||_max=', dT_max
    call TLS_info (string)

    if (iter > max_iter) call TLS_info ('  WARNING: face temperatures not converged')

  contains

    !! This auxiliary routine computes an approximate Jacobian matrix for the
    !! time-independent heat equation (diffusion operator). This is a cell and
    !! face system matrix, but we are only interested in the (2,2) block (face-
    !! face).  That part is approximate in that it ignores the dependence of
    !! the radiosity on face temperatures (assuming the radiosities have been
    !! solved for). Moreover the Jacobian would be constant except for the T^4
    !! terms present in radiation BC.
    !!
    !! This mirrors the preconditioner matrix construction; keep the two paths
    !! aligned when changing thermal flux or radiation Jacobian contributions.

    subroutine make_matrix (matrix)

      type(mfd_diff_matrix), intent(inout) :: matrix

      real(r8) :: D(this%mesh%ncell), Tface_wide(this%mesh%nface)

      !! Compute the diffusion coefficient.
      call this%thermal%conductivity%compute_value(temp, conc, D)
      if (allocated(this%void_cell)) then
        where (this%void_cell) D = 0.0_r8
      end if

      !! Compute the basic diffusion matrix
      call matrix%compute (D)

      !! Dirichlet boundary condition fixups.
      if (allocated(this%thermal%bc_dir)) then
        call this%thermal%bc_dir%compute(t)
        call matrix%set_dir_faces(this%thermal%bc_dir%index)
      end if

      Tface_wide(:) = u%tf
      call this%mesh%face_imap%gather_offp(Tface_wide)

      !! Flux-type thermal boundary/interface condition Jacobian contributions.
      call this%thermal%add_flux_bc_jacobian(t, Tface_wide, this%mesh%area, matrix, &
          this%void_face)

      !! Dirichlet fix-ups for void faces.
      if (allocated(this%void_dir_faces)) call matrix%set_dir_faces(this%void_dir_faces)

      !! Enclosure radiation contributions to the preconditioner.
      !! This captures T^4 emitted heat (simple) but ignores absorbed heat
      !! from other faces (complex and non-local).
      if (this%vf_rad%is_active()) then
        call this%vf_rad%add_heat_precon_deriv(t, Tface_wide, this%mesh%area, matrix)
      end if

    end subroutine make_matrix

  end subroutine compute_face_temp

  !! This auxiliary procedure solves for the consistent face concentrations
  !! for one species component, given all cell temperatures and concentrations.

  subroutine compute_face_conc(this, index, t, u, params)

    use mfd_diff_matrix_type
    use pcsr_precon_class
    use pcsr_precon_factory
    use parameter_list_type
    use nka_type
    use parallel_communication, only: global_maxval, global_all

    type(htsd_model), intent(inout), target :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: t
    type(htsd_vector), intent(inout) :: u
    type(parameter_list) :: params

    integer :: stat, max_iter, iter
    type(mfd_diff_matrix), target :: matrix
    real(r8), allocatable :: atol(:), rtol(:), z(:)
    real(r8) :: error0, error, dC_max
    type(htsd_vector) :: udot, f
    type(nka) :: nka_obj
    procedure(pardp), pointer :: dp
    character(80) :: string
    character(:), allocatable :: errmsg
    class(pcsr_precon), allocatable :: precon

    call TLS_info('  Computing consistent face concentrations for species ' // i_to_c(index) // ' ...')

    !! Initial residual of the face concentration equations.
    call udot%init(u)
    call udot%setval(0.0_r8)
    call f%init(u)
    call this%residual(t, u, udot, f)
    associate (r => f%cf(1:this%mesh%nface_onP,index))
      error0 = sqrt(global_dot_product(r, r))
    end associate
    if (error0 == 0.0_r8) return
    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      write(string,'(4x,a,i0,a,es9.2)') '||Rface(', index, ',0)||_2 = ', error0
      call TLS_info(string)
    end if

    !! Form the Jacobian (approx) of the face concentration system.
    call matrix%init(this%disc)
    call make_matrix(matrix)

    call alloc_pcsr_precon(precon, matrix%a22, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_FACE_CONC: ' // errmsg)
    call precon%compute

    call params%get('atol-conc', atol)
    call params%get('rtol-conc', rtol)
    call params%get('max-iter', max_iter)
    error = error0
    dC_max = 0.0_r8

    !! Simple Picard iteration.
    allocate(z(this%mesh%nface_onP))
    call nka_obj%init(size(z), 10)
    dp => pardp
    call nka_obj%set_dot_prod(dp)

    associate (Cface => u%cf(1:this%mesh%nface_onP,index), &
               r => f%cf(1:this%mesh%nface_onP,index))
      iter = 0
      do ! until converged

        iter = iter + 1
        if (iter > max_iter) exit

        z = r
        call precon%apply(z)
        call nka_obj%accel_update(z) ! nonlinear krylov acceleration (NKA)

        Cface = Cface - z

        dC_max = global_maxval(abs(z))

        !! Recompute the face concentration residual.
        call this%residual(t, u, udot, f)
        error = sqrt(global_dot_product(r, r))
        if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
          write(string,'(4x,a,i0,a,i0,2(a,es9.2))') &
              '||Rface(', index, ',', iter, ')||=', error, ', ||ΔCface||_max=', dC_max
          call TLS_info(string)
        end if

        !! Check for convergence.
        if (global_all(abs(z) <= atol(index) + rtol(index) * abs(Cface))) exit
      end do
    end associate

    write(string,'(4x,a,i0,a,i0,3(a,es9.2))') &
        '||Rface(', index, ',', iter, ')||=', error, ', ||Rface(0)||=', error0, &
        ', ||ΔCface||_max=', dC_max
    call TLS_info(string)

    if (iter > max_iter) call TLS_info('    WARNING: face concentrations not converged')

  contains

    subroutine make_matrix(matrix)

      type(mfd_diff_matrix), intent(inout) :: matrix

      real(r8) :: D(this%mesh%ncell), Tface_wide(this%mesh%nface), &
          Cface_wide(this%mesh%nface)

      !! Compute the diffusion coefficient.
      call this%species(index)%diffusivity%compute_value(u%tc, u%cc, D)
      if (allocated(this%void_cell)) where (this%void_cell) D = 0.0_r8

      !! Compute the basic diffusion matrix.
      call matrix%compute(D)

      !! Dirichlet boundary condition fixups.
      if (allocated(this%species(index)%bc_dir)) then
        call this%species(index)%bc_dir%compute(t)
        call matrix%set_dir_faces(this%species(index)%bc_dir%index)
      end if
      if (allocated(this%void_dir_faces)) call matrix%set_dir_faces(this%void_dir_faces)

      Tface_wide(:) = u%tf
      Cface_wide(:) = u%cf(:,index)
      call this%mesh%face_imap%gather_offp(Tface_wide)
      call this%mesh%face_imap%gather_offp(Cface_wide)

      !! Flux-type species boundary/interface condition Jacobian contributions.
      call this%species(index)%add_flux_bc_jacobian(t, Tface_wide, &
          Cface_wide, matrix, this%void_face)

    end subroutine make_matrix

  end subroutine compute_face_conc

  function pardp (a, b) result (dp)
    use parallel_communication, only: global_dot_product
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    dp = global_dot_product(a, b)
  end function pardp

end module htsd_ic_solver_type
