!!
!! SD_IC_SOLVER_TYPE
!!
!! This module defines the initial-condition solver used by the species
!! transport solver. It takes user-facing cell concentration initial
!! conditions and computes the full species vector state and initial time
!! derivative needed by the implicit integrator.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The species face concentrations are solved to satisfy their algebraic
!! residual equations. Face concentration time derivatives are computed by
!! finite differencing advanced states with similarly solved face values.
!!

#include "f90_assert.fpp"

module sd_ic_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use sd_model_type
  use sd_vector_type
  use truchas_logging_services
  use parallel_communication
  use parameter_list_type
  use string_utilities, only: i_to_c
  implicit none
  private

  type, public :: sd_ic_solver
    private
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(sd_model), pointer :: model => null() ! unowned reference
    type(parameter_list), pointer :: params => null() ! unowned reference
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_udot
  end type

contains

  subroutine init(this, model, params)
    class(sd_ic_solver), intent(out) :: this
    type(sd_model), intent(in), target :: model
    type(parameter_list), intent(inout), target :: params
    this%model => model
    this%mesh => model%mesh
    this%params => params
  end subroutine

  subroutine compute(this, t, conc, u, udot)

    class(sd_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: conc(:,:)
    type(sd_vector), intent(inout) :: u, udot

    integer :: n
    real(r8), allocatable :: Cface(:)
    type(sd_vector) :: f

    call TLS_info('')
    call TLS_info('Computing consistent initial state for SD solver ...')

    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      ASSERT(size(conc,dim=1) == ncell_onP)
      ASSERT(size(conc,dim=2) == this%model%num_comp)

      call u%setval(0.0_r8)
      do n = 1, this%model%num_comp
        u%cc(:ncell_onP,n) = conc(:,n)
        if (allocated(this%model%void_cell)) &
            where (this%model%void_cell(:ncell_onP)) u%cc(:ncell_onP,n) = 0.0_r8
      end do
      call u%gather_offp

      !! Initial guess for the algebraic face concentration unknowns.
      allocate(Cface(this%mesh%nface))
      do n = 1, this%model%num_comp
        call average_to_faces(this%mesh, u%cc(:,n), Cface, this%model%void_cell)
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

      call TLS_info('')
      call TLS_info('Computing consistent initial state derivative for SD solver ...')

      call udot%setval(0.0_r8)
      call f%init(u)
      call this%model%residual(t, u, udot, f)

      do n = 1, this%model%num_comp
        udot%cc(:ncell_onP,n) = -f%cc(:ncell_onP,n) / this%mesh%volume(:ncell_onP)
        if (allocated(this%model%void_cell)) &
            where (this%model%void_cell(:ncell_onP)) udot%cc(:ncell_onP,n) = 0.0_r8
      end do
      call udot%gather_offp

      call compute_face_conc_deriv(this%model, t, u, udot, this%params)
    end associate

  end subroutine compute


  subroutine compute_udot(this, t, u, udot)

    class(sd_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    type(sd_vector), intent(inout) :: u
    type(sd_vector), intent(inout) :: udot

    integer :: n
    type(sd_vector) :: f

    call TLS_info('')
    call TLS_info('Computing consistent initial state derivative for SD solver ...')

    associate (ncell_onP => this%mesh%ncell_onP)
      call udot%setval(0.0_r8)
      call f%init(u)
      call this%model%residual(t, u, udot, f)

      do n = 1, this%model%num_comp
        udot%cc(:ncell_onP,n) = -f%cc(:ncell_onP,n) / this%mesh%volume(:ncell_onP)
        if (allocated(this%model%void_cell)) &
            where (this%model%void_cell(:ncell_onP)) udot%cc(:ncell_onP,n) = 0.0_r8
      end do
      call udot%gather_offp

      call compute_face_conc_deriv(this%model, t, u, udot, this%params)
    end associate

  end subroutine compute_udot


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

  !! Compute face concentration time derivatives by advancing the cell
  !! concentrations by DT, solving the advanced face concentration equations,
  !! and forming the discrete difference.

  subroutine compute_face_conc_deriv(this, t, u, udot, params)

    use parameter_list_type

    type(sd_model), intent(inout), target :: this
    real(r8), intent(in) :: t
    type(sd_vector), intent(inout) :: u    ! data is intent(in)
    type(sd_vector), intent(inout) :: udot
    type(parameter_list) :: params

    integer :: n
    real(r8) :: dt
    type(sd_vector) :: f

    associate (nface_onP => this%mesh%nface_onP)
      call params%get('dt', dt)
      call f%init(u)
      call f%copy(u)
      call f%update(dt, udot)
      call f%gather_offp

      do n = 1, this%num_comp
        call compute_face_conc(this, n, t+dt, f, params)
      end do

      call f%update(-1.0_r8, u)
      call f%scale(1.0_r8/dt)
      udot%cf(:nface_onP,:) = f%cf(:nface_onP,:)
      call udot%gather_offp
    end associate

  end subroutine compute_face_conc_deriv

  !! This auxiliary procedure solves for the consistent face concentrations
  !! for one species component, given all cell concentrations.

  subroutine compute_face_conc(this, index, t, u, params)

    use mfd_diff_matrix_type
    use pcsr_precon_class
    use pcsr_precon_factory
    use parameter_list_type
    use nka_type
    use parallel_communication, only: global_maxval, global_all

    type(sd_model), intent(inout), target :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: t
    type(sd_vector), intent(inout) :: u
    type(parameter_list) :: params

    integer :: stat, max_iter, iter
    type(mfd_diff_matrix), target :: matrix
    real(r8), allocatable :: atol(:), rtol(:), z(:)
    real(r8) :: error0, error, dC_max
    type(sd_vector) :: udot, f
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

      real(r8) :: D(this%mesh%ncell), Cface_wide(this%mesh%nface)

      !! Compute the diffusion coefficient.
      call this%species(index)%diffusivity%compute_value(u%cc, D)
      if (allocated(this%void_cell)) where (this%void_cell) D = 0.0_r8

      !! Compute the basic diffusion matrix.
      call matrix%compute(D)

      !! Dirichlet boundary condition fixups.
      if (allocated(this%species(index)%bc_dir)) then
        call this%species(index)%bc_dir%compute(t)
        call matrix%set_dir_faces(this%species(index)%bc_dir%index)
      end if
      if (allocated(this%void_dir_faces)) call matrix%set_dir_faces(this%void_dir_faces)

      Cface_wide(:) = u%cf(:,index)
      call this%mesh%face_imap%gather_offp(Cface_wide)

      !! Flux-type species boundary/interface condition Jacobian contributions.
      call this%species(index)%add_flux_bc_jacobian(t, Cface_wide, matrix, this%void_face)

    end subroutine make_matrix

  end subroutine compute_face_conc

  function pardp(a, b) result(dp)
    use parallel_communication, only: global_dot_product
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    dp = global_dot_product(a, b)
  end function pardp

end module sd_ic_solver_type
