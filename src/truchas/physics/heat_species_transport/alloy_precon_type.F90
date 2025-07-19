!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module alloy_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use alloy_vector_type
  use alloy_model_type
  use unstr_mesh_type
  use mfd_diff_precon_type
  use mfd_diff_matrix_type
  use truchas_timers
  implicit none
  private

  type, public :: alloy_precon
    type(alloy_model), pointer :: model => null()
    type(unstr_mesh), pointer :: mesh  => null()
    real(r8) :: dt
    real(r8), allocatable :: dHdT(:), dHdg(:), dgdT(:)
    real(r8), allocatable :: drdg(:), drdH(:), b(:)
    type(mfd_diff_precon) :: pc
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type alloy_precon

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(alloy_precon), intent(out) :: this
    type(alloy_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(mfd_diff_matrix), allocatable :: matrix

    this%model => model
    this%mesh  => model%mesh

    allocate(this%dHdT(this%mesh%ncell), this%dHdg(this%mesh%ncell), this%dgdT(this%mesh%ncell))
    allocate(this%drdg(this%mesh%ncell), this%drdH(this%mesh%ncell), this%b(this%mesh%ncell))
    allocate(matrix)
    call matrix%init(model%disc)
    call this%pc%init(matrix, params, stat, errmsg)

  end subroutine init


  subroutine compute(this, t, u, dt)

    class(alloy_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    type(alloy_vector), intent(inout) :: u
    target :: u

    real(r8), pointer :: state(:,:) => null()

    integer :: index, j, n, n1, n2
    real(r8) :: D(this%mesh%ncell), A(this%mesh%ncell), term
    real(r8), allocatable :: values(:), values2(:,:)
    type(mfd_diff_matrix), pointer :: dm

    ASSERT(dt > 0.0_r8)

    call start_timer('alloy-precon-compute')

    !TODO: The existing prop_mesh_func%compute_value function expects a rank-2
    !      state array. This is a workaround until prop_mesh_func is redesigned.
    state(1:this%mesh%ncell,1:1) => u%tc

    !call this%mesh%cell_imap%gather_offp(u%lf)
    !call this%mesh%cell_imap%gather_offp(u%hc)
    !call this%mesh%cell_imap%gather_offp(u%tc)
    !call this%mesh%face_imap%gather_offp(u%tf)
    call u%gather_offp

    call this%model%alloy%compute_dHdg(u%tc, u%lf, this%dHdg)
    call this%model%alloy%compute_dHdT(u%tc, u%lf, this%dHdT)
    call this%model%alloy%compute_drdg(u%hc, u%lf, this%drdg)
    call this%model%alloy%compute_drdH(u%hc, u%lf, this%drdH)

    this%B = 1 + this%dHdg*this%drdH/this%drdg

    this%dt = dt

    !! Hardwired assumption that T is the first component of state -- FIXME!
    !call this%model%H_of_T%compute_deriv(state, 1, this%dHdT)
    call this%model%conductivity%compute_value(state, D)
    A = this%mesh%volume * (this%dHdT/this%B) / dt

    !! Jacobian of the basic heat equation that ignores nonlinearities
    !! in the conductivity.  This has the H/T relation eliminated.
    dm => this%pc%matrix_ref()
    call dm%compute(D)
    call dm%incr_cell_diag(A)

    !! Dirichlet boundary condition fixups.
    if (allocated(this%model%bc_dir)) then
      call this%model%bc_dir%compute(t)
      call dm%set_dir_faces(this%model%bc_dir%index)
    end if

    !! External HTC boundary condition contribution.
    if (allocated(this%model%bc_htc)) then
      call this%model%bc_htc%compute_deriv(t, u%tf)
      associate (index => this%model%bc_htc%index, &
                 deriv => this%model%bc_htc%deriv)
        call dm%incr_face_diag(index, deriv)
      end associate
    end if

    !! Simple radiation boundary condition contribution.
    if (allocated(this%model%bc_rad)) then
      call this%model%bc_rad%compute(t, u%tf)
      associate (index => this%model%bc_rad%index, &
                 deriv => this%model%bc_rad%deriv)
        call dm%incr_face_diag(index, deriv)
      end associate
    end if

    !! Experimental evaporation heat flux contribution.
    if (allocated(this%model%evap_flux)) then
      call this%model%evap_flux%compute_deriv(t, u%tf)
      associate (index => this%model%evap_flux%index, &
                 deriv => this%model%evap_flux%deriv)
        call dm%incr_face_diag(index, this%mesh%area(index)*deriv)
      end associate
    endif

    !! Internal HTC interface condition contribution.
    if (allocated(this%model%ic_htc)) then
      call this%model%ic_htc%compute(t, u%tf)
      associate (index => this%model%ic_htc%index, &
                 deriv => this%model%ic_htc%deriv)
        call dm%incr_interface_flux3(index, deriv) !TODO: rename these methods
      end associate
    end if

    !! Internal gap radiation condition contribution.
    if (allocated(this%model%ic_rad)) then
      call this%model%ic_rad%compute(t, u%tf)
      associate (index => this%model%ic_rad%index, &
                 deriv => this%model%ic_rad%deriv)
        call dm%incr_interface_flux3(index, deriv) !TODO: rename these methods
      end associate
    end if

    !! The matrix is now complete; re-compute the preconditioner.
    call this%pc%compute

    call stop_timer('alloy-precon-compute')

  end subroutine compute


  subroutine apply(this, t, u, f)

    class(alloy_precon), intent(in) :: this
    real(r8), intent(in) :: t
    type(alloy_vector), intent(inout) :: u   ! data is intent(in)
    type(alloy_vector), intent(inout) :: f   ! data is intent(inout)

    integer :: index, j, n
    real(r8), allocatable :: z(:)

    call start_timer('alloy-precon-apply')

    ! Elimination of algebraic constraints
    f%hc = f%hc + (this%dHdg/this%drdg) * f%lf
    f%tc = f%tc - (this%mesh%volume/this%dt)*f%hc/this%B

    ! Precondition the heat equation
    call this%mesh%cell_imap%gather_offp(f%tc)
    call this%mesh%face_imap%gather_offp(f%tf)
    call this%pc%apply(f%tc, f%tf)

    !! Back-substitute to precondition the algebraic constraints
    f%hc = (f%hc + this%dHdT*f%tc)/this%B
    f%lf = (f%lf - this%drdH*f%hc)/this%drdg

    call stop_timer('alloy-precon-apply')

  end subroutine apply

end module alloy_precon_type
