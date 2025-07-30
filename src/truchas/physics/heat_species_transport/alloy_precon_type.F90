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
  use alloy_back_diff_type, only: alloy_back_diff_jac
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
    type(alloy_back_diff_jac), allocatable :: jac(:)
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

    select case (this%model%model_type)
    case (2) ! Wang-Beckermann
      block !FIXME: THIS DATA STRUCTURE IS AWFUL!
        integer :: j
        allocate(this%jac(this%mesh%ncell))
        do j = 1, this%mesh%ncell
          call this%jac(j)%init(this%model%pd%num_comp)
        end do
      end block
    end select

  end subroutine init


  subroutine compute(this, t, u, udot, dt)

    class(alloy_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    type(alloy_vector), intent(inout) :: u, udot
    target :: u

    integer :: j
    real(r8) :: D(this%mesh%ncell), A(this%mesh%ncell)
    type(mfd_diff_matrix), pointer :: dm

    ASSERT(dt > 0.0_r8)

    call start_timer('alloy-precon-compute')

    call u%gather_offp
    call udot%gather_offp

    select case (this%model%model_type)
    case (1) ! lever rule
      call this%model%alloy%compute_g_jac(this%model%C, u%lf, u%hc, this%drdg, this%drdH) !TODO: rename result arrays
      call this%model%alloy%compute_H_jac(u%lf, u%hc, u%tc, this%dHdg, this%B, this%dHdT)
      this%B = this%B - this%dHdg*this%drdH/this%drdg
      A = -this%mesh%volume * (this%dHdT/this%B) / dt
    case (2) ! Wang-Beckermann
      do j = 1, this%mesh%ncell
        call this%model%pd%compute_f_jac(this%model%C(:,j), this%model%Cdot(:,j), &
            u%lsf(:,j), u%lf(j), u%hc(j), u%tc(j), udot%lsf(:,j), udot%lf(j), dt, this%jac(j))
        call this%jac(j)%lu_factor
        this%B(j) = (1/dt)*this%mesh%volume(j)/this%jac(j)%dfHdH
        A(j) = -this%B(j)*this%jac(j)%dfHdT
      end do
    end select

    this%dt = dt

    call this%model%get_conductivity(u%lf, u%tc, D)

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

    !! The matrix is now complete; re-compute the preconditioner.
    call this%pc%compute

    call stop_timer('alloy-precon-compute')

  end subroutine compute


  subroutine apply(this, t, u, f)

    class(alloy_precon), intent(in) :: this
    real(r8), intent(in) :: t
    type(alloy_vector), intent(inout) :: u   ! data is intent(in)
    type(alloy_vector), intent(inout) :: f   ! data is intent(inout)

    integer :: j

    call start_timer('alloy-precon-apply')

    ! Elimination of algebraic constraints
    select case (this%model%model_type)
    case (1) ! lever rule
      f%hc = f%hc - (this%dHdg/this%drdg) * f%lf
      f%tc = f%tc - (this%mesh%volume/this%dt)*f%hc/this%B
    case (2) ! Wang-Beckermann
      do j = 1, this%mesh%ncell
        call this%jac(j)%lower_solve(f%lsf(:,j), f%lf(j), f%hc(j))
        f%tc(j) = f%tc(j) - this%B(j) * f%hc(j)
      end do
    end select

    ! Precondition the heat equation
    call this%mesh%cell_imap%gather_offp(f%tc)
    call this%mesh%face_imap%gather_offp(f%tf)
    call this%pc%apply(f%tc, f%tf)

    !! Back-substitute to precondition the algebraic constraints
    select case (this%model%model_type)
    case (1) ! lever rule
      f%hc = (f%hc - this%dHdT*f%tc)/this%B
      f%lf = (f%lf - this%drdH*f%hc)/this%drdg
    case (2) ! Wang-Beckermann
      do j = 1, this%mesh%ncell
        call this%jac(j)%upper_solve(f%lsf(:,j), f%lf(j), f%hc(j), f%tc(j))
      end do
    end select

    call stop_timer('alloy-precon-apply')

  end subroutine apply

end module alloy_precon_type
