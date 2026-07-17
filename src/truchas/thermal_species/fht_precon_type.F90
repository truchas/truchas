!!
!! FHT_PRECON_TYPE
!!
!! This module defines the preconditioner used by the nonadaptive integrator
!! for the nonlinear flow-coupled thermal transport time-step system. It
!! assembles and applies an approximate Jacobian for the cell and face
!! temperature residuals, including thermal boundary contributions and optional
!! view-factor radiation coupling.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module fht_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fht_model_type
  use fht_vector_type
  use unstr_mesh_type
  use mfd_diff_matrix_type
  use mfd_diff_precon_type
  use thermal_view_factor_coupling_type, only: VFR_JAC, VFR_FGS, VFR_BGS, VFR_FAC, &
      vfr_precon_coupling_from_string
  use truchas_timers
  implicit none
  private
  
  type, public :: fht_precon
    type(fht_model), pointer :: model => null() ! unowned reference
    type(unstr_mesh), pointer :: mesh  => null() ! unowned reference
    integer :: vfr_precon_coupling = VFR_BGS
    type(mfd_diff_precon) :: pc
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type
  
contains

  subroutine init(this, model, params, stat, errmsg)
  
    use parameter_list_type

    class(fht_precon), intent(out) :: this
    type(fht_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    
    character(:), allocatable :: method
    type(mfd_diff_matrix), allocatable :: matrix

    this%model => model
    this%mesh  => model%mesh
    
    !! Create the preconditioner for the heat equation.
    allocate(matrix)
    call matrix%init(model%disc)
    call this%pc%init(matrix, params, stat, errmsg)
    if (stat /= 0) return
    
    !! Heat equation/view factor radiation preconditioner coupling method.
    if (model%vf_rad%is_active()) then
      call params%get('vfr-precon-coupling', method, default='backward-gs')
      this%vfr_precon_coupling = vfr_precon_coupling_from_string(method)
    end if

  end subroutine init
  
  subroutine apply(this, t, u, f)

    class(fht_precon), intent(inout) :: this
    real(r8), intent(in) :: t
    type(fht_vector), intent(in) :: u
    type(fht_vector), intent(inout) :: f

    real(r8), allocatable :: z(:)

    call start_timer('precon apply')

    !! Precondition the radiosity components:
    !! Factorization and forward Gauss-Seidel coupling methods.
    if (this%model%vf_rad%is_active()) then
      if (this%vfr_precon_coupling == VFR_FGS .or. &
          this%vfr_precon_coupling == VFR_FAC) then
        z = f%qrad
        call this%model%vf_rad%apply_rad_precon(t, z)
        if (this%vfr_precon_coupling == VFR_FGS) f%qrad(:) = z
        !! Update the heat equation face residual.
        call this%model%vf_rad%apply_rad_precon_matvec1(t, z)
        call this%model%vf_rad%add_heat_flux_to_residual(this%mesh%area, z, f%tf)
        deallocate(z)
      end if
    end if

    !! Precondition the heat equation. The diffusion preconditioner owns any
    !! off-process gathers it needs; only on-process results are significant
    !! after the call.
    call this%pc%apply(f%tc, f%tf)

    !! Precondition the radiosity components:
    !! Factorization, backward Gauss-Seidel and Jacobi coupling methods.
    if (this%model%vf_rad%is_active()) then
      if (this%vfr_precon_coupling == VFR_JAC .or. &
          this%vfr_precon_coupling == VFR_BGS .or. &
          this%vfr_precon_coupling == VFR_FAC) then
        !! Update the radiosity residual components.
        if (this%vfr_precon_coupling /= VFR_JAC) then
          call this%model%vf_rad%add_rhs_deriv_times_face_vector(t, u%tf, f%tf, f%qrad)
        end if
        call this%model%vf_rad%apply_rad_precon(t, f%qrad)
      end if
    end if

    call stop_timer('precon apply')

  end subroutine apply

  subroutine compute(this, t, u, h)

    class(fht_precon), intent(inout) :: this
    real(r8), intent(in) :: t, h
    type(fht_vector), intent(inout) :: u ! data is intent(in)

    integer :: n, j
    real(r8) :: D(this%mesh%ncell), A(this%mesh%ncell), dHdT(this%mesh%ncell)
    type(mfd_diff_matrix), pointer :: matrix
    integer, allocatable :: more_dir_faces(:)

    ASSERT(h > 0.0_r8)
    ASSERT(associated(this%model%void_cell))
    ASSERT(associated(this%model%void_face))

    call start_timer('precon compute')

    !! Preconditioner assembly evaluates material properties and boundary
    !! condition derivatives over the off-process extended temperature field.
    call u%gather_offp

    !! Derivative with respect to temperature, the first material-function argument.
    call this%model%thermal%H_of_T%compute_deriv(u%tc, 1, dHdT)
    call this%model%thermal%conductivity%compute_value(u%tc, D)
    A = this%mesh%volume * dHdT / h

    !! Correct data on void cells.
    where (this%model%void_cell)
      dHdT = 0.0_r8
      D = 0.0_r8
      A = 1.0_r8
    end where

    matrix => this%pc%matrix_ref()

    !! Jacobian of the basic heat equation that ignores nonlinearities
    !! in the conductivity.  This has the algebraic H/T relation eliminated.
    call matrix%compute(D)
    call matrix%incr_cell_diag(A)

    !! Dirichlet boundary condition fixups.
    if (allocated(this%model%thermal%bc_dir)) then
      call this%model%thermal%bc_dir%compute(t)
      call matrix%set_dir_faces(this%model%thermal%bc_dir%index)
    end if

    !! Flux-type thermal boundary/interface condition Jacobian contributions.
    call this%model%thermal%add_flux_bc_jacobian(t, u%tf, this%mesh%area, matrix, &
        this%model%void_face)

    n = count(this%model%void_face)
    if (n > 0) then
      allocate(more_dir_faces(n))
      n = 0
      do j = 1, this%mesh%nface
        if (this%model%void_face(j)) then
          n = n + 1
          more_dir_faces(n) = j
        end if
      end do
      call matrix%set_dir_faces(more_dir_faces)
      deallocate(more_dir_faces)
    end if

    !! Heat-side view-factor radiation contribution to the preconditioner.
    !! The block coupling method is handled in apply.
    if (this%model%vf_rad%is_active()) then
      call this%model%vf_rad%add_heat_precon_deriv(t, u%tf, this%mesh%area, matrix)
    end if

    !! The matrix is now complete; re-compute the preconditioner.
    call this%pc%compute

    call stop_timer('precon compute')

  end subroutine compute

end module fht_precon_type
