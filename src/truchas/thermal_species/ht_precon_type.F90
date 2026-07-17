!!
!! HT_PRECON_TYPE
!!
!! This module defines the preconditioner used by the implicit integrator for
!! the nonlinear thermal transport time-step system. It assembles and applies
!! an approximate Jacobian for the coupled enthalpy-temperature residual,
!! including thermal boundary contributions and optional view-factor radiation
!! coupling.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ht_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ht_vector_type
  use ht_model_type
  use unstr_mesh_type
  use mfd_diff_precon_type
  use mfd_diff_matrix_type
  use thermal_view_factor_coupling_type, only: VFR_JAC, VFR_FGS, VFR_BGS, VFR_FAC
  implicit none
  private

  type, public :: ht_precon
    type(ht_model), pointer :: model => null() ! unowned reference
    type(unstr_mesh), pointer :: mesh  => null() ! unowned reference
    real(r8) :: dt
    real(r8), allocatable :: dHdT(:)
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
    use thermal_view_factor_coupling_type, only: vfr_precon_coupling_from_string

    class(ht_precon), intent(out) :: this
    type(ht_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(mfd_diff_matrix), allocatable :: matrix
    character(:), allocatable :: method

    this%model => model
    this%mesh  => model%mesh

    !! Enthalpy density/temperature relation derivative.
    allocate(this%dHdT(this%mesh%ncell))

    !! Create the heat equation preconditioner; compute defines its values.
    allocate(matrix)
    call matrix%init(model%disc)
    call this%pc%init(matrix, params, stat, errmsg)

    !! Heat equation/view factor radiation preconditioner coupling method.
    if (model%vf_rad%is_active()) then
      call params%get('vfr-precon-coupling', method, default='backward-gs')
      this%vfr_precon_coupling = vfr_precon_coupling_from_string(method)
    end if

  end subroutine init


  subroutine compute(this, t, u, dt)

    class(ht_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    type(ht_vector), intent(inout) :: u

    real(r8) :: D(this%mesh%ncell), A(this%mesh%ncell)
    type(mfd_diff_matrix), pointer :: matrix

    ASSERT(dt > 0.0_r8)

    this%dt = dt

    !! Preconditioner assembly evaluates material properties and boundary
    !! condition derivatives over the off-process extended temperature field.
    call u%gather_offp

    !! Derivative with respect to temperature, the first material-function argument.
    call this%model%thermal%H_of_T%compute_deriv(u%tc, 1, this%dHdT)
    call this%model%thermal%conductivity%compute_value(u%tc, D)
    A = this%mesh%volume * this%dHdT / dt

    !! Correct data on void cells.
    if (allocated(this%model%void_cell)) then
      where (this%model%void_cell)
        this%dHdT = 0.0_r8
        D = 0.0_r8
        A = 1.0_r8
      end where
    end if

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

    !! Dirichlet fix-ups for void faces.
    if (allocated(this%model%void_dir_faces)) call matrix%set_dir_faces(this%model%void_dir_faces)

    !! Enclosure radiation contributions to the preconditioner.
    !! TODO: what about factorization coupling?  Is this still correct?
    if (this%model%vf_rad%is_active()) then
      call this%model%vf_rad%add_heat_precon_deriv(t, u%tf, this%mesh%area, matrix)
    end if

    !! The matrix is now complete; re-compute the preconditioner.
    call this%pc%compute

  end subroutine compute


  subroutine apply(this, t, u, f)

    class(ht_precon), intent(in) :: this
    real(r8), intent(in) :: t
    type(ht_vector), intent(inout) :: u   ! data is intent(in)
    type(ht_vector), intent(inout) :: f   ! data is intent(inout)

    integer :: ncell_onP
    real(r8), allocatable :: z(:)

    ncell_onP = this%mesh%ncell_onP

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

    !! Heat equation cell residual with the H/T relation residual eliminated.
    if (allocated(this%model%void_cell)) then
      where (.not.this%model%void_cell(:ncell_onP))
        f%tc(:ncell_onP) = f%tc(:ncell_onP) &
            - (this%mesh%volume(:ncell_onP)/this%dt)*f%hc(:ncell_onP)
      end where
    else
      f%tc(:ncell_onP) = f%tc(:ncell_onP) &
          - (this%mesh%volume(:ncell_onP)/this%dt)*f%hc(:ncell_onP)
    end if

    !! Precondition the heat equation. The diffusion preconditioner owns any
    !! off-process gathers it needs; only on-process results are significant
    !! after the call.
    call this%pc%apply(f%tc, f%tf)

    !! Backsubstitute to obtain the preconditioned H/T-relation residual.
    f%hc(:ncell_onP) = f%hc(:ncell_onP) + this%dHdT(:ncell_onP)*f%tc(:ncell_onP)

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

  end subroutine apply

end module ht_precon_type
