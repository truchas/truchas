!!
!! HT_2D_PRECON_TYPE
!!
!! This module defines the preconditioner used by the implicit 2D thermal
!! transport solver. It assembles and applies an approximate Jacobian for the
!! coupled enthalpy-temperature residual.
!!
!! David Neill-Asanza <davidhneill@gmail.com>, April 2020
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! Property evaluation is restricted to on-process cells because material
!! composition is stored there; the resulting coefficient fields are gathered
!! before diffusion-matrix assembly.
!!

#include "f90_assert.fpp"

module ht_2d_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ht_2d_model_type
  use ht_2d_vector_type
  use unstr_2d_mesh_type
  use mfd_2d_diff_precon_type
  use mfd_2d_diff_matrix_type
  implicit none
  private

  type, public :: ht_2d_precon
    type(ht_2d_model),   pointer :: model => null() ! unowned reference
    type(unstr_2d_mesh), pointer :: mesh  => null() ! unowned reference
    real(r8) :: dt ! time step
    real(r8), allocatable :: dHdT(:) ! derivative of the enthalpy/temperature relation
    type(mfd_2d_diff_precon) :: pc   ! heat equation preconditioner
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(ht_2d_precon), intent(out), target :: this
    type(ht_2d_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(mfd_2d_diff_matrix), allocatable :: dm

    this%model => model
    this%mesh  => model%mesh

    !! Heat density/temperature relation derivative.
    allocate(this%dHdT(this%mesh%ncell))

    !! Create the preconditioner for the heat equation.
    !! The preconditioner assumes ownership of the matrix.
    allocate(dm)
    call dm%init(model%disc)
    call this%pc%init(dm, params, stat, errmsg)

  end subroutine init


  subroutine compute(this, t, u, dt)

    class(ht_2d_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    type(ht_2d_vector), intent(inout) :: u

    real(r8) :: coef(this%mesh%ncell)
    type(mfd_2d_diff_matrix), pointer :: dm

    ASSERT(dt > 0.0_r8)

    this%dt = dt
    dm => this%pc%matrix_ref()
    call this%mesh%face_imap%gather_offp(u%tf)
    call this%model%conductivity%compute_value(u%tc, coef(:this%mesh%ncell_onP))
    call this%mesh%cell_imap%gather_offp(coef)
    call dm%compute(coef)
    call this%model%H_of_T%compute_deriv(u%tc, 1, this%dHdT(:this%mesh%ncell_onP))
    call this%mesh%cell_imap%gather_offp(this%dHdT)
    call dm%incr_cell_diag(this%mesh%volume*this%dHdT/dt)

    if (allocated(this%model%bc_dir)) then
      call this%model%bc_dir%compute(t)
      call dm%set_dir_faces(this%model%bc_dir%index)
    end if
    if (allocated(this%model%bc_htc)) then
      call this%model%bc_htc%compute(t, u%tf)
      call dm%incr_face_diag(this%model%bc_htc%index, this%model%bc_htc%deriv)
    end if
    if (allocated(this%model%bc_rad)) then
      call this%model%bc_rad%compute(t, u%tf)
      call dm%incr_face_diag(this%model%bc_rad%index, this%model%bc_rad%deriv)
    end if
    call this%pc%compute

  end subroutine compute


  subroutine apply(this, t, u, r)

    class(ht_2d_precon), intent(in) :: this
    real(r8), intent(in) :: t
    type(ht_2d_vector), intent(inout) :: u
    type(ht_2d_vector), intent(inout) :: r

    associate (mesh => this%mesh)
      !! Eliminate the enthalpy-temperature residual from the heat equation.
      r%tc(:mesh%ncell_onP) = r%tc(:mesh%ncell_onP) &
                             - (mesh%volume(:mesh%ncell_onP)/this%dt)*r%hc(:mesh%ncell_onP)
      call this%pc%apply(r%tc, r%tf)

      !! Backsubstitute the enthalpy-temperature relation residual.
      r%hc(:mesh%ncell_onP) = r%hc(:mesh%ncell_onP) &
                             + this%dHdT(:mesh%ncell_onP)*r%tc(:mesh%ncell_onP)
    end associate

  end subroutine apply

end module ht_2d_precon_type
