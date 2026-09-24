!!
!! SD_PRECON_TYPE
!!
!! This module defines the preconditioner used by the implicit integrator for
!! the nonlinear species transport time-step system. It assembles and applies
!! approximate Jacobians for the species transport residuals, including species
!! boundary contributions.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module sd_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use sd_model_type
  use sd_vector_type
  use unstr_mesh_type
  use mfd_diff_precon_type
  use mfd_diff_matrix_type
  implicit none
  private

  type, public :: sd_precon
    type(sd_model), pointer :: model => null() ! unowned reference
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(mfd_diff_precon), allocatable :: pc(:)
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(sd_precon), intent(out) :: this
    type(sd_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    type(mfd_diff_matrix), allocatable, target :: matrix
    type(mfd_diff_matrix), pointer :: mold_matrix => null()

    this%model => model
    this%mesh => model%mesh

    !! Create the species transport preconditioner; compute defines its values.
    allocate(this%pc(model%num_comp))
    do n = 1, model%num_comp
      allocate(matrix)
      if (associated(mold_matrix)) then
        call matrix%init(mold=mold_matrix)
      else
        call matrix%init(model%disc)
      end if
      call this%pc(n)%init(matrix, params, stat, errmsg)
      if (stat /= 0) return
      if (.not.associated(mold_matrix)) mold_matrix => this%pc(n)%matrix_ref()
    end do

  end subroutine init

  subroutine compute(this, t, u, dt)

    class(sd_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    type(sd_vector), intent(inout) :: u

    integer :: n
    real(r8) :: values(this%mesh%ncell)

    ASSERT(dt > 0.0_r8)

    !! Preconditioner assembly evaluates material properties over the
    !! off-process extended concentration fields.
    call u%gather_offp

    !! Compute the species transport preconditioner.
    do n = 1, size(this%pc)
      call compute_species_precon(n)
    end do

  contains

    subroutine compute_species_precon(index)

      integer, intent(in) :: index

      type(mfd_diff_matrix), pointer :: matrix

      matrix => this%pc(index)%matrix_ref()

      !! Jacobian of the diffusion operator that ignores nonlinearities.
      call this%model%species(index)%diffusivity%compute_value(u%cc, values)
      if (allocated(this%model%void_cell)) where (this%model%void_cell) values = 0.0_r8
      call matrix%compute(values)

      !! Time derivative contribution to the diffusion equation Jacobian.
      values = this%mesh%volume / dt
      if (allocated(this%model%void_cell)) where (this%model%void_cell) values = 1.0_r8
      call matrix%incr_cell_diag(values)

      !! Dirichlet BC fixups.
      if (allocated(this%model%species(index)%bc_dir)) then
        call this%model%species(index)%bc_dir%compute(t)
        call matrix%set_dir_faces(this%model%species(index)%bc_dir%index)
      end if
      if (allocated(this%model%void_dir_faces)) call matrix%set_dir_faces(this%model%void_dir_faces)

      !! Flux-type species boundary/interface condition Jacobian contributions.
      call this%model%species(index)%add_flux_bc_jacobian(t, u%cf(:,index), &
          matrix, this%model%void_face)

      !! The matrix is now complete; re-compute the preconditioner.
      call this%pc(index)%compute

    end subroutine compute_species_precon

  end subroutine compute

  subroutine apply(this, t, u, f)
    class(sd_precon), intent(in) :: this
    real(r8), intent(in) :: t
    type(sd_vector), intent(inout) :: u  ! data is intent(in)
    type(sd_vector), intent(inout) :: f
    integer :: index
    do index = 1, this%model%num_comp
      !! The component preconditioner owns any off-process gathers it needs;
      !! only on-process entries are significant on return.
      call this%pc(index)%apply(f%cc(:,index), f%cf(:,index))
    end do
  end subroutine

end module sd_precon_type
