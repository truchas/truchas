!!
!! MFD_2D_DIFF_PRECON_TYPE
!!
!! This module implements the Schur-complement preconditioner for the local
!! 2D mimetic finite difference diffusion operator. It owns a frozen-
!! coefficient diffusion matrix and its face Schur-complement preconditioner.
!!
!! David Neill <davidhneill@gmail.com>, April 2020
!! Neil Carlson <neil.n.carlson@gmai.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! APPLY accepts valid on-process cell and face entries. It gathers the
!! off-process workspace needed by the local matrix operations internally.
!!

#include "f90_assert.fpp"

module mfd_2d_diff_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfd_2d_diff_matrix_type
  use pcsr_matrix_type
  use pcsr_precon_class
  use index_map_type
  implicit none
  private

  type, public :: mfd_2d_diff_precon
    private
    type(mfd_2d_diff_matrix), allocatable :: dm
    type(pcsr_matrix) :: Sff
    class(pcsr_precon), allocatable :: Sff_precon
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
    procedure :: matrix_ref
  end type

contains

  subroutine init(this, dm, params, stat, errmsg)

    use pcsr_precon_factory
    use parameter_list_type

    class(mfd_2d_diff_precon), intent(out), target :: this
    type(mfd_2d_diff_matrix), allocatable, intent(inout) :: dm
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    call this%Sff%init(mold=dm%a22)
    call alloc_pcsr_precon(this%Sff_precon, this%Sff, params, stat, errmsg)
    if (stat /= 0) return
    call move_alloc(dm, this%dm)

  end subroutine init

  function matrix_ref(this) result(matrix)
    class(mfd_2d_diff_precon), intent(in), target :: this
    type(mfd_2d_diff_matrix), pointer :: matrix
    matrix => this%dm
  end function

  subroutine compute(this)
    class(mfd_2d_diff_precon), intent(inout) :: this
    call this%dm%compute_face_schur_matrix(this%Sff)
    call this%Sff_precon%compute
  end subroutine

  subroutine apply(this, r1, r2)

    class(mfd_2d_diff_precon), intent(in) :: this
    real(r8), intent(inout) :: r1(:), r2(:)

    ASSERT(size(r1) == this%dm%mesh%ncell)
    ASSERT(size(r2) == this%dm%mesh%nface)

    !! Callers provide valid on-process residuals. The off-process tails are
    !! workspace and are synchronized here before operations that read
    !! neighboring cell or face values.
    call this%dm%mesh%cell_imap%gather_offp(r1)
    call this%dm%mesh%face_imap%gather_offp(r2)

    !! Eliminate the cell unknowns.
    call forward_elimination(this%dm, r1, r2)

    !! Approximately solve the Schur complement system for the face unknowns.
    call this%Sff_precon%apply(r2)
    call this%dm%mesh%face_imap%gather_offp(r2)

    !! Solve for the cell unknowns by back substitution.
    call backward_substitution(this%dm, r1, r2)

  end subroutine apply

  subroutine forward_elimination(dm, b1, b2)

    type(mfd_2d_diff_matrix), intent(in) :: dm
    real(r8), intent(in) :: b1(:)
    real(r8), intent(inout) :: b2(:)

    integer :: j
    real(r8) :: s
    real(r8), allocatable :: b2_dir(:)

    ASSERT(size(b1) == dm%mesh%ncell)
    ASSERT(size(b2) == dm%mesh%nface)

    if (allocated(dm%dir_faces)) then
      allocate(b2_dir(size(dm%dir_faces)))
      b2_dir = b2(dm%dir_faces)
    end if

    do j = 1, dm%mesh%ncell
      s = b1(j) / dm%a11(j)
      associate (cface => dm%mesh%cface(dm%mesh%cstart(j):dm%mesh%cstart(j+1)-1), &
                 a12 => dm%a12_val(dm%mesh%cstart(j):dm%mesh%cstart(j+1)-1))
        b2(cface) = b2(cface) - a12 * s
      end associate
    end do

    if (allocated(dm%dir_faces)) then
      b2(dm%dir_faces) = b2_dir
      deallocate(b2_dir)
    end if

  end subroutine forward_elimination


  subroutine backward_substitution(dm, b1, u2)

    type(mfd_2d_diff_matrix), intent(in) :: dm
    real(r8), intent(inout) :: b1(:), u2(:)

    integer :: j
    real(r8), allocatable :: u2_dir(:)

    if (allocated(dm%dir_faces)) then
      allocate(u2_dir(size(dm%dir_faces)))
      u2_dir = u2(dm%dir_faces)
      u2(dm%dir_faces) = 0.0_r8
    end if

    do j = 1, dm%mesh%ncell_onP
      associate (cface => dm%mesh%cface(dm%mesh%cstart(j):dm%mesh%cstart(j+1)-1), &
                 a12 => dm%a12_val(dm%mesh%cstart(j):dm%mesh%cstart(j+1)-1))
        b1(j) = (b1(j) - dot_product(a12, u2(cface))) / dm%a11(j)
      end associate
    end do

    if (allocated(dm%dir_faces)) then
      u2(dm%dir_faces) = u2_dir
      deallocate(u2_dir)
    end if

  end subroutine backward_substitution

end module mfd_2d_diff_precon_type
