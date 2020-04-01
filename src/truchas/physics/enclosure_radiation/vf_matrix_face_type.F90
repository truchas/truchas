!!
!! VF_MATRIX_FACE_TYPE
!!
!! A concrete implementation for the abstract base class VF_MATRIX.
!! This implementation operates on face-based view factor data.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 18 July 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#include "f90_assert.fpp"

module vf_matrix_face_type

  use kinds, only: r8
  use vf_matrix_class, only: vf_matrix
  use rad_encl_type
  use rad_encl_file_type
  use parallel_communication
  implicit none

  type, extends(vf_matrix), public :: vf_matrix_face
    contains
      procedure, public  :: init => init_VFM_FACE
      procedure, public  :: partition_ER_faces => partition_ER_faces_VFM_FACE
      procedure, public  :: load_view_factors => load_view_factors_VFM_FACE
      procedure, public  :: phi_x => phi_x_VFM_FACE
  end type

contains


  !! Initialize VF_MATRIX_FACE and define a distribution for the matrix rows.
  subroutine init_VFM_FACE (this, file)

    use,intrinsic :: iso_c_binding, only: c_size_t

    class(vf_matrix_face), intent(out) :: this
    type(rad_encl_file), intent(in) :: file

    integer :: nface_tot, npatch_tot
    integer(c_size_t) :: nnonz

    if (is_IOP) call file%get_vf_dims(nface_tot, npatch_tot, nnonz)

    call broadcast(nface_tot)
    this%nface_tot = nface_tot

    !! Partition matrix rows in blocks
    this%nface = this%nface_tot/nPE
    if (this_PE <= modulo(this%nface_tot,nPE)) this%nface = 1 + this%nface

  end subroutine init_VFM_FACE


  !! Partitions the ER mesh according to the distribution of matrix rows
  !! defined by THIS%INIT.
  subroutine partition_ER_faces_VFM_FACE (this, color)

    class(vf_matrix_face), intent(inout) :: this
    integer, allocatable, intent(out) :: color(:)

    integer, allocatable :: color_l(:)
    integer :: n

    !! PARTITION ER FACES. For now we do something really simple and just
    !! equidistribute the faces, dividing up the faces like a salami (no face
    !! reordering).  What we really want to do is partition the faces so that
    !! nonzeros of the (row) distributed view factor matrix are approximately
    !! equidistributed (computational cost) balanced against the communication
    !! cost of moving data between the HC and ER partitions.

    !! Block coloring of the enclosure faces.
    n = merge(this%nface_tot, 0, is_IOP)
    allocate(color(n), color_l(this%nface))
    color_l = this_PE
    call collate (color, color_l)
    deallocate(color_l)

  end subroutine partition_ER_faces_VFM_FACE


  !! Read and distribute the VF matrix data according to the distribution
  !! defined in THIS%INIT.
  subroutine load_view_factors_VFM_FACE (this, file, encl)

    class(vf_matrix_face), intent(inout) :: this
    type(rad_encl_file), intent(in) :: file
    type(rad_encl), intent(in) :: encl

    logical :: renumbered
    integer :: i, offset

    !! The current implementation of THIS%DISTRIBUTE_VF_ROWS doesn't allow
    !! renumbering the enclosure faces, so check that no renumbering has
    !! occurred.  This should be the case if COLOR is a blocked coloring.
    renumbered = .false.
    offset = encl%face_ip%first_index() - 1
    do i = 1, encl%nface_onP
      if (encl%face_map(i) /= i+offset) renumbered = .true.
    end do
    INSIST(.not.global_any(renumbered))

    call this%distribute_vf_rows(file, this%nface, this%nface_tot)

  end subroutine load_view_factors_VFM_FACE


  !! Computes the global product PHI*x, where x is a local vector.
  subroutine phi_x_VFM_FACE (this, lhs, x)

    class(vf_matrix_face), intent(in) :: this
    real(r8), intent(out) :: lhs(:)
    real(r8), intent(in) :: x(:)

    real(r8) :: global_x(this%nface_tot)
    integer :: i, j

    ASSERT(size(lhs) == this%nface)
    ASSERT(size(x) == this%nface)

    call collate (global_x, x)
    call broadcast (global_x)

    lhs = 0.0_r8
    do i = 1, this%nface
      do j = this%ia(i), this%ia(i+1)-1
        lhs(i) = lhs(i) + this%vf(j) * global_x(this%ja(j))
      end do
    end do

  end subroutine phi_x_VFM_FACE


end module vf_matrix_face_type
