!!
!! FACET_VF_TYPE
!!
!! A concrete implementation for the abstract base class STATIC_VF.
!! This implementation operates on enclosure mesh face-based view factor data.
!!
!! David Neill-Asanza <dhna@lanl.gov>, July 2019
!! Neil N. Carlson <nnc@lanl.gov>, refactoring June 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module facet_vf_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use static_vf_class
  use vf_data_type
  use parallel_communication
  implicit none
  private

  type, extends(static_vf), public :: facet_vf
    private
    type(vf_data) :: vf
    ! In parent classes:
    ! integer, public :: nface, nface_tot
    ! integer, allocatable :: face_map(:)
    ! real, allocatable, public :: amb_vf(:)
    ! logical :: has_ambient
  contains
    procedure :: init
    procedure :: matvec
  end type

contains

  subroutine init(this, file)

    use rad_encl_file_type

    class(facet_vf), intent(out) :: this
    type(rad_encl_file), intent(in) :: file  ! only referenced on IOP

    integer :: j, n, offset, bsize(nPE)

    call this%vf%init(file)

    if (is_IOP) call file%get_patch_dims(this%nface_tot, n)
    call broadcast(this%nface_tot)

    !! Distributed face renumbering (new-to-old face map) subordinate to the
    !! block-row partitioning of the view factor matrix. This is the identity!
    this%nface = this%vf%npatch
    allocate(this%face_map(this%nface))
    call collate(bsize, this%nface)
    call broadcast(bsize)
    offset = sum(bsize(:this_PE-1))
    do j = 1, this%nface
      this%face_map(j) = offset + j
    end do

    !! Move the ambient view factors to the public vector
    this%has_ambient = allocated(this%vf%amb_vf)
    if (this%has_ambient) then
      call move_alloc(this%vf%amb_vf, this%amb_vf)
    else
      allocate(this%amb_vf(this%nface))
      this%amb_vf = 0
    end if

  end subroutine init

  !! Computes the matrix-vector product Y = PHI*X.
  subroutine matvec(this, x, y)
    class(facet_vf), intent(in) :: this
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: y(:)
    call this%vf%matvec(x, y)
  end subroutine

end module facet_vf_type
