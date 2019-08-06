!!
!! provides a derived type for packing together cell geometry data
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module cell_geom_2d_vof_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  integer, parameter, public :: &
      CELL_TRI3 = 1, CELL_QUA4 = 2

  type, public :: cell_geom
    real(r8), allocatable :: node(:,:), face_area(:), face_normal(:,:)
    real(r8) :: volume
    integer :: nfc, cell_type
  contains
    procedure :: init
  end type cell_geom

contains

  subroutine init (this, node, volume, face_area, face_normal)

    use truchas_logging_services

    class(cell_geom), intent(out) :: this
    real(r8), intent(in) :: node(:,:), volume, face_area(:), face_normal(:,:)

    this%node = node
    this%face_area = face_area
    this%face_normal = face_normal

    this%volume = volume

    ! get the number of faces for this cell type
    select case (size(node, dim=2))
    case (3) ! triangle
      this%nfc = 3
      this%cell_type = CELL_TRI3
    case (4) ! quadrilateral
      this%nfc = 4
      this%cell_type = CELL_QUA4
    case default
      INSIST(.false.)
    end select

  end subroutine init

end module cell_geom_2d_vof_type
