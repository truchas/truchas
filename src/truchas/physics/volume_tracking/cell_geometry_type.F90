!!
!! provides a derived type for packing together cell geometry data
!!

#include "f90_assert.fpp"

module cell_geom_type

  use kinds, only: r8
  implicit none
  private

  type, public :: cell_geom
    real(r8) :: node(3,8), volume
    real(r8) :: face_area(6), face_normal(3,6)
  contains
    procedure :: init
  end type cell_geom

  ! ! cube vertex positions for unit testing
  ! real(r8), parameter, public :: cube_v(3,8) = reshape([ &
  !      0.0_r8, 0.0_r8, 0.0_r8, &
  !      1.0_r8, 0.0_r8, 0.0_r8, &
  !      1.0_r8, 1.0_r8, 0.0_r8, &
  !      0.0_r8, 1.0_r8, 0.0_r8, &
  !      0.0_r8, 0.0_r8, 1.0_r8, &
  !      1.0_r8, 0.0_r8, 1.0_r8, &
  !      1.0_r8, 1.0_r8, 1.0_r8, &
  !      0.0_r8, 1.0_r8, 1.0_r8],&
  !      shape(cube_v))

contains

  subroutine init (this, node, volume, face_area, face_normal)

    class(cell_geom), intent(out) :: this
    real(r8), intent(in) :: node(:,:), volume, face_area(:), face_normal(:,:)

    this%node = node
    this%volume = volume
    this%face_area = face_area
    this%face_normal = face_normal

  end subroutine init

end module cell_geom_type
