!!
!! provides a derived type for packing together cell geometry data
!!

#include "f90_assert.fpp"

module cell_geom_type

  use kinds, only: r8
  implicit none
  private

  integer, parameter, public :: &
      CELL_TET4 = 1, &
      CELL_PYR5 = 2, &
      CELL_WED6 = 3, &
      CELL_HEX8 = 4

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
    case (4) ! tet
      this%nfc = 4
      this%cell_type = CELL_TET4
    case (5) ! pyramid
      this%nfc = 5
      this%cell_type = CELL_PYR5
    case (6) ! wedge
      this%nfc = 5
      this%cell_type = CELL_WED6
    case (8) ! hex
      this%nfc = 6
      this%cell_type = CELL_HEX8
    case default
      call TLS_fatal('unaccounted topology in truncation_volume_type')
    end select

  end subroutine init

end module cell_geom_type
