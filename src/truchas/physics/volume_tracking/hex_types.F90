!!
!! hex_types
!!
!! note: this could be derived from the polyhedron type
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

#include "f90_assert.fpp"

module hex_types

  use kinds, only: r8
  use plane_type
  use cell_topology
  implicit none
  private

  ! hex type to make divide and conquer algorithm simpler
  type, public :: base_hex
    real(r8) :: node(3,8), volume
  contains
    procedure :: calc_volume
    ! procedure         :: cell_center
    ! procedure         :: face_centers

    ! procedure         :: edge_centers
    ! procedure         :: contains_interface
  end type base_hex

  type, extends(base_hex), public :: cell_data
    real(r8) :: face_area(6), face_normal(3,6)
  contains
    procedure :: init => init_cell_data
    procedure :: calc_face_areas_and_normals
  end type cell_data

  type, extends(base_hex), public :: reconstruction_hex
    type(plane) :: P
    real(r8)    :: int_area  ! area of the interface for materials
    real(r8)    :: vof
    ! contains
    !   procedure              :: locate_plane
  end type reconstruction_hex

  ! cube vertex positions for unit testing
  real(r8), parameter, public :: cube_v(3,8) = reshape([ &
       0.0_r8, 0.0_r8, 0.0_r8, &
       1.0_r8, 0.0_r8, 0.0_r8, &
       1.0_r8, 1.0_r8, 0.0_r8, &
       0.0_r8, 1.0_r8, 0.0_r8, &
       0.0_r8, 0.0_r8, 1.0_r8, &
       1.0_r8, 0.0_r8, 1.0_r8, &
       1.0_r8, 1.0_r8, 1.0_r8, &
       0.0_r8, 1.0_r8, 1.0_r8],&
       shape(cube_v))

contains

  subroutine init_cell_data (this, node, volume, face_area, face_normal)

    class(cell_data), intent(out) :: this
    real(r8), intent(in)  :: node(:,:)
    real(r8), intent(in)  :: volume, face_area(:), face_normal(:,:)

    this%node = node
    this%volume = volume
    this%face_area = face_area
    this%face_normal = face_normal

  end subroutine init_cell_data

  ! calculates the volume of a hex
  real(r8) function calc_volume (this)
    use cell_geometry, only: eval_hex_volumes
    class(base_hex), intent(in) :: this
    real(r8) :: cvol_tmp(8)
    call eval_hex_volumes(this%node, calc_volume, cvol_tmp)
  end function calc_volume

  subroutine calc_face_areas_and_normals (this)
    use cell_geometry, only: quad_face_normal, vector_length

    class(cell_data), intent(inout) :: this

    integer :: f

    do f = 1,6
      this%face_normal(:,f) = quad_face_normal(this%node(:,HEX8_FACES(HEX8_XFACE(f):HEX8_XFACE(f+1)-1)))
      this%face_area(f) = vector_length(this%face_normal(:,f))
      this%face_normal(:,f) = this%face_normal(:,f) / sqrt(sum(this%face_normal(:,f)**2))

      ! ensure the normals are outward facing
      this%face_normal(:,f) = calculate_outward_normal (this%face_normal(:,f), sum(this%node, dim=2)/8.0_r8, &
           sum(this%node(:,HEX8_FACES(HEX8_XFACE(f):HEX8_XFACE(f+1)-1)), dim=2)/4.0_r8)
    end do

  end subroutine calc_face_areas_and_normals

  function calculate_outward_normal (normal, cell_center, face_center) result(outward_normal)
    real(r8), intent(in) :: normal(:), cell_center(:), face_center(:)
    real(r8)             :: outward_normal(3)

    real(r8) :: outward_dir(3)

    outward_dir = face_center - cell_center

    if ( sum(normal*outward_dir)>0.0_r8 ) then
      outward_normal = normal
    else
      outward_normal = -normal
    end if

  end function calculate_outward_normal

end module hex_types
