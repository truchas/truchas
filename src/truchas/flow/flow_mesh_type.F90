#include "f90_assert.fpp"

module flow_mesh_type

  use kinds, only: r8
  use unstr_mesh_type
  use index_partitioning
  use parallel_communication
  use bitfield_type
  implicit none
  private

  type, public :: flow_mesh
    type(unstr_mesh), pointer :: mesh => null() ! read only
    real(r8), allocatable :: face_centroid(:,:)
    real(r8), allocatable :: cell_centroid(:,:)
  contains
    procedure :: init
  end type flow_mesh

contains

  subroutine init(this, m)
    class(flow_mesh), intent(inout) :: this
    type(unstr_mesh), pointer, intent(in) :: m
    !-
    integer :: i

    this%mesh = m
    allocate(this%face_centroid(3, m%nface))
    allocate(this%cell_centroid(3, m%ncell))

    do i = 1, m%nface
      associate(fn => m%fnode(m%xfnode(j):m%xfnode(j+1)-1))
        this%face_centroid(:,i) = sum(m%x(:,fn))/real(size(fn),r8)
      end associate
    end do

    ! TODO: the center-of-mass is not necessarily the arithmetic mean of the
    ! vertices.  This is a simple approximation which may or may not matter
    do i = 1, m%ncell
      associate(cn => m%cnode(m%xcnode(j):m%xcnode(j+1)-1))
        this%cell_centroid(:,i) = sum(m%x(:,cn))/real(size(cn),r8)
      end associate
    end do
  end subroutine init
end module flow_mesh_type
