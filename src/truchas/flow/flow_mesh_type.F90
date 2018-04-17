#include "f90_assert.fpp"

module flow_mesh_type

  use kinds, only: r8
  use unstr_mesh_type
  implicit none
  private

  ! fcell(1:2,i) -> contains the 2 cellids of the cells which share this face
  !                 In general, only the entries for i in [1,nface_onP] can be relied on
  !                 The first entry is always the cell for which the face normal points inward
  !                 The second entry is always the cell for which the face normal points outward
  ! for i in [1,nface_onP], the first entry will only be 0 at the boundary
  !                         the second entry will never be 0.

  type, public :: flow_mesh
    type(unstr_mesh), pointer :: mesh => null() ! read only
    integer, allocatable :: fcell(:,:) ! face cells
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
    integer :: i, j, s

    this%mesh => m
    allocate(this%face_centroid(3, m%nface))
    allocate(this%cell_centroid(3, m%ncell))

    do i = 1, m%nface
      associate(fn => m%fnode(m%xfnode(i):m%xfnode(i+1)-1))
        this%face_centroid(:,i) = sum(m%x(:,fn))/real(size(fn),r8)
      end associate
    end do

    ! TODO: the center-of-mass is not necessarily the arithmetic mean of the
    ! vertices.  This is a simple approximation which may or may not matter
    do i = 1, m%ncell
      associate(cn => m%cnode(m%xcnode(i):m%xcnode(i+1)-1))
        this%cell_centroid(:,i) = sum(m%x(:,cn))/real(size(cn),r8)
      end associate
    end do

    allocate(this%fcell(2,m%nface))
    this%fcell = 0

    do i = 1, m%ncell
      associate (fn => m%cface(m%xcface(i):m%xcface(i+1)-1))
        do j = 1, size(fn)
          if (fn(j) == 0) cycle
          if (btest(m%cfpar(i,pos=j))) then
            this%fcell(1,fn(j)) = i
          else
            this%fcell(2,fn(j)) = i
          endif
        end do
      end associate
    end do

    !! DEBUGGING CHECK
    do i = 1, m%ncell_onP
      associate (cn => m%cnhbr(m%xcnhbr(i):m%xcnhbr(i+1)-1), &
          fi => m%cface(m%xcface(i):m%xcface(i+1)-1))
        do j = 1, size(cn)
          if (cn(j) > 0) then
            assert(all(this%fcell(:,fi(j))) > 0)
            if (btest(m%cfpar(i,pos=j))) then
              assert(this%fcell(1,fi(j)) == i)
              assert(this%fcell(2,fi(j)) == cn(j))
            else
              assert(this%fcell(1,fi(j)) == cn(j))
              assert(this%fcell(2,fi(j)) == i)
            end if
          else
            assert(this%fcell(1,fi(j)) == 0)
            assert(this%fcell(2,fi(j)) == i)
          endif
        end do
      end associate
    end do
  end subroutine init
end module flow_mesh_type
