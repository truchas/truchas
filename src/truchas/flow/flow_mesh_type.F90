#include "f90_assert.fpp"

module flow_mesh_type

  use kinds, only: r8
  use unstr_mesh_type
  implicit none
  private

  type, public :: flow_mesh
    type(unstr_mesh), pointer :: mesh => null() ! read only
    integer, allocatable :: xfcell(:), fcell(:) ! face cells
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
    integer, allocatable :: fp(:), fm(:)

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


    ! build face->cell indexing
    allocate(fp(m%nface))
    allocate(fm(m%nface))
    fp = 0
    fm = 0

    do i = 1, m%ncell
      associate (fn => m%cface(m%xcface(i):m%xcface(i+1)-1))
        do j = 1, size(fn)
          if (btest(m%cfpar(i,pos=j))) then
            fp(fn(j)) = i
          else
            fm(fn(j)) = i
          endif
        end do
      end associate
    end do
    s = count(fp > 0) + count(fm > 0)
    allocate(this%fcell(s))
    allocate(this%xfcell(m%nface+1))

    j = 1
    this%xfcell(i) = j
    do i = 1, m%nface
      if (fp(i) > 0 .and. fm(i) > 0) then
        this%fcell(j) = fp(i)
        this%fcell(j+1) = fm(i)
        j = j+2
      elseif (fp(i) > 0) then
        this%fcell(j) = fp(i)
        j = j+1
      elseif (fm(i) > 0) then
        this%fcell(j) = fm(i)
        j = j+1
      end if
      this%xfcell(i+1) = j
    end do

    !! DEBUGGING CHECK
    do i = 1, m%ncell_onP
      associate (cn => m%cnhbr(m%xcnhbr(i):m%xcnhbr(i+1)-1), &
          fi => m%cface(m%xcface(i):m%xcface(i+1)-1))

        do j = 1, size(cn)
          associate(fc => this%fcell(this%xfcell(fi(j)):this%xfcell(fi(j)+1)-1))
            if (cn(j) > 0) then
              assert(size(fc) == 2)
              if (fc(1) == i) then
                assert(fc(2) == cn(j))
              else
                assert(fc(1) == cn(j))
                assert(fc(2) == i)
              end if
            else
              assert(size(fc) == 1)
              assert(fc(1) == i)
            endif
          end do
        end do
      end associate
    end do

    deallocate(fp)
    deallocate(fm)
  end subroutine init
end module flow_mesh_type
