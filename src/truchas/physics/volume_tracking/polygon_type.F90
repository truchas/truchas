!!
!! polygon_type
!!
!! This module defines an arbitrary polygon type, along with routines for
!! calculating area, splitting, locating intersections, etc.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! References:
!!     1. Mirtich. Fast and Accurate Computation of Polehedral Mass Properties.
!!        Journal of Graphics Tools, 1996.
!!

#include "f90_assert.fpp"

module polygon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  implicit none
  private

  type, public :: polygon
    integer :: nVerts
    real(r8), allocatable :: x(:,:) !, norm(:)
    real(r8) :: norm(3)
  contains
    procedure :: init => init_polygon
    procedure :: centroid
    procedure :: order
    procedure :: basis
    procedure :: update_plane_normal
    procedure :: print_data
  end type polygon

contains

  subroutine init_polygon (this, x, norm)
    class(polygon),     intent(out)   :: this
    real(r8),           intent(in)    :: x(:,:)
    real(r8), optional, intent(inout) :: norm(:)

    this%nVerts = size(x, dim=2)
    this%x = x

    call this%update_plane_normal (norm)

  end subroutine init_polygon

  ! calculate the normal via cross product from vectors defined by 3 vertices
  subroutine update_plane_normal (this,norm)

    use cell_geometry, only: cross_product, normalized
    use near_zero_function

    class(polygon),     intent(inout) :: this
    real(r8), optional, intent(inout) :: norm(:) ! return the newly calculated norm if it wasn't known

    integer :: i,j

    ! the direction of the normal is assumed from the node ordering (and assuming convex)

    ! if the polygon has >3 verteces, it could be non-planar and we should subdivide
    if (present(norm)) then
      ASSERT(size(norm)==3)
      this%norm = norm
    else
      !if (allocated(this%norm)) deallocate(this%norm)
      !if (.not.allocated(this%norm)) allocate(this%norm(3))
      this%norm = 0.0_r8
    end if

    if (all(near_zero (this%norm))) then
      i = 3
      ! make sure we pick 3 vertices that don't all lie in the same plane
      do while (all(near_zero (this%norm)) .and. i<=this%nVerts)
        this%norm = normalized(cross_product (this%x(:,2) - this%x(:,1), this%x(:,i) - this%x(:,1)))
        i = i + 1
      end do
      if (i>this%nVerts .and. all(near_zero(this%norm))) then
        call this%print_data ()
        call TLS_fatal ("polygon only consists of a line")
      end if
      if (present(norm)) norm = this%norm
    end if

  end subroutine update_plane_normal

  function centroid (this)
    class(polygon), intent(in) :: this
    real(r8)                   :: centroid(3)
    centroid = sum(this%x, dim=2) / this%nVerts
  end function centroid

  ! order the vertices of a convex polygon
  !
  ! this is done by calculating the vector between each vertex and the polygon centroid,
  ! then the angle of that vector with respect to the x-axis in that space.
  ! this angle is used to sort the vertices

  subroutine order (this,array)

    use cell_geometry, only: normalized
    use sort_utilities, only: insertion_sort
    use permutations, only: invert_perm

    class(polygon),    intent(inout) :: this
    integer, optional, intent(inout) :: array(:) ! an array that gets sorted along with the polygon

    real(r8), allocatable :: q(:,:)
    real(r8) :: xc(3), t(this%nVerts), t2(this%nVerts), prjx(2), xl(3,this%nVerts), tmp
    integer :: i,ind(size(this%x,dim=2))

    ! calculate the location of the centroid, and the vertex coordinates with respect to the centroid
    ! normalize to avoid floating point cutoffs, since the polygon may be very tiny
    xc = this%centroid ()
    do i = 1,this%nVerts
      xl(:,i) = normalized(this%x(:,i) - xc(:))
    end do

    ! the projection coordinate directions
    q = orthonormalBasis(xl)

    ! calculate the rotation angle
    do i = 1,this%nVerts
      ! get coordinates for the vertex in the 2D plane defined by the polygon
      prjx(1) = dot_product(xl(:,i),q(:,1))
      prjx(2) = dot_product(xl(:,i),q(:,2))

      ! find the angle made by this vertex with respect to the first vertex
      t(i) = atan2(prjx(2),prjx(1))
    end do

    ! sort based on angle
    if (present(array)) then
      t2 = t
      ind = [(i, i=1,size(this%x,dim=2))]
      call insertion_sort (ind,t2)
      call invert_perm (ind)
      do i = 1,size(array)
        if (array(i)>0) array(i) = ind(array(i))
      end do
    end if

    call insertion_sort (this%x,t)

    call this%update_plane_normal ()

  end subroutine order

  ! generate an orthogonal basis for the polygon, approximately scaled to the size of the polygon
  ! the first element is the shortest dimension, the second the longest
  function basis (this)

    class(polygon), intent(in) :: this
    real(r8) :: basis(3,2)

    integer :: i, iN
    real(r8) :: xc(3), xl(3), minlen, maxlen, length

    ! find the minimum and maximum directions
    maxlen = 0.0_r8
    minlen = huge(1.0_r8)
    xc = this%centroid ()
    ! do i = 1,this%nVerts
    !   xl = xc - this%x(:,i)
    !   length = norm2(xl)

    !   if (length < minlen) then
    !     minlen = length
    !     basis(:,1) = xl
    !   else if (length > maxlen) then
    !     maxlen = length
    !     basis(:,2) = xl
    !   end if
    ! end do

    ! find closest and farthest points to the centroid on the polygon edge
    do i = 1,this%nVerts
      iN = modulo(i,this%nVerts) + 1
      xl = this%x(:,i) + projectOnto(xc - this%x(:,i), this%x(:,iN) - this%x(:,i)) - xc
      length = norm2(xl)

      if (length < minlen) then
        minlen = length
        basis(:,1) = xl
      else if (length > maxlen) then
        maxlen = length
        basis(:,2) = xl
      end if
    end do

    ! align basis(:,2) such that the two vectors are orthogonal
    ! it should be almost orthogonal as is
    basis(:,2) = basis(:,2) - projectOnto(basis(:,2),basis(:,1))

    call this%print_data ()
    print '(a,3es15.5)', 'b: ',basis(:,1)
    print '(a,3es15.5)', 'b: ',basis(:,2)
    print *

  end function basis

  subroutine print_data (this)

    class(polygon), intent(in) :: this

    integer :: v

    print *, 'POLYGON DATA:'
    if (allocated(this%x)) then
      do v = 1,this%nVerts
        print '(a,i3,a,3es15.5)', 'x ',v,':  ',this%x(:,v)
      end do
      write(*,*)
    end if

    print '(a,3es15.5)', 'norm ',this%norm

  end subroutine print_data

  ! given a set of vectors x, return an orthonormal basis q for the same space
  function orthonormalBasis (x)

    use truchas_logging_services

    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable :: orthonormalBasis(:,:)

    integer :: m, n, k, lwork, ierr
    real(r8), allocatable :: Q(:,:), tau(:), work(:)

    ! copy the input variables to a matrix to be modified in-place
    m = size(x, dim=1)
    n = size(x, dim=2)
    k = min(m,n)
    lwork = 8*n
    allocate(work(lwork), tau(k))
    Q = x

    ! compute the QR factorization
    call dgeqrf (m, n, Q, m, tau, work, lwork, ierr)
    if (ierr/=0) call TLS_fatal ("failed dgeqrf in orthonormalBasis")

    call dorgqr (m, m, k, Q, m, tau, work, lwork, ierr)
    if (ierr/=0) call TLS_fatal ("failed dorgqr in orthonormalBasis")

    ! return the first m columns of Q as the orthonormal basis with the same span as x
    orthonormalBasis = Q(:,1:m)

  end function orthonormalBasis

  ! project vector x1 into the direction of vector x2
  pure function projectOnto (x1, x2)
    use cell_geometry, only: normalized
    real(r8), intent(in) :: x1(:), x2(:)
    real(r8)             :: projectOnto(size(x1))
    real(r8)             :: n2(size(x1))
    n2 = normalized(x2)
    projectOnto = dot_product(x1, n2) * n2
  end function projectOnto

end module polygon_type
