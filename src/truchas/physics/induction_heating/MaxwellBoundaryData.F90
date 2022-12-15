!!
!! The MaxwellBoundaryData Module
!!
!! Neil N Carlson <nnc@newmexico.com> 10 Jul 2003
!! Last revised 19 Jul 2003
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module MaxwellBoundaryData

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  implicit none
  private

  type, public :: efield_bndry_func
    type(simpl_mesh), pointer :: mesh => null()
    integer :: nebgroup = 0
    integer, allocatable :: ebedge(:)
    real(r8), allocatable :: eb(:)
  contains
    procedure :: init
    procedure :: set_Eb_function, set_Eb_values
  end type

  type, public :: hfield_bndry_func
    type(simpl_mesh), pointer :: mesh => null()
    integer :: nhbgroup = 0
    integer, allocatable :: hbface(:)
    real(r8), allocatable :: hb(:,:)
  contains
    procedure :: init => init_hfield_bndry_func
    procedure :: set_Hb_function, get_Hb_source
  end type

contains

  subroutine init(this, mesh, setids)

    use bitfield_type

    class(efield_bndry_func), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    integer, intent(in) :: setids(:)

    integer :: i, j, k, stat
    type(bitfield) :: bitmask

    this%mesh => mesh ! we'll need this for subsequent initialization steps

    allocate(this%ebedge(mesh%nedge), this%eb(mesh%nedge))
    this%ebedge = 0
    this%eb = 0.0_r8
    this%nebgroup = size(setids)
    do i = 1, this%nebgroup
      call this%mesh%get_face_set_bitmask([setids(i)], bitmask, stat)
      ASSERT(stat == 0)
      do j = 1, mesh%nface
        if (popcnt(iand(bitmask, this%mesh%face_set_mask(j))) /= 0) then
        !if (btest(this%mesh%face_set_mask(j), pos=ebgroup(i))) then
          do k = 1, 3
            if (this%ebedge(abs(mesh%fedge(k,j))) == 0) this%ebedge(abs(mesh%fedge(k,j))) = i
          end do
        end if
      end do
    end do

    !! Ensure consistency of values for edges on the partition boundary.
    call mesh%edge_imap%gather_offp(this%ebedge)

  end subroutine init

  subroutine set_Eb_function(this, n, f, order)

    class(efield_bndry_func), intent(inout) :: this
    integer, intent(in) :: n
    integer, intent(in), optional :: order

    interface
      function f (x) result (fx)
        use,intrinsic :: iso_fortran_env, only: r8 => real64
        real(r8), intent(in) :: x(:)
        real(r8) :: fx(3)
      end function f
    end interface

    integer :: j, ord
    real(r8), pointer :: x1(:), x2(:)

    ASSERT(allocated(this%ebedge))
    ! Assertion fails for 0-sized ebedge array
    !ASSERT(n > 0 .and. n <= maxval(this%ebedge))

    ord = 1 ! use trapezoid rule quadrature by default
    if (present(order)) ord = order

    do j = 1, this%mesh%nedge
      if (this%ebedge(j) /= n) cycle
      x1 => this%mesh%x(:,this%mesh%enode(1,j))
      x2 => this%mesh%x(:,this%mesh%enode(2,j))
      this%eb(j) = dot_product(x2 - x1, vector_average(x1, x2, f, ord))
    end do

  end subroutine set_Eb_function

  subroutine set_Eb_values(this, coef, e)

    class(efield_bndry_func), intent(in) :: this
    real(r8), intent(in) :: coef(:)
    real(r8), intent(inout) :: e(:)

    integer :: j

    ASSERT(size(coef) == this%nebgroup)
    ASSERT(size(e) == size(this%eb) )

    do j = 1, size(e)
      if (this%ebedge(j) == 0) cycle
      e(j) = coef(this%ebedge(j)) * this%eb(j)
    end do

  end subroutine set_Eb_values

  function vector_average (x1, x2, f, order) result (avg)
    real(r8), intent(in) :: x1(:), x2(:)
    integer, intent(in) :: order
    real(r8) :: avg(size(x1))
    interface
      function f (x) result (fx)
        use,intrinsic :: iso_fortran_env, only: r8 => real64
        real(r8), intent(in) :: x(:)
        real(r8) :: fx(3)
      end function f
    end interface
    select case (order)
    case (:0) ! midpoint rule
      avg = f(0.5_r8*(x1+x2))
    case (1)  ! trapezoid rule
      avg = 0.5_r8*(f(x1)+f(x2))
    case (2:) ! Simpson's rule
      avg = (f(x1)+4.0_r8*f(0.5_r8*(x1+x2))+f(x2))/6.0_r8
    end select
  end function vector_average

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_hfield_bndry_func(this, mesh, setids)

    use bitfield_type

    class(hfield_bndry_func), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    integer, intent(in) :: setids(:)

    integer :: i, j, stat
    type(bitfield) :: bitmask

    this%mesh => mesh

    allocate(this%hbface(mesh%nface), this%hb(mesh%nedge,size(setids)))
    this%hbface = 0
    this%nhbgroup = size(setids)
    do i = 1, this%nhbgroup
      call this%mesh%get_face_set_bitmask([setids(i)], bitmask, stat)
      ASSERT(stat == 0)
      do j = 1, mesh%nface
        if (popcnt(iand(bitmask, this%mesh%face_set_mask(j))) /= 0) then
          if (this%hbface(j) == 0) this%hbface(j) = merge(i, -i, this%mesh%fcell(2,j) == 0)
        end if
      end do
    end do

  end subroutine init_hfield_bndry_func

  subroutine set_Hb_function(this, n, f)

    class(hfield_bndry_func), intent(inout) :: this
    integer, intent(in) :: n

    integer :: j, edge(3)
    real(r8) :: h(3), c

    interface
      function f(x) result(fx)
        import r8
        real(r8), intent(in) :: x(:)
        real(r8) :: fx(3)
      end function
    end interface

    ASSERT(n > 0 .and. n <= this%nhbgroup)

    this%hb(:,n) = 0.0_r8
    do j = 1, this%mesh%nface
      if (abs(this%hbface(j)) /= n) cycle
      edge = this%mesh%fedge(:,j)
      c = sign(1, this%hbface(j)) / 6.0_r8
      call project_on_face_edges(this%mesh%x(:,this%mesh%fnode(:,j)), f, h)
      this%hb(edge(1),n) = this%hb(edge(1),n) + c * (h(3) - h(2))
      this%hb(edge(2),n) = this%hb(edge(2),n) - c * (h(1) - h(3))
      this%hb(edge(3),n) = this%hb(edge(3),n) + c * (h(2) - h(1))
    end do

  end subroutine set_Hb_function

  subroutine project_on_face_edges(x, f, p)

    real(r8), intent(in) :: x(:,:)
    real(r8), intent(out):: p(:)

    interface
      function f (x) result (fx)
        use,intrinsic :: iso_fortran_env, only: r8 => real64
        real(r8), intent(in) :: x(:)
        real(r8) :: fx(3)
      end function f
    end interface

    integer :: k
    real(r8) :: fv(3,3)

    ASSERT( size(x,dim=1) == 3 )
    ASSERT( size(x,dim=2) == 3 )
    ASSERT( size(p) == 3 )

    !! Length-weighted edge projections using trapezoid rule
    do k = 1, 3
      fv(:,k) = f(x(:,k))
    end do
    p(1) = 0.5_r8 * dot_product(x(:,3)-x(:,2), fv(:,3)+fv(:,2))
    p(2) = 0.5_r8 * dot_product(x(:,1)-x(:,3), fv(:,1)+fv(:,3))
    p(3) = 0.5_r8 * dot_product(x(:,2)-x(:,1), fv(:,2)+fv(:,1))

  end subroutine project_on_face_edges


  subroutine get_Hb_source(this, coef, bsrc)

    class(hfield_bndry_func), intent(in) :: this
    real(r8), intent(in) :: coef(:)
    real(r8), intent(out) :: bsrc(:)

    integer :: j, k

    ASSERT(size(coef) == size(this%hb,dim=2))
    ASSERT(size(bsrc) == size(this%hb,dim=1))

    bsrc = 0.0_r8
    do k = 1, size(coef)
      do j = 1, size(bsrc)
        bsrc(j) = bsrc(j) + coef(k) * this%hb(j,k)
      end do
    end do

  end subroutine get_Hb_source

end module MaxwellBoundaryData
