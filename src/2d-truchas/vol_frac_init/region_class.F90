module region_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  private

  type, abstract, public :: region(dim)
    integer,len :: dim
  contains
    procedure(encloses), deferred :: encloses
  end type

  abstract interface
    logical function encloses(this, x)
      import region, r8
      class(region(*)), intent(in) :: this
      real(r8), intent(in) :: x(this%dim)
    end function
  end interface
  
  type, extends(region), public :: sphere
    real(r8) :: center(dim), radsq
  contains
    procedure :: encloses => sphere_encloses
  end type

!  type, extends(region2d) :: cell_set_region
!    private
!    integer :: region_mask
!    logical :: fill_outside = .false.
!  contains
!    procedure :: ifunc
!  end type
!
!  type, public :: region_indicator(num_reg)
!    private
!    integer,len :: num_reg
!    type(region_box) :: reg_array(num_reg)
!  contains
!    procedure :: region_at_cell_point
!  end type

contains

!  integer function 
!
!  logical function contains_2d(this, x) result(contains)
!    class(region2d), intent(in) :: this
!    real(r8), intent(in) :: x(:)
!    if (size(x) == 2) then
!      contains = this%contains_(x)
!    else
!      error stop 'foo'
!    end if
!  end function

  logical function sphere_encloses(this, x)
    class(sphere(*)), intent(in) :: this
    real(r8), intent(in) :: x(this%dim)
    sphere_encloses = (sum((x-this%center)**2) <= this%radsq)
  end function

!end module region_class
!
!module foo
!
!  type :: tri_cell
!    real(r8) :: x(2,3)
!    integer  :: vertex_region(3)
!    type(foo), pointer :: reg_id
!  contains
!    procedure :: subdivide
!  end type
!    
!contains
!
!
!
!  !! Subdivide the triangular cell THIS into 4 congruent triangular cells
!  !! and return the result in the array SUBTRI.
!
!  subroutine subdivide(this, subtri)
!    class(tri_cell), intent(in) :: this
!    type(tri_cell), intent(out) :: subtri(4)
!    real(r8) :: xmid1(2), xmid2(2), xmid(3) ! edge midpoints
!    xmid1 = 0.5_r8*(this%x(:,2) + this%x(:,3))
!    xmid2 = 0.5_r8*(this%x(:,1) + this%x(:,3))
!    xmid3 = 0.5_r8*(this%x(:,1) + this%x(:,2))
!    subtri(1)%x(:,1) = this%x(:,1)
!    subtri(1)%x(:,2) = xmid3
!    subtri(1)%x(:,3) = xmid2
!    subtri(2)%x(:,1) = xmid3
!    subtri(2)%x(:,2) = this%x(:,2)
!    subtri(2)%x(:,3) = xmid1
!    subtri(3)%x(:,1) = xmid2
!    subtri(3)%x(:,2) = xmid1
!    subtri(3)%x(:,3) = this%x(:,3)
!    subtri(4)%x(:,1) = xmid1
!    subtri(4)%x(:,2) = xmid2
!    subtri(4)%x(:,3) = xmid3
!    
!    regid1 = this%reg_id%eval(xmid1)
!    regid2 = this%reg_id%eval(xmid2)
!    regid3 = this%reg_id%eval(xmid3)
!    subtri(1)%reg_id(1) = this%reg_id(1)
!    subtri(1)%reg_id(2) = regid3
!    subtri(1)%reg_id(3) = regid2
!    subtri(2)%reg_id(1) = regid3
!    subtri(2)%reg_id(2) = this%reg_id(2)
!    subtri(2)%reg_id(3) = regid1
!    subtri(2)%reg_id(1) = regid3
!    subtri(2)%reg_id(2) = this%reg_id(2)
!    subtri(2)%reg_id(3) = regid1
!    subtri(4)%reg_id(1) = regid1
!    subtri(4)%reg_id(2) = regid2
!    subtri(4)%reg_id(3) = regid3
!  end subroutine
!
!  !! Subdivide an oriented convex polygon with nodes x(:,j), j = 1, n,
!  !! into an array of triangular cells.
!
!  subroutine triangulate(x, subtri)
!    real(r8), intent(in) :: x(:,:)
!    type(tri_cell), allocatable, intent(out) :: subtri(:)
!    integer :: n, i
!    n = size(x,dim=2)
!    ASSERT(n >= 3)
!    allocate(subtri(n-2))
!    call triangulate_aux(x, [(i, i=1,n)], subtri)
!  end subroutine
!
!  !! Recursively subdivide a convex polygon into an array of triangular cells.
!  !! The nodes of the oriented polygon are x(:,cnode(i)), i = 1, size(cnode).
!
!  recursive subroutine triangulate_aux(x, cnode, subtri)
!
!    real(r8), intent(in) :: x(:,:)
!    integer, intent(in) :: cnode(:)
!    type(tri_cell), intent(inout) :: subtri(:)
!
!    integer :: i, j, n, m, p, q
!    integer, allocatable :: subcnode(:)
!    real(r8) :: tmp, minlen
!
!    n = size(cnode)
!    ASSERT(n >= 3)
!    ASSERT(size(subtri) == n-2)
!
!    if (n == 3) then ! it's a triangle; nothing much to do
!      subtri(1)%x = x(:,cnode)
!      ! other stuff?
!    else ! split polygon along shortest diagonal and recurse
!      minlen = huge(minlen)
!      do i = 1, n-2
!        do j = i+2, merge(n-1, n, i==1)
!          tmp = norm2(x(:,cnode(i))-x(:,cnode(j)))
!          if (tmp < minlen) then
!            minlen = tmp
!            p = i
!            q = j
!          end if
!        end do
!      end do
!      subcnode = [(cnode(i), i = p, q)]
!      m = size(subcnode) - 2
!      call triangulate_aux(x, subcnode, subtri(:m))
!      subcnode = [(cnode(i), i = 1, p), (cnode(i), i = q, n)]
!      call triangulate_aux(x, subcnode, subtri(m+1:))
!    end if
!
!  end subroutine triangulate_aux
!
!  recursive function volumes(this, rlev)
!    class(tri_cell), intent(in) :: this
!    integer, intent(in) :: rlev
!    real(r8) :: volumes(???)
!    
!    volumes = 0.0_r8
!    n = this%reg_id(1)
!    if (any(this%reg_id(2:3) /= n)) then ! multi-region cell
!      if (rlev > 0) then ! recurse
!        block
!          type(tri_cell) :: subtri(4)
!          call subdivide(this, subtri)
!          do i = 1, 4
!            volumes = volumes + subtri(i)%volumes(rlev-1)
!          end do
!          volumes = 0.25_r8 * volumes
!        end block
!      else
!        n = region-at-centroid
!        volumes(n) = 1.0_r8
!      end if
!    else
!      n = this%reg_id(1)
!      volumes(n) = 1.0_r8
!    end if
!
!  end function
!  
!  subroutine get_vol_frac(x, ???)
!  
!    real(r8), intent(in) :: x(:,:)
!    type(region_indicator), intent(in) :: reg_id
!    
!    
!  end subroutine

end module
