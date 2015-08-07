!!
!! FIELD_PROBES
!!
!! This module provides some tools for probing the discrete scalar and vector
!! fields in mimetic discretization methods on tetrahedral meshes at specific
!! points, and recording the values over time.
!!
!! This module formed part of my existing EM solver that was imported into
!! Truchas.
!!
!! Neil N. Carlson <nnc@newmexico.com>
!! Last revised 15 Mar 2004.
!!
!! IMPLEMENTATION NOTES
!!
!! This module has grown over time, and the higher level structures really need
!! to be redesigned before any further significant development work is done.
!!

#include "f90_assert.fpp"

module field_probes

  use simplex_geometry, only: tet_volume, tet_face_normal
  use simpl_mesh_type
  use mimetic_discretization
  implicit none
  private
  
  public  :: new_probe_array, delete_probe_array, set_probe_point, initialize_probe_array, &
             update_probe_array, num_w1_probes, num_w2_probes, num_w3_probes, get_probe_times, &
             get_w1_probe_values, get_w2_probe_values, get_w3_probe_values
  private :: initialize_w1_probe, initialize_w2_probe, initialize_w3_probe, &
             extend_w1_probe, extend_w2_probe, extend_w3_probe, &
             w1_probe_value, w2_probe_value, w3_probe_value, cell, outside
  
  integer, parameter, private :: r8 = selected_real_kind(10,50)
  integer, parameter, private :: DIMEN = 3  ! spatial dimension
  integer, parameter, private :: NVERT = 4  ! number of vertices
  integer, parameter, private :: NEDGE = 6  ! number of edges
  integer, parameter, private :: NFACE = 4  ! number of faces
  
  type, public :: ProbeArray
    private
    type(W1Probe), dimension(:), pointer :: w1_probe => null()
    type(W2Probe), dimension(:), pointer :: w2_probe => null()
    type(W3Probe), dimension(:), pointer :: w3_probe => null()
    real, dimension(:), pointer :: time => null()
    integer :: n = 0
    integer :: max = 0
    integer :: n_init = 100 ! initial size
    integer :: n_incr = 100
  end type ProbeArray
  
  type, private :: W0Probe
    logical       :: defined = .false.
    real(kind=r8) :: point(DIMEN)
    real(kind=r8) :: mult(NVERT)
    integer       :: index(NVERT)
    real, pointer :: value(:) => null()
  end type W0Probe
  
  type, private :: W1Probe
    logical       :: defined = .false.
    real(kind=r8) :: point(DIMEN)       ! probe point
    real(kind=r8) :: mult(DIMEN,NEDGE)  ! interpolation multipliers
    integer       :: index(NEDGE)       ! interpolation edges
    real, pointer :: value(:,:) => null()
  end type W1Probe
  
  type, private :: W2Probe
    logical       :: defined = .false.
    real(kind=r8) :: point(DIMEN)       ! probe point
    real(kind=r8) :: mult(DIMEN,NFACE)  ! interpolation multipliers
    integer       :: index(NFACE)       ! interpolation faces
    real, pointer :: value(:,:) => null()
  end type W2Probe
  
  type, private :: W3Probe
    logical       :: defined = .false.
    real(kind=r8) :: point(DIMEN)       ! probe point
    real(kind=r8) :: mult(1)            ! interpolation multiplier
    integer       :: index(1)           ! interpolation cell
    real, pointer :: value(:) => null()
  end type W3Probe
  
contains

  subroutine new_probe_array (this, n)
    type(ProbeArray), intent(inout) :: this
    integer, intent(in) :: n
    integer :: j
    allocate(this%w1_probe(n), this%w2_probe(n), this%w3_probe(n))
    this%max = this%n_init
    allocate(this%time(this%max))
    do j = 1, n
      call extend_w1_probe (this%w1_probe(j), this%max)
      call extend_w2_probe (this%w2_probe(j), this%max)
      call extend_w3_probe (this%w3_probe(j), this%max)
    end do
  end subroutine new_probe_array
  
  subroutine delete_probe_array (this)
    type(ProbeArray), intent(inout) :: this
    integer :: j
    if (associated(this%w1_probe)) then
      do j = 1, size(this%w1_probe)
        deallocate(this%w1_probe(j)%value)
      end do
      deallocate(this%w1_probe)
    end if
    if (associated(this%w2_probe)) then
      do j = 1, size(this%w2_probe)
        deallocate(this%w2_probe(j)%value)
      end do
      deallocate(this%w2_probe)
    end if
    if (associated(this%w3_probe)) then
      do j = 1, size(this%w3_probe)
        deallocate(this%w3_probe(j)%value)
      end do
      deallocate(this%w3_probe)
    end if
    if (associated(this%time)) deallocate(this%time)
    this%n = 0
    this%max = 0
  end subroutine delete_probe_array
    
  subroutine set_probe_point (this, n, point)
    type(ProbeArray), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=r8), dimension(:), intent(in) :: point
    ASSERT( size(point) == DIMEN )
    this%w1_probe(n)%point = point
    this%w2_probe(n)%point = point
    this%w3_probe(n)%point = point
  end subroutine set_probe_point
  
  subroutine initialize_probe_array (this, d, mask)
    type(ProbeArray),  intent(inout) :: this
    type(simpl_mesh),   intent(in)    :: d
    logical, optional, intent(out)   :: mask(:)
    integer :: j
    !! Each type of field probe is at the same point => one mask.
    do j = 1, size(this%w1_probe)
      call initialize_w1_probe (this%w1_probe(j), d)
    end do
    do j = 1, size(this%w2_probe)
      call initialize_w2_probe (this%w2_probe(j), d)
    end do
    do j = 1, size(this%w3_probe)
      call initialize_w3_probe (this%w3_probe(j), d)
    end do
    if (present(mask)) then
      ASSERT( size(mask) == this%n )
      mask = this%w1_probe%defined .and. this%w2_probe%defined .and. this%w3_probe%defined
    else if (.not.all(this%w1_probe%defined .and. &
                      this%w2_probe%defined .and. &
                      this%w3_probe%defined)) then
      print *, 'initialize_probe_array: PANIC! probe point belongs to no cell!'
      stop
    end if
  end subroutine initialize_probe_array

  subroutine initialize_w1_probe (probe, d)
    type(W1Probe), intent(inout) :: probe
    type(simpl_mesh), intent(in) :: d
    integer :: n
    n = cell(d, probe%point) ! Locate probe point in the grid
    if (n == 0) return
    probe%defined = .true.
    call eval_w1_interp_coef (d, probe%point, n, probe%mult, probe%index)
  end subroutine initialize_w1_probe
  
  subroutine initialize_w2_probe (probe, d)
    type(W2Probe), intent(inout) :: probe
    type(simpl_mesh), intent(in) :: d
    integer :: n
    n = cell(d, probe%point) ! Locate probe point in the grid
    if (n == 0) return
    probe%defined = .true.
    call eval_w2_interp_coef (d, probe%point, n, probe%mult, probe%index)
  end subroutine initialize_w2_probe
  
  subroutine initialize_w3_probe (probe, d)
    type(W3Probe), intent(inout) :: probe
    type(simpl_mesh), intent(in) :: d
    integer :: n
    n = cell(d, probe%point) ! Locate probe point in the grid
    if (n == 0) return
    probe%defined = .true.
    call eval_w3_interp_coef (d, probe%point, n, probe%mult, probe%index)
  end subroutine initialize_w3_probe
  
  subroutine update_probe_array (this, time, w1_field, w2_field, w3_field)
  
    type(ProbeArray), intent(inout) :: this
    real, intent(in) :: time
    real(kind=r8), intent(in) :: w1_field(:), w2_field(:), w3_field(:)
    
    integer :: j, n, old_size
    real, dimension(:), pointer :: old_time
    
    n = 1 + this%n
    
    !! Verify space exists; extend arrays if necessary
    if (n > this%max) then  ! out of space; allocate more
      if (associated(this%time)) then
        old_size = size(this%time)
        this%max = old_size + this%n_incr
        old_time => this%time
        allocate(this%time(this%max))
        this%time(1:old_size) = old_time
        deallocate(old_time)
      else  ! initial allocation
        this%max = this%n_init
        allocate(this%time(this%max))
      end if
      call extend_w1_probe (this%w1_probe, this%max)
      call extend_w2_probe (this%w2_probe, this%max)
      call extend_w3_probe (this%w3_probe, this%max)
    end if

    this%n = n
    this%time(n) = time
    
    do j = 1, size(this%w1_probe)
      if (this%w1_probe(j)%defined) then
        this%w1_probe(j)%value(:,n) = w1_probe_value(this%w1_probe(j), w1_field)
      end if
    end do
    
    do j = 1, size(this%w2_probe)
      if (this%w2_probe(j)%defined) then
        this%w2_probe(j)%value(:,n) = w2_probe_value(this%w2_probe(j), w2_field)
      end if
    end do
 
    do j = 1, size(this%w3_probe)
      if (this%w3_probe(j)%defined) then
        this%w3_probe(j)%value(n) = w3_probe_value(this%w3_probe(j), w3_field)
      end if
    end do
 
  end subroutine update_probe_array
  
  elemental subroutine extend_w1_probe (probe, max)
  
    type(W1Probe), intent(inout) :: probe
    integer, intent(in) :: max
    
    integer :: old_size
    real, pointer :: old_value(:,:)
    
    !! ASSERT( size(probe%value,dim=2) < max )
    
    if (associated(probe%value)) then
      old_size = size(probe%value,dim=2)
      old_value => probe%value
      allocate(probe%value(DIMEN,max))
      probe%value(:,1:old_size) = old_value
      deallocate(old_value)
    else
      allocate(probe%value(DIMEN,max))
    end if
    
  end subroutine extend_w1_probe
  
  elemental subroutine extend_w2_probe (probe, max)
  
    type(W2Probe), intent(inout) :: probe
    integer, intent(in) :: max
    
    integer :: old_size
    real, pointer :: old_value(:,:)
    
    !! ASSERT( size(probe%value,dim=2) < max )
    
    if (associated(probe%value)) then
      old_size = size(probe%value,dim=2)
      old_value => probe%value
      allocate(probe%value(DIMEN,max))
      probe%value(:,1:old_size) = old_value
      deallocate(old_value)
    else
      allocate(probe%value(DIMEN,max))
    end if
      
  end subroutine extend_w2_probe
    
  elemental subroutine extend_w3_probe (probe, max)
  
    type(W3Probe), intent(inout) :: probe
    integer, intent(in) :: max
    
    integer :: old_size
    real, pointer :: old_value(:)
    
    !! ASSERT( size(probe%value) < max )
    
    if (associated(probe%value)) then
      old_size = size(probe%value)
      old_value => probe%value
      allocate(probe%value(max))
      probe%value(1:old_size) = old_value
      deallocate(old_value)
    else
      allocate(probe%value(max))
    end if
      
  end subroutine extend_w3_probe
    
  
  function w1_probe_value (this, field) result (value)
    type(W1Probe), intent(in) :: this
    real(kind=r8), intent(in) :: field(:)
    real, dimension(DIMEN) :: value
    value = matmul(this%mult, field(this%index))
  end function w1_probe_value
  
  function w2_probe_value (this, field) result (value)
    type(W2Probe), intent(in) :: this
    real(kind=r8), intent(in) :: field(:)
    real :: value(DIMEN)
    value = matmul(this%mult, field(this%index))
  end function w2_probe_value
  
  function w3_probe_value (this, field) result (value)
    type(W3Probe), intent(in) :: this
    real(kind=r8), intent(in) :: field(:)
    real :: value
    value = this%mult(1) * field(this%index(1))
  end function w3_probe_value
  

  integer function num_w1_probes (this)
    type(ProbeArray), intent(in) :: this
    if (associated(this%w1_probe)) then
      num_w1_probes = size(this%w1_probe)
    else
      num_w1_probes = 0
    end if
  end function num_w1_probes
  
  integer function num_w2_probes (this)
    type(ProbeArray), intent(in) :: this
    if (associated(this%w2_probe)) then
      num_w2_probes = size(this%w2_probe)
    else
      num_w2_probes = 0
    end if
  end function num_w2_probes
  
  integer function num_w3_probes (this)
    type(ProbeArray), intent(in) :: this
    if (associated(this%w3_probe)) then
      num_w3_probes = size(this%w3_probe)
    else
      num_w3_probes = 0
    end if
  end function num_w3_probes
  
  
  
  subroutine get_probe_times (this, time)
    type(ProbeArray), intent(in) :: this
    real, dimension(:), pointer :: time
    time => this%time(1:this%n)
  end subroutine get_probe_times
  
  subroutine get_w1_probe_values (this, n, value)
    type(ProbeArray), intent(in) :: this
    integer, intent(in) :: n
    real, pointer :: value(:,:)
    if (this%w1_probe(n)%defined) then
      value => this%w1_probe(n)%value(:,1:this%n)
    else
      value => null()
    end if
  end subroutine get_w1_probe_values
  
  subroutine get_w2_probe_values (this, n, value)
    type(ProbeArray), intent(in) :: this
    integer, intent(in) :: n
    real, pointer :: value(:,:)
    if (this%w2_probe(n)%defined) then
      value => this%w2_probe(n)%value(:,1:this%n)
    else
      value => null()
    end if
  end subroutine get_w2_probe_values
    
  subroutine get_w3_probe_values (this, n, value)
    type(ProbeArray), intent(in) :: this
    integer, intent(in) :: n
    real, pointer :: value(:)
    if (this%w3_probe(n)%defined) then
      value => this%w3_probe(n)%value(1:this%n)
    else
      value => null()
    end if
  end subroutine get_w3_probe_values
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! FIND THE LOCATION OF A POINT WITHIN THE MESH
!!
!! This is EXTREMELY simple-minded and SLOW.  We have no 'fuzzy' measure, so that points
!! on the boundary may appear not to belong to any cell.  We simply search cell-by-cell.
!! A much more sophisticated procedure -- e.g. Kd-trees, etc. -- would be desirable, but
!! since we only need to locate a few points this is hopefully adequate.
!!
    
  integer function cell (d, point)
    type(simpl_mesh), intent(in) :: d
    real(kind=r8), intent(in) :: point(:)
    integer :: j
    cell = 0
    do j = 1, d%ncell
      if (outside(d%x(:,d%cnode(:,j)), point)) cycle
      cell = j
      exit
    end do
  end function cell

  logical function outside (x, point)
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8), intent(in) :: point(:)
    integer :: k, kp1
    real(kind=r8) :: p(3,4)
    outside = .true.
    p = tet_face_normal(x)  ! area-weighted outward face normals
    if (tet_volume(x) < 0.0_r8) p = -p
    do k = 1, 4
      kp1 = 1 + modulo(k,4) ! any index other than k will do
      if (dot_product(p(:,k), point - x(:,kp1)) > 0.0_r8) return
    end do
    outside = .false.
  end function outside

end module field_probes
