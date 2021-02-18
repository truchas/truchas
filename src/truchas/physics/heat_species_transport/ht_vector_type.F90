
module ht_vector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_class
  use unstr_mesh_type
  use parallel_communication, only: global_sum
  implicit none
  private

  type, extends(vector), public :: ht_vector
    type(unstr_mesh), pointer :: mesh => null()
    real(r8), allocatable :: hc(:), tc(:), tf(:)
  contains
    !! Deferred base class procedures
    procedure :: clone1
    procedure :: clone2
    procedure :: copy_
    procedure :: setval
    procedure :: scale
    procedure :: update1_
    procedure :: update2_
    procedure :: update3_
    procedure :: update4_
    procedure :: dot_
    procedure :: norm2 => norm2_
    procedure :: checksum
    !! Additional procedures specific to this type
    generic :: init => init_mesh, init_mold
    procedure, private :: init_mesh, init_mold
    procedure :: gather_boundary => ht_gather_boundary
  end type

contains

  !! Specific subroutine for the generic INIT. Initialize a HT_VECTOR object
  !! for the given unstructured MESH. The elements are initialized to 0.

  subroutine init_mesh(this, mesh)
    class(ht_vector), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%hc(mesh%ncell), this%tc(mesh%ncell), this%tf(mesh%nface))
    call setval(this, 0.0_r8)
  end subroutine

  !! Specific subroutine for the generic INIT. Initialize a HT_VECTOR object
  !! to be a clone of MOLD. The elements are initialized to zero.

  subroutine init_mold(this, mold)
    class(ht_vector), intent(out) :: this
    class(ht_vector), intent(in)  :: mold
    call init_mesh(this, mold%mesh)
  end subroutine

  subroutine ht_gather_boundary(this)
    use index_partitioning, only: gather_boundary
    class(ht_vector), intent(inout) :: this
    call gather_boundary(this%mesh%cell_ip, this%hc)
    call gather_boundary(this%mesh%cell_ip, this%tc)
    call gather_boundary(this%mesh%face_ip, this%tf)
  end subroutine

  subroutine clone1(this, clone)
    class(ht_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone
    type(ht_vector), allocatable :: tmp
    allocate(tmp)
    tmp%mesh => this%mesh
    allocate(tmp%hc, mold=this%hc)
    allocate(tmp%tc, mold=this%tc)
    allocate(tmp%tf, mold=this%tf)
    call move_alloc(tmp, clone)
  end subroutine

  subroutine clone2(this, clone, n)
    class(ht_vector), intent(in)  :: this
    class(vector), allocatable, intent(out) :: clone(:)
    type(ht_vector), allocatable :: tmp(:)
    integer, intent(in) :: n
    integer :: j
    allocate(tmp(n))
    do j = 1, n
      tmp(j)%mesh => this%mesh
      allocate(tmp(j)%hc, mold=this%hc)
      allocate(tmp(j)%tc, mold=this%tc)
      allocate(tmp(j)%tf, mold=this%tf)
    end do
    call move_alloc(tmp, clone)
  end subroutine

  subroutine copy_(dest, src)
    class(ht_vector), intent(inout) :: dest
    class(vector), intent(in) :: src
    select type (src)
    class is (ht_vector)
      dest%hc(:) = src%hc
      dest%tc(:) = src%tc
      dest%tf(:) = src%tf
    end select
  end subroutine

  subroutine setval(this, val)
    class(ht_vector), intent(inout) :: this
    real(r8), intent(in) :: val
    this%hc = val
    this%tc = val
    this%tf = val
  end subroutine

  subroutine scale(this, a)
    class(ht_vector), intent(inout) :: this
    real(r8), intent(in) :: a
    this%hc = a * this%hc
    this%tc = a * this%tc
    this%tf = a * this%tf
  end subroutine

  !! Conventional SAXPY procedure: y <-- a*x + y
  subroutine update1_(this, a, x)
    class(ht_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a
    select type (x)
    class is (ht_vector)
      this%hc = a * x%hc + this%hc
      this%tc = a * x%tc + this%tc
      this%tf = a * x%tf + this%tf
    end select
  end subroutine

  !! SAXPY-like procedure: y <-- a*x + b*y
  subroutine update2_(this, a, x, b)
    class(ht_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a, b
    select type (x)
    class is (ht_vector)
      this%hc = a * x%hc + b * this%hc
      this%tc = a * x%tc + b * this%tc
      this%tf = a * x%tf + b * this%tf
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + z
  subroutine update3_(this, a, x, b, y)
    class(ht_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b
    select type (x)
    class is (ht_vector)
      select type (y)
      class is (ht_vector)
        this%hc = a * x%hc + b * y%hc + this%hc
        this%tc = a * x%tc + b * y%tc + this%tc
        this%tf = a * x%tf + b * y%tf + this%tf
      end select
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + c*z
  subroutine update4_(this, a, x, b, y, c)
    class(ht_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b, c
    select type (x)
    class is (ht_vector)
      select type (y)
      class is (ht_vector)
        this%hc = a * x%hc + b * y%hc + c * this%hc
        this%tc = a * x%tc + b * y%tc + c * this%tc
        this%tf = a * x%tf + b * y%tf + c * this%tf
      end select
    end select
  end subroutine

  function dot_(x, y) result(dp)
    class(ht_vector), intent(in) :: x
    class(vector), intent(in) :: y
    real(r8) :: dp
    integer :: j
    select type (y)
    class is (ht_vector)
      dp = 0
      do j = 1, x%mesh%ncell_onP
        dp = dp + x%hc(j) * y%hc(j)
      end do
      do j = 1, x%mesh%ncell_onP
        dp = dp + x%tc(j) * y%tc(j)
      end do
      do j = 1, x%mesh%nface_onP
        dp = dp + x%tf(j) * y%tf(j)
      end do
      dp = global_sum(dp)
    end select
  end function

  function norm2_(this)
    use parallel_communication, only: global_sum
    class(ht_vector), intent(in) :: this
    real(r8) :: norm2_
    norm2_ = norm2(this%hc(:this%mesh%ncell_onP))**2 + &
             norm2(this%tc(:this%mesh%ncell_onP))**2 + &
             norm2(this%tf(:this%mesh%nface_onP))**2
    norm2_ = sqrt(global_sum(norm2_))
  end function

  function checksum(this, full) result(string)
    use md5_hash_type
    class(ht_vector), intent(in) :: this
    logical, intent(in), optional :: full ! default is FALSE
    character(:), allocatable :: string
    type(md5_hash) :: hash
    logical :: strict
    strict = .true.
    if (present(full)) strict = .not.full
    if (strict) then
      call hash%update(this%hc(:this%mesh%ncell_onP))
      call hash%update(this%tc(:this%mesh%ncell_onP))
      call hash%update(this%tf(:this%mesh%nface_onP))
    else
      call hash%update(this%hc(:this%mesh%ncell))
      call hash%update(this%tc(:this%mesh%ncell))
      call hash%update(this%tf(:this%mesh%nface))
    end if
    string = hash%hexdigest()
  end function

end module ht_vector_type
