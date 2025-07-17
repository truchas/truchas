
module alloy_vector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_class
  use unstr_mesh_type
  use parallel_communication, only: global_sum
  implicit none
  private

  type, extends(vector), public :: alloy_vector
    type(unstr_mesh), pointer :: mesh => null()
    real(r8), allocatable :: lf(:), hc(:), tc(:), tf(:)
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
    procedure :: gather_offp
  end type

contains

  !! Specific subroutine for the generic INIT. Initialize a alloy_vector object
  !! for the given unstructured MESH. The elements are initialized to 0.

  subroutine init_mesh(this, mesh)
    class(alloy_vector), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%lf(mesh%ncell), this%hc(mesh%ncell), this%tc(mesh%ncell), this%tf(mesh%nface))
    call this%setval(0.0_r8)
  end subroutine

  !! Specific subroutine for the generic INIT. Initialize a alloy_vector object
  !! to be a clone of MOLD. The elements are initialized to zero.

  subroutine init_mold(this, mold)
    class(alloy_vector), intent(out) :: this
    class(alloy_vector), intent(in)  :: mold
    call init_mesh(this, mold%mesh)
  end subroutine

  subroutine gather_offp(this)
    class(alloy_vector), intent(inout) :: this
    call this%mesh%cell_imap%gather_offp(this%lf)
    call this%mesh%cell_imap%gather_offp(this%hc)
    call this%mesh%cell_imap%gather_offp(this%tc)
    call this%mesh%face_imap%gather_offp(this%tf)
  end subroutine

  subroutine clone1(this, clone)
    class(alloy_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone
    allocate(clone, source=this) ! easy, but with unwanted copy of data too
  end subroutine

  subroutine clone2(this, clone, n)
    class(alloy_vector), intent(in)  :: this
    class(vector), allocatable, intent(out) :: clone(:)
    integer, intent(in) :: n
    allocate(clone(n), source=this) ! easy, but with unwanted copy of data too
  end subroutine

  subroutine copy_(dest, src)
    class(alloy_vector), intent(inout) :: dest
    class(vector), intent(in) :: src
    select type (src)
    class is (alloy_vector)
      dest%lf(:) = src%lf
      dest%hc(:) = src%hc
      dest%tc(:) = src%tc
      dest%tf(:) = src%tf
    end select
  end subroutine

  subroutine setval(this, val)
    class(alloy_vector), intent(inout) :: this
    real(r8), intent(in) :: val
    this%lf = val
    this%hc = val
    this%tc = val
    this%tf = val
  end subroutine

  subroutine scale(this, a)
    class(alloy_vector), intent(inout) :: this
    real(r8), intent(in) :: a
    this%lf = a * this%lf
    this%hc = a * this%hc
    this%tc = a * this%tc
    this%tf = a * this%tf
  end subroutine

  !! Conventional SAXPY procedure: y <-- a*x + y
  subroutine update1_(this, a, x)
    class(alloy_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a
    select type (x)
    class is (alloy_vector)
      this%lf = a * x%lf + this%lf
      this%hc = a * x%hc + this%hc
      this%tc = a * x%tc + this%tc
      this%tf = a * x%tf + this%tf
    end select
  end subroutine

  !! SAXPY-like procedure: y <-- a*x + b*y
  subroutine update2_(this, a, x, b)
    class(alloy_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a, b
    select type (x)
    class is (alloy_vector)
      this%lf = a * x%lf + b * this%lf
      this%hc = a * x%hc + b * this%hc
      this%tc = a * x%tc + b * this%tc
      this%tf = a * x%tf + b * this%tf
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + z
  subroutine update3_(this, a, x, b, y)
    class(alloy_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b
    select type (x)
    class is (alloy_vector)
      select type (y)
      class is (alloy_vector)
        this%lf = a * x%lf + b * y%lf + this%lf
        this%hc = a * x%hc + b * y%hc + this%hc
        this%tc = a * x%tc + b * y%tc + this%tc
        this%tf = a * x%tf + b * y%tf + this%tf
      end select
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + c*z
  subroutine update4_(this, a, x, b, y, c)
    class(alloy_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b, c
    select type (x)
    class is (alloy_vector)
      select type (y)
      class is (alloy_vector)
        this%lf = a * x%lf + b * y%lf + c * this%lf
        this%hc = a * x%hc + b * y%hc + c * this%hc
        this%tc = a * x%tc + b * y%tc + c * this%tc
        this%tf = a * x%tf + b * y%tf + c * this%tf
      end select
    end select
  end subroutine

  function dot_(x, y) result(dp)
    class(alloy_vector), intent(in) :: x
    class(vector), intent(in) :: y
    real(r8) :: dp
    integer :: j
    select type (y)
    class is (alloy_vector)
      dp = 0
      do j = 1, x%mesh%ncell_onP
        dp = dp + x%lf(j) * y%lf(j)
      end do
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
    class(alloy_vector), intent(in) :: this
    real(r8) :: norm2_
    norm2_ = norm2(this%lf(:this%mesh%ncell_onP))**2 + &
             norm2(this%hc(:this%mesh%ncell_onP))**2 + &
             norm2(this%tc(:this%mesh%ncell_onP))**2 + &
             norm2(this%tf(:this%mesh%nface_onP))**2
    norm2_ = sqrt(global_sum(norm2_))
  end function

  function checksum(this, full) result(string)
    use md5_hash_type
    class(alloy_vector), intent(in) :: this
    logical, intent(in), optional :: full ! default is FALSE
    character(:), allocatable :: string
    type(md5_hash) :: hash
    logical :: strict
    strict = .true.
    if (present(full)) strict = .not.full
    if (strict) then
      call hash%update(this%lf(:this%mesh%ncell_onP))
      call hash%update(this%hc(:this%mesh%ncell_onP))
      call hash%update(this%tc(:this%mesh%ncell_onP))
      call hash%update(this%tf(:this%mesh%nface_onP))
    else
      call hash%update(this%lf(:this%mesh%ncell))
      call hash%update(this%hc(:this%mesh%ncell))
      call hash%update(this%tc(:this%mesh%ncell))
      call hash%update(this%tf(:this%mesh%nface))
    end if
    string = hash%hexdigest()
  end function

end module alloy_vector_type
