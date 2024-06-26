
module fdme_vector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_class
  use simpl_mesh_type
  use parallel_communication, only: global_sum
  implicit none
  private

  type, extends(vector), public :: fdme_vector
    type(simpl_mesh), pointer :: mesh => null()
    real(r8), allocatable :: array(:,:)
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
    procedure :: norm1 => norm1_
    procedure :: norm2 => norm2_
    procedure :: norm_max => norm_max_
    procedure :: checksum
    !! Additional procedures specific to this type
    generic :: init => init_mesh, init_mold
    procedure, private :: init_mesh, init_mold
    procedure :: gather_offp
  end type

contains

  !! Specific subroutine for the generic INIT. Initialize a FDME_VECTOR object
  !! for the given unstructured MESH. The elements are initialized to 0.

  subroutine init_mesh(this, mesh)
    class(fdme_vector), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%array(2,mesh%nedge), source=0.0_r8)
  end subroutine

  !! Specific subroutine for the generic INIT. Initialize a FDME_VECTOR object
  !! to be a clone of MOLD. The elements are initialized to zero.

  subroutine init_mold(this, mold)
    class(fdme_vector), intent(out) :: this
    class(fdme_vector), intent(in)  :: mold
    call init_mesh(this, mold%mesh)
  end subroutine

  subroutine gather_offp(this)
    class(fdme_vector), intent(inout) :: this
    call this%mesh%edge_imap%gather_offp(this%array)
  end subroutine

  subroutine clone1(this, clone)
    class(fdme_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone
    allocate(clone, source=this) ! easy, but with unwanted copy of data too
  end subroutine

  subroutine clone2(this, clone, n)
    class(fdme_vector), intent(in)  :: this
    class(vector), allocatable, intent(out) :: clone(:)
    integer, intent(in) :: n
    allocate(clone(n), source=this) ! easy, but with unwanted copy of data too
  end subroutine

  subroutine copy_(dest, src)
    class(fdme_vector), intent(inout) :: dest
    class(vector), intent(in) :: src
    select type (src)
    class is (fdme_vector)
      dest%array(:,:) = src%array
    end select
  end subroutine

  subroutine setval(this, val)
    class(fdme_vector), intent(inout) :: this
    real(r8), intent(in) :: val
    this%array = val
  end subroutine

  subroutine scale(this, a)
    class(fdme_vector), intent(inout) :: this
    real(r8), intent(in) :: a
    this%array = a * this%array
  end subroutine

  !! Conventional SAXPY procedure: y <-- a*x + y
  subroutine update1_(this, a, x)
    class(fdme_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a
    select type (x)
    class is (fdme_vector)
      this%array = a * x%array + this%array
    end select
  end subroutine

  !! SAXPY-like procedure: y <-- a*x + b*y
  subroutine update2_(this, a, x, b)
    class(fdme_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a, b
    select type (x)
    class is (fdme_vector)
      this%array = a * x%array + b * this%array
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + z
  subroutine update3_(this, a, x, b, y)
    class(fdme_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b
    select type (x)
    class is (fdme_vector)
      select type (y)
      class is (fdme_vector)
        this%array = a * x%array + b * y%array + this%array
      end select
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + c*z
  subroutine update4_(this, a, x, b, y, c)
    class(fdme_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b, c
    select type (x)
    class is (fdme_vector)
      select type (y)
      class is (fdme_vector)
        this%array = a * x%array + b * y%array + c * this%array
      end select
    end select
  end subroutine

  function dot_(x, y) result(dp)
    class(fdme_vector), intent(in) :: x
    class(vector), intent(in) :: y
    real(r8) :: dp
    integer :: j
    select type (y)
    class is (fdme_vector)
      dp = 0.0_r8
      do j = 1, x%mesh%nedge_onP
        dp = dp + x%array(1,j) * y%array(1,j) + x%array(2,j) * y%array(2,j)
      end do
      dp = global_sum(dp)
    end select
  end function

  real(r8) function norm1_(this) result(norm)
    use parallel_communication, only: global_sum
    class(fdme_vector), intent(in) :: this
    norm = sum(abs(this%array(:,:this%mesh%nedge_onP)))
    norm = global_sum(norm)
  end function

  real(r8) function norm2_(this) result(norm)
    use parallel_communication, only: global_sum
    class(fdme_vector), intent(in) :: this
    norm = norm2(this%array(:,:this%mesh%nedge_onP))**2
    norm = sqrt(global_sum(norm))
  end function

  real(r8) function norm_max_(this) result(norm)
    use parallel_communication, only: global_maxval
    class(fdme_vector), intent(in) :: this
    norm = maxval(abs(this%array(:,:this%mesh%nedge_onP)))
    norm = global_maxval(norm)
  end function

  function checksum(this, full) result(string)
    use md5_hash_type
    class(fdme_vector), intent(in) :: this
    logical, intent(in), optional :: full ! default is FALSE
    character(:), allocatable :: string
    type(md5_hash) :: hash
    logical :: strict
    strict = .true.
    if (present(full)) strict = .not.full
    if (strict) then
      call hash%update(this%array(:,:this%mesh%nedge_onP))
    else
      call hash%update(this%array)
    end if
    string = hash%hexdigest()
  end function

end module fdme_vector_type
