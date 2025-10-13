!!
!! IMAP_ZVECTOR
!!
!! This defines an implementation of the abstract ZVECTOR class that stores
!! a single distributed complex vector whose distribution is described by an
!! INDEX_MAP object.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module imap_zvector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use zvector_class
  use index_map_type
  implicit none
  private

  type, extends(zvector), public :: imap_zvector
    type(index_map), pointer :: imap => null() ! reference only
    complex(r8), allocatable :: v(:)
  contains
    !! Deferred base class procedures
    procedure :: clone1
    procedure :: clone2
    procedure :: copy_
    procedure :: setvalr_
    procedure :: setvalc_
    procedure :: conjg1_
    procedure :: conjg2_
    procedure :: scale
    procedure :: update1_
    procedure :: update2_
    procedure :: update3_
    procedure :: update4_
    procedure :: dotc_
    procedure :: dotu_
    procedure :: norm1 => norm1_
    procedure :: norm2 => norm2_
    procedure :: norm_max => norm_max_
    procedure :: checksum
    !! Additional procedures specific to this type
    generic :: init => init_imap, init_mold
    procedure, private :: init_imap, init_mold
    procedure :: gather_offp
  end type

contains

  !! Specific subroutine for the generic INIT. Initialize a IMAP_ZVECTOR object
  !! for the given index_map object IMAP. The elements are initialized to 0.

  subroutine init_imap(this, imap)
    class(imap_zvector), intent(out) :: this
    type(index_map), intent(in), target :: imap
    this%imap => imap
    allocate(this%v(imap%local_size))
    this%v = 0
  end subroutine

  !! Specific subroutine for the generic INIT. Initialize a IMAP_VECTOR object
  !! to be a clone of MOLD. The elements are initialized to zero.

  subroutine init_mold(this, mold)
    class(imap_zvector), intent(out) :: this
    class(imap_zvector), intent(in)  :: mold
    call init_imap(this, mold%imap)
  end subroutine

  subroutine gather_offp(this)
    class(imap_zvector), intent(inout) :: this
    call this%imap%gather_offp(this%v)
  end subroutine

  subroutine clone1(this, clone)
    class(imap_zvector), intent(in) :: this
    class(zvector), allocatable, intent(out) :: clone
    allocate(clone, source=this)
  end subroutine

  subroutine clone2(this, clone, n)
    class(imap_zvector), intent(in)  :: this
    class(zvector), allocatable, intent(out) :: clone(:)
    integer, intent(in) :: n
    allocate(clone(n), source=this)
  end subroutine

  subroutine copy_(dest, src)
    class(imap_zvector), intent(inout) :: dest
    class(zvector), intent(in) :: src
    select type (src)
    class is (imap_zvector)
      dest%v(:) = src%v
    end select
  end subroutine

  subroutine setvalr_(this, val)
    class(imap_zvector), intent(inout) :: this
    real(r8), intent(in) :: val
    this%v = val
  end subroutine

  subroutine setvalc_(this, val)
    class(imap_zvector), intent(inout) :: this
    complex(r8), intent(in) :: val
    this%v = val
  end subroutine

  subroutine conjg1_(this)
    class(imap_zvector), intent(inout) :: this
    integer :: j
    do j = 1, this%imap%onp_size
      this%v(j)%im = -this%v(j)%im
    end do
  end subroutine

  subroutine conjg2_(this, src)
    class(imap_zvector), intent(inout) :: this
    class(zvector), intent(in) :: src
    integer :: j
    select type (src)
    type is (imap_zvector)
      do j = 1, this%imap%onp_size
        this%v(j) = conjg(src%v(j))
      end do
    end select
  end subroutine

  subroutine scale(this, a)
    class(imap_zvector), intent(inout) :: this
    complex(r8), intent(in) :: a
    integer :: j
    do j = 1, this%imap%onp_size
      this%v(j) = a * this%v(j)
    end do
  end subroutine

  !! Conventional SAXPY procedure: y <-- a*x + y
  subroutine update1_(this, a, x)
    class(imap_zvector), intent(inout) :: this
    class(zvector), intent(in) :: x
    complex(r8), intent(in) :: a
    integer :: j
    select type (x)
    class is (imap_zvector)
      do j = 1, this%imap%onp_size
        this%v(j) = a * x%v(j) + this%v(j)
      end do
    end select
  end subroutine

  !! SAXPY-like procedure: y <-- a*x + b*y
  subroutine update2_(this, a, x, b)
    class(imap_zvector), intent(inout) :: this
    class(zvector), intent(in) :: x
    complex(r8), intent(in) :: a, b
    integer :: j
    select type (x)
    class is (imap_zvector)
      if (b == 0) then
        do j = 1, this%imap%onp_size
          this%v(j) = a * x%v(j)
        end do
      else
        do j = 1, this%imap%onp_size
          this%v(j) = a * x%v(j) + b * this%v(j)
        end do
      end if
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + z
  subroutine update3_(this, a, x, b, y)
    class(imap_zvector), intent(inout) :: this
    class(zvector), intent(in) :: x, y
    complex(r8), intent(in) :: a, b
    integer :: j
    select type (x)
    class is (imap_zvector)
      select type (y)
      class is (imap_zvector)
        do j = 1, this%imap%onp_size
          this%v(j) = a * x%v(j) + b * y%v(j) + this%v(j)
        end do
      end select
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + c*z
  subroutine update4_(this, a, x, b, y, c)
    class(imap_zvector), intent(inout) :: this
    class(zvector), intent(in) :: x, y
    complex(r8), intent(in) :: a, b, c
    integer :: j
    select type (x)
    class is (imap_zvector)
      select type (y)
      class is (imap_zvector)
        if (c == 0) then
          do j = 1, this%imap%onp_size
            this%v(j) = a * x%v(j) + b * y%v(j)
          end do
        else
          do j = 1, this%imap%onp_size
            this%v(j) = a * x%v(j) + b * y%v(j) + c * this%v(j)
          end do
        end if
      end select
    end select
  end subroutine

  function dotc_(x, y) result(dp)
    use parallel_communication, only: global_sum
    class(imap_zvector), intent(in) :: x
    class(zvector), intent(in) :: y
    complex(r8) :: dp
    integer :: j
    select type (y)
    class is (imap_zvector)
      dp = 0.0_r8
      do j = 1, x%imap%onp_size
        dp = dp + conjg(x%v(j)) * y%v(j)
      end do
      dp = global_sum(dp)
    end select
  end function

  function dotu_(x, y) result(dp)
    use parallel_communication, only: global_sum
    class(imap_zvector), intent(in) :: x
    class(zvector), intent(in) :: y
    complex(r8) :: dp
    integer :: j
    select type (y)
    class is (imap_zvector)
      dp = 0.0_r8
      do j = 1, x%imap%onp_size
        dp = dp + x%v(j) * y%v(j)
      end do
      dp = global_sum(dp)
    end select
  end function

  real(r8) function norm1_(this) result(norm)
    use parallel_communication, only: global_sum
    class(imap_zvector), intent(in) :: this
    integer :: j
    norm = 0
    do j = 1, this%imap%onp_size
      norm = norm + abs(this%v(j))
    end do
    norm = global_sum(norm)
  end function

  real(r8) function norm2_(this) result(norm)
    use parallel_communication, only: global_sum
    class(imap_zvector), intent(in) :: this
    integer :: j
    norm = 0
    do j = 1, this%imap%onp_size
      norm = norm + this%v(j)%re**2 + this%v(j)%im**2
    end do
    norm = sqrt(global_sum(norm))
  end function

  real(r8) function norm_max_(this) result(norm)
    use parallel_communication, only: global_maxval
    class(imap_zvector), intent(in) :: this
    integer :: j
    norm = 0
    do j = 1, this%imap%onp_size
      norm = max(norm, abs(this%v(j)))
    end do
    norm = global_maxval(norm)
  end function

  function checksum(this, full) result(string)
    use md5_hash_type
    class(imap_zvector), intent(in) :: this
    logical, intent(in), optional :: full ! default is FALSE
    character(:), allocatable :: string
    type(md5_hash) :: hash
    logical :: strict
    strict = .true.
    if (present(full)) strict = .not.full
    if (strict) then
      call hash%update(this%v(:this%imap%onp_size))
    else
      call hash%update(this%v)
    end if
    string = hash%hexdigest()
  end function

end module imap_zvector_type
