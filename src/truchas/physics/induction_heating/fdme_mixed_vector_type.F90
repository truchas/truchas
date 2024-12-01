
module fdme_mixed_zvector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use zvector_class
  use simpl_mesh_type
  use parallel_communication, only: global_sum
  implicit none
  private

  type, extends(zvector), public :: fdme_mixed_zvector
    type(simpl_mesh), pointer :: mesh => null()
    complex(r8), allocatable :: u(:)  ! edge-based unknowns
    complex(r8), allocatable :: p(:)  ! node-based lagrange multipliers
  contains
    !! Deferred base class procedures
    procedure :: clone1
    procedure :: clone2
    procedure :: copy_
    procedure :: setval
    procedure :: conjg1_
    procedure :: conjg2_
    procedure :: scale
    procedure :: update1_
    procedure :: update2_
    procedure :: update3_
    procedure :: update4_
    procedure :: dotc_
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

  !! Specific subroutine for the generic INIT. Initialize a FDME_MIXED_ZVECTOR
  !! object for the given simplicial MESH. The elements are initialized to 0.

  subroutine init_mesh(this, mesh)
    class(fdme_mixed_zvector), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%u(mesh%nedge), source=cmplx(0,0,kind=r8))
    allocate(this%p(mesh%nnode), source=cmplx(0,0,kind=r8))
  end subroutine

  !! Specific subroutine for the generic INIT. Initialize a FDME_MIXED_ZVECTOR
  !! object to be a clone of MOLD. The elements are initialized to 0.

  subroutine init_mold(this, mold)
    class(fdme_mixed_zvector), intent(out) :: this
    class(fdme_mixed_zvector), intent(in)  :: mold
    call init_mesh(this, mold%mesh)
  end subroutine

  subroutine gather_offp(this)
    class(fdme_mixed_zvector), intent(inout) :: this
    call this%mesh%edge_imap%gather_offp(this%u)
    call this%mesh%node_imap%gather_offp(this%p)
  end subroutine

  subroutine clone1(this, clone)
    class(fdme_mixed_zvector), intent(in) :: this
    class(zvector), allocatable, intent(out) :: clone
    allocate(clone, source=this)
  end subroutine

  subroutine clone2(this, clone, n)
    class(fdme_mixed_zvector), intent(in)  :: this
    class(zvector), allocatable, intent(out) :: clone(:)
    integer, intent(in) :: n
    allocate(clone(n), source=this)
  end subroutine

  subroutine copy_(dest, src)
    class(fdme_mixed_zvector), intent(inout) :: dest
    class(zvector), intent(in) :: src
    select type (src)
    class is (fdme_mixed_zvector)
      dest%u(:) = src%u
      dest%p(:) = src%p
    end select
  end subroutine

  subroutine setval(this, val)
    class(fdme_mixed_zvector), intent(inout) :: this
    complex(r8), intent(in) :: val
    this%u = val
    this%p = val
  end subroutine

  subroutine conjg1_(this)
    class(fdme_mixed_zvector), intent(inout) :: this
    integer :: j
    do j = 1, this%mesh%nedge_onP
      this%u(j)%im = -this%u(j)%im
    end do
    do j = 1, this%mesh%nnode_onP
      this%p(j)%im = -this%p(j)%im
    end do
  end subroutine

  subroutine conjg2_(this, src)
    class(fdme_mixed_zvector), intent(inout) :: this
    class(zvector), intent(in) :: src
    integer :: j
    select type (src)
    type is (fdme_mixed_zvector)
      do j = 1, this%mesh%nedge_onP
        this%u(j) = conjg(src%u(j))
      end do
      do j = 1, this%mesh%nnode_onP
        this%p(j) = conjg(src%p(j))
      end do
    end do
  end subroutine

  subroutine scale(this, a)
    class(fdme_mixed_zvector), intent(inout) :: this
    complex(r8), intent(in) :: a
    integer :: j
    do j = 1, this%mesh%nedge_onP
      this%u(j) = a * this%u(j)
    end do
    do j = 1, this%mesh%nnode_onP
      this%p(j) = a * this%p(j)
    end do
  end subroutine

  !! Conventional SAXPY procedure: y <-- a*x + y
  subroutine update1_(this, a, x)
    class(fdme_mixed_zvector), intent(inout) :: this
    class(zvector), intent(in) :: x
    complex(r8), intent(in) :: a
    integer :: j
    select type (x)
    class is (fdme_mixed_zvector)
      do j = 1, this%mesh%nedge_onP
        this%u(j) = a * x%u(j) + this%u(j)
      end do
      do j = 1, this%mesh%nnode_onP
        this%p(j) = a * x%p(j) + this%p(j)
      end do
    end select
  end subroutine

  !! SAXPY-like procedure: y <-- a*x + b*y
  subroutine update2_(this, a, x, b)
    class(fdme_mixed_zvector), intent(inout) :: this
    class(zvector), intent(in) :: x
    complex(r8), intent(in) :: a, b
    integer :: j
    select type (x)
    class is (fdme_mixed_zvector)
      if (b == 0) then
        do j = 1, this%mesh%nedge_onP
          this%u(j) = a * x%u(j)
        end do
        do j = 1, this%mesh%nnode_onP
          this%p(j) = a * x%p(j)
        end do
      else
        do j = 1, this%mesh%nedge_onP
          this%u(j) = a * x%u(j) + b * this%u(j)
        end do
        do j = 1, this%mesh%nnode_onP
          this%p(j) = a * x%p(j) + b * this%p(j)
        end do
      end if
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + z
  subroutine update3_(this, a, x, b, y)
    class(fdme_mixed_zvector), intent(inout) :: this
    class(zvector), intent(in) :: x, y
    complex(r8), intent(in) :: a, b
    select type (x)
    class is (fdme_mixed_zvector)
      select type (y)
      class is (fdme_mixed_zvector)
        do j = 1, this%mesh%nedge_onP
          this%u(j) = a * x%u(j) + b * y%u(j) + this%u(j)
        end do
        do j = 1, this%mesh%nnode_onP
          this%p(j) = a * x%p(j) + b * y%p(j) + this%p(j)
        end do
      end select
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + c*z
  subroutine update4_(this, a, x, b, y, c)
    class(fdme_mixed_zvector), intent(inout) :: this
    class(zvector), intent(in) :: x, y
    complex(r8), intent(in) :: a, b, c
    select type (x)
    class is (fdme_mixed_zvector)
      select type (y)
      class is (fdme_mixed_zvector)
        if (c == 0) then
          do j = 1, this%mesh%nedge_onP
            this%u(j) = a * x%u(j) + b * y%u(j)
          end do
          do j = 1, this%mesh%nnode_onP
            this%p(j) = a * x%p(j) + b * y%p(j)
          end do
        else
          do j = 1, this%mesh%nedge_onP
            this%u(j) = a * x%u(j) + b * y%u(j) + c * this%u(j)
          end do
          do j = 1, this%mesh%nnode_onP
            this%p(j) = a * x%p(j) + b * y%p(j) + c * this%p(j)
          end do
        end if
      end select
    end select
  end subroutine

  function dotc_(x, y) result(dp)
    class(fdme_mixed_zvector), intent(in) :: x
    class(zvector), intent(in) :: y
    complex(r8) :: dp
    integer :: j
    select type (y)
    class is (fdme_mixed_zvector)
      dp = 0
      do j = 1, x%mesh%nedge_onP
        dp = dp + conjg(x%u(j)) * y%u(j)
      end do
      do j = 1, x%mesh%nnode_onP
        dp = dp + conjg(x%p(j)) * y%p(j)
      end do
      dp = global_sum(dp)
    end select
  end function

  real(r8) function norm1_(this) result(norm)
    use parallel_communication, only: global_sum
    class(fdme_mixed_zvector), intent(in) :: this
    integer :: j
    norm = 0
    do j = 1, this%mesh%nedge_onP
      norm = norm + abs(this%u(j))
    end do
    do j = 1, this%mesh%nnode_onP
      norm = norm + abs(this%p(j))
    end do
    norm = global_sum(norm)
  end function

  real(r8) function norm2_(this) result(norm)
    use parallel_communication, only: global_sum
    class(fdme_mixed_zvector), intent(in) :: this
    integer :: j
    norm = 0
    do j = 1, this%mesh%nedge_onP
      norm = norm + this%u(j)%re**2 + this%u(j)%im**2
    end do
    do j = 1, this%mesh%nnode_onP
      norm = norm + this%p(j)%re**2 + this%p(j)%im**2
    end do
    norm = sqrt(global_sum(norm))
  end function

  real(r8) function norm_max_(this) result(norm)
    use parallel_communication, only: global_maxval
    class(fdme_mixed_zvector), intent(in) :: this
    integer :: j
    norm = 0
    do j = 1, this%mesh%nedge_onP
      norm = max(norm, abs(this%u(j))
    end do
    do j = 1, this%mesh%nnode_onP
      norm = max(norm, abs(this%p(j))
    end do
    norm = global_maxval(norm)
  end function

  function checksum(this, full) result(string)
    use md5_hash_type
    class(fdme_mixed_zvector), intent(in) :: this
    logical, intent(in), optional :: full ! default is FALSE
    character(:), allocatable :: string
    type(md5_hash) :: hash
    logical :: strict
    strict = .true.
    if (present(full)) strict = .not.full
    if (strict) then
      !NB: newer versions of md5_hash_type support complex data directly
      call hash%update(this%u(:this%mesh%nedge_onP)%re)
      call hash%update(this%u(:this%mesh%nedge_onP)%im)
      call hash%update(this%p(:this%mesh%nnode_onP)%re)
      call hash%update(this%p(:this%mesh%nnode_onP)%im)
    else
      call hash%update(this%u%re)
      call hash%update(this%u%im)
      call hash%update(this%p%re)
      call hash%update(this%p%im)
    end if
    string = hash%hexdigest()
  end function

end module fdme_mixed_zvector_type
