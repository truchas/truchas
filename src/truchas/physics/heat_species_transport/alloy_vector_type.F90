
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
    real(r8), allocatable :: lsf(:,:)
    integer :: num_comp
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

  !! Specific subroutine for the generic INIT. Initialize a alloy_vector object
  !! for the given unstructured MESH. The elements are initialized to 0.

  subroutine init_mesh(this, mesh, num_comp)
    class(alloy_vector), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in), optional :: num_comp
    this%mesh => mesh
    allocate(this%lf(mesh%ncell), this%hc(mesh%ncell), this%tc(mesh%ncell), this%tf(mesh%nface))
    this%num_comp = 0
    if (present(num_comp)) this%num_comp = num_comp
    if (this%num_comp > 0) allocate(this%lsf(this%num_comp,mesh%ncell))
    call this%setval(0.0_r8)
  end subroutine

  !! Specific subroutine for the generic INIT. Initialize a alloy_vector object
  !! to be a clone of MOLD. The elements are initialized to zero.

  subroutine init_mold(this, mold)
    class(alloy_vector), intent(out) :: this
    class(alloy_vector), intent(in)  :: mold
    call init_mesh(this, mold%mesh, mold%num_comp)
  end subroutine

  subroutine gather_offp(this)
    class(alloy_vector), intent(inout) :: this
    call this%mesh%cell_imap%gather_offp(this%lf)
    call this%mesh%cell_imap%gather_offp(this%hc)
    call this%mesh%cell_imap%gather_offp(this%tc)
    call this%mesh%face_imap%gather_offp(this%tf)
    if (allocated(this%lsf)) call this%mesh%cell_imap%gather_offp(this%lsf)
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
      if (allocated(dest%lsf)) dest%lsf(:,:) = src%lsf
    end select
  end subroutine

  subroutine setval(this, val)
    class(alloy_vector), intent(inout) :: this
    real(r8), intent(in) :: val
    this%lf = val
    this%hc = val
    this%tc = val
    this%tf = val
    if (allocated(this%lsf)) this%lsf = val
  end subroutine

  subroutine scale(this, a)
    class(alloy_vector), intent(inout) :: this
    real(r8), intent(in) :: a
    this%lf = a * this%lf
    this%hc = a * this%hc
    this%tc = a * this%tc
    this%tf = a * this%tf
    if (allocated(this%lsf)) this%lsf = a * this%lsf
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
      if (allocated(this%lsf)) this%lsf = a * x%lsf + this%lsf
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
      if (allocated(this%lsf)) this%lsf = a * x%lsf + b * this%lsf
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
        if (allocated(this%lsf)) this%lsf = a * x%lsf + b * y%lsf + this%lsf
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
        if (allocated(this%lsf)) this%lsf = a * x%lsf + b * y%lsf + c * this%lsf
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
      if (allocated(x%lsf)) then
        do j = 1, x%mesh%ncell_onP
          dp = dp + dot_product(x%lsf(:,j), y%lsf(:,j))
        end do
      end if
      dp = global_sum(dp)
    end select
  end function

  real(r8) function norm1_(this) result(norm)
    use parallel_communication, only: global_sum
    class(alloy_vector), intent(in) :: this
    integer :: j, k
    norm = 0
    do j = 1, this%mesh%ncell_onP
      norm = norm + abs(this%lf(j))
    end do
    do j = 1, this%mesh%ncell_onP
      norm = norm + abs(this%hc(j))
    end do
    do j = 1, this%mesh%ncell_onP
      norm = norm + abs(this%tc(j))
    end do
    do j = 1, this%mesh%nface_onP
      norm = norm + abs(this%tf(j))
    end do
    if (allocated(this%lsf)) then
      do j = 1, this%mesh%ncell_onP
        do k = 1, size(this%lsf,dim=1)
          norm = norm + abs(this%lsf(k,j))
        end do
      end do
    end if
    norm = global_sum(norm)
  end function

  function norm2_(this)
    use parallel_communication, only: global_sum
    class(alloy_vector), intent(in) :: this
    real(r8) :: norm2_
    norm2_ = norm2(this%lf(:this%mesh%ncell_onP))**2 + &
             norm2(this%hc(:this%mesh%ncell_onP))**2 + &
             norm2(this%tc(:this%mesh%ncell_onP))**2 + &
             norm2(this%tf(:this%mesh%nface_onP))**2
    if (allocated(this%lsf)) then
      norm2_ = norm2_ + norm2(this%lsf(:,:this%mesh%ncell_onP))**2
    end if
    norm2_ = sqrt(global_sum(norm2_))
  end function

  real(r8) function norm_max_(this) result(norm)
    use parallel_communication, only: global_maxval
    class(alloy_vector), intent(in) :: this
    integer :: j, k
    norm = 0
    do j = 1, this%mesh%ncell_onP
      norm = max(norm, abs(this%lf(j)))
    end do
    do j = 1, this%mesh%ncell_onP
      norm = max(norm, abs(this%hc(j)))
    end do
    do j = 1, this%mesh%ncell_onP
      norm = max(norm, abs(this%tc(j)))
    end do
    do j = 1, this%mesh%nface_onP
      norm = max(norm, abs(this%tf(j)))
    end do
    if (allocated(this%lsf)) then
      do j = 1, this%mesh%ncell_onP
        do k = 1, size(this%lsf,dim=1)
          norm = max(norm, abs(this%lsf(k,j)))
        end do
      end do
    end if
    norm = global_maxval(norm)
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
      if (allocated(this%lsf)) call hash%update(this%lsf(:,:this%mesh%ncell_onP))
    else
      call hash%update(this%lf(:this%mesh%ncell))
      call hash%update(this%hc(:this%mesh%ncell))
      call hash%update(this%tc(:this%mesh%ncell))
      call hash%update(this%tf(:this%mesh%nface))
      if (allocated(this%lsf)) call hash%update(this%lsf)
    end if
    string = hash%hexdigest()
  end function

end module alloy_vector_type
