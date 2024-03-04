module ih_source_factory_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_class
  use vector_func_class
  use ih_source_func_type
  implicit none
  private

  type, public :: ih_source_factory
    private
    real(r8) :: t
    real(r8), allocatable :: times(:), freq(:), const_src(:), current(:,:)
    type(solenoid), allocatable :: coil(:)
    character(1) :: axis
  contains
    procedure :: init
    procedure :: set_time
    procedure :: source_frequency
    procedure :: alloc_src_mod_func
    procedure :: alloc_ih_source_func
    procedure :: source_is_zero
    generic :: source_differs => source_differs1, source_differs2
    procedure, private :: source_differs1, source_differs2
    generic :: source_is_scaled => source_is_scaled1, source_is_scaled2
    procedure, private :: source_is_scaled1, source_is_scaled2
    procedure :: source_data
    procedure :: coil_geom_fingerprint
    procedure, private :: data_index
  end type

  type, extends(scalar_func), private :: src_mod_func
  contains
    procedure :: eval => src_mod
  end type

contains

  subroutine init(this, params, t, stat, errmsg)

    use parameter_list_type
    use string_utilities, only: i_to_c

    class(ih_source_factory), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    real(r8), intent(in) :: t
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, n
    real(r8), allocatable :: array(:)
    type(parameter_list), pointer :: plist
    type(parameter_list_iterator) :: piter
    character(:), allocatable :: axis

    this%t = t

    !TODO: For case where source-times is not specified we should expect scalar
    ! values for source-frequency, uniform-source, and coil current rather than
    ! size-1 arrays.

    if (params%is_parameter('source-times')) then
      call params%get('source-times', this%times, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      n = size(this%times)
      if (n > 1) then
        if (any(this%times(2:n) <= this%times(:n-1))) then
          stat = 1
          errmsg = 'source-times values not strictly increasing'
          return
        end if
      end if
    else
      allocate(this%times(0))
      n = 0
    end if

    call params%get('source-frequency', this%freq, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (any(this%freq <= 0.0_r8)) then
      stat = 1
      errmsg = 'source-frequency is <= 0.0'
      return
    else if (size(this%freq) /= n+1) then
      stat = 1
      errmsg = 'expect ' // i_to_c(n+1) // ' values for source-frequency'
      return
    endif

    call params%get('uniform-source', this%const_src, stat=stat, errmsg=errmsg, default=spread(0.0_r8,1,n+1))
    if (stat /= 0) return
    if (size(this%const_src) /= n+1) then
      stat = 1
      errmsg = 'expect ' // i_to_c(n+1) // ' values for uniform-source'
      return
    endif

    plist => params%sublist('coils')
    piter = parameter_list_iterator(plist, sublists_only=.true.)
    allocate(this%coil(piter%count()), this%current(piter%count(),n+1))
    do i = 1, size(this%coil)
      plist => piter%sublist()
      associate (coil => this%coil(i))
        call plist%get('center', array, stat=stat, errmsg=errmsg, default=spread(0.0_r8,1,3))
        if (stat /= 0) return
        if (size(array) /= 3) then
          stat = 1
          errmsg = piter%name() //': center requires a 3-vector value'
          return
        else
          coil%center = array
        end if
        call plist%get('radius', coil%radius, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        if (coil%radius <= 0.0_r8) then
          stat = 1
          errmsg = piter%name() // ': radius is <= 0.0'
          return
        end if
        call plist%get('num-turns', coil%nloop, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        if (coil%nloop <= 0) then
          stat = 1
          errmsg = piter%name() // ': num-turns is <= 0'
          return
        end if
        if (coil%nloop > 1) then
          call plist%get('length', coil%length, stat=stat, errmsg=errmsg)
          if (stat /= 0) return
          if (coil%length <= 0) then
            stat = 1
            errmsg = piter%name() // ': length is <= 0'
            return
          end if
        else
          coil%length = 0.0_r8
        end if
        call plist%get('current', array, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        if (size(array) /= n+1) then
          stat = 1
          errmsg = piter%name() // ': expect ' // i_to_c(n+1) // ' values for current'
          return
        end if
        this%current(i,:) = array
      end associate
      call piter%next
    end do

    call params%get('symmetry-axis', axis, stat=stat, errmsg=errmsg, default='Z')
    if (stat /= 0) return
    select case (axis)
    case ('x','X','y','Y','z','Z')
      this%axis = axis
    case default
      stat = 1
      errmsg = 'invalid value for symmetry-axis: ' // axis
      return
    end select

  end subroutine

  subroutine set_time(this, t)
    class(ih_source_factory), intent(inout) :: this
    real(r8), intent(in) :: t
    this%t = t
  end subroutine

  real(r8) function source_frequency(this) result(freq)
    class(ih_source_factory), intent(in) :: this
    freq = this%freq(this%data_index(this%t))
  end function

  subroutine alloc_src_mod_func(this, f)
    class(ih_source_factory), intent(in) :: this
    class(scalar_func), allocatable, intent(out) :: f
    allocate(src_mod_func :: f)
  end subroutine

  subroutine alloc_ih_source_func(this, f, scf)
    class(ih_source_factory), intent(in) :: this
    class(vector_func), allocatable, intent(out) :: f
    real(r8), intent(in), optional :: scf
    type(ih_source_func), allocatable :: h
    integer :: i
    i = this%data_index(this%t)
    allocate(h)
    call h%init(this%axis, this%const_src(i), this%coil, this%current(:,i), scf)
    call move_alloc(h, f)
  end subroutine

  integer function data_index(this, t) result(n)
    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in) :: t
    n = size(this%times)
    do while (n > 0)
      if (t >= this%times(n)) exit
      n = n - 1
    end do
    n = n + 1
  end function

  function source_data(this, t) result(src_data)
    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), allocatable :: src_data(:)
    integer :: n
    n = this%data_index(t)
    src_data = [this%freq(n), this%const_src(n), this%current(:,n)]
  end function

  logical function source_is_zero(this, t)
    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), allocatable :: d(:)
    d = this%source_data(t)
    source_is_zero = all(d(2:) == 0.0_r8)
  end function

  logical function source_differs1(this, t1, t2) result(differs)
    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in) :: t1, t2
    differs = this%source_differs(this%source_data(t1), t2)
  end function

  logical function source_differs2(this, d1, t2) result(differs)
    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in) :: d1(:), t2
    real(r8), allocatable :: d2(:)
    d2 = this%source_data(t2)
    associate (f1 => d1(1), u1 => d1(2:), f2 => d2(1), u2 => d2(2:))
      differs = .true.
      if (any(u1 /= u2)) return
      if (any(u1 /= 0.0_r8) .and. f1 /= f2) return
      differs = .false.
    end associate
  end function

  logical function source_is_scaled1(this, t1, t2, s) result(scaled)
    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in)  :: t1, t2
    real(r8), intent(out) :: s
    scaled = this%source_is_scaled(this%source_data(t1), t2, s)
  end function source_is_scaled1

  logical function source_is_scaled2(this, d1, t2, s) result(scaled)

    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in)  :: d1(:), t2
    real(r8), intent(out) :: s

    real(r8) :: a, err
    real(r8), allocatable :: d2(:)

    d2 = this%source_data(t2)

    associate (f1 => d1(1), u1 => d1(2:), f2 => d2(1), u2 => d2(2:))
      !! Best scale factor in least-squares sense
      a = norm2(u1)
      if (a > 0.0_r8) then
        s = dot_product(u1,u2) / a**2
      else
        s = 0.0_r8
      end if
      err = norm2(u2-s*u1) ! l2 error in best scaling
      scaled = (f1 == f2) .and. (err <= a*1.0e-6)
    end associate

  end function source_is_scaled2

  function coil_geom_fingerprint(this) result(fp)
    use secure_hash_factory
    class(ih_source_factory), intent(in) :: this
    character(:), allocatable :: fp
    integer :: n
    class(secure_hash), allocatable :: hash
    call new_secure_hash(hash, 'md5')
    do n = 1, size(this%coil)
      call hash%update(this%coil(n)%length)
      call hash%update(this%coil(n)%radius)
      call hash%update(this%coil(n)%center)
      call hash%update(this%coil(n)%nloop)
    end do
    call hash%update(this%axis)
    fp = hash%hexdigest()
  end function

  function src_mod(this, x) result(fx)
    class(src_mod_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    real(r8), parameter :: PI =    3.1415926535897932385_r8
    real(r8), parameter :: TWOPI = 6.2831853071795864769_r8
    associate (t => x(1))
      select case (0)
      case (1)  ! Truncated l2 fit to a square wave
        fx = (sin(TWOPI*t)+sin(3*TWOPI*t)/3.0+sin(5*TWOPI*t)/5.0+sin(7*TWOPI*t)/7.0)*(4.0/PI)
      case (2)  ! Non-oscillatory 'square' wave
        fx = (1225*sin(TWOPI*t)+245*sin(3*TWOPI*t)+49*sin(5*TWOPI*t)+5*sin(7*TWOPI*t))/1024.0
      case default ! Basic sinusoidal wave form.
        fx = sin(TWOPI*t)
      end select
      fx = (1.0_r8 - exp(-2.0_r8*t**2))*fx
    end associate
  end function

end module ih_source_factory_type
