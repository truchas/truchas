!!
!! IH_SOURCE_FACTORY_TYPE
!!
!! This module defines a derived type that encapsulates the parameters defining
!! an external magnetic field, with factory methods for creating abstract
!! function objects for the evaluation of the field. This magnetic field will
!! serve as the driving force in the computation of Joule heat in induction
!! heating (IH) simulations, entering the equations in the nxH boundary condition.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! Refactored February 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The external source is the superposition of a uniform periodic background
!! field and the periodic fields due to a collection of current loops carrying
!! alternating currents, all having the same frequency and phase. Over the
!! much longer time scale of heat transfer, this rapidly varying external
!! source may change at discrete times in "piecewise-constant" manner. Hence
!! the need for this type with factory methods capable of creating abstract
!! field functions on demand at any given induction heating simulation time.
!!

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
    type(induction_coil), allocatable :: coil(:)
    character(1) :: axis
  contains
    procedure :: init
    procedure :: set_time
    procedure :: H_freq
    procedure :: alloc_H_profile_func
    procedure :: alloc_H_waveform_func
    procedure :: source_is_zero
    generic :: source_differs => source_differs1, source_differs2
    procedure, private :: source_differs1, source_differs2
    generic :: source_is_scaled => source_is_scaled1, source_is_scaled2
    procedure, private :: source_is_scaled1, source_is_scaled2
    procedure :: source_data
    procedure :: coil_geom_fingerprint
    procedure, private :: data_index
  end type

  !! Custom SCALAR_FUNC extension that implements the time-periodic
  !! waveform factor for the magnetic field source.
  type, extends(scalar_func), private :: waveform_func
  contains
    procedure :: eval => waveform
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

    !TODO: For case where source-times is not specified should we expect scalar
    ! values for source-frequency, uniform-source, and coil current rather than
    ! size-1 arrays?

    if (params%is_parameter('times')) then
      call params%get('times', this%times, stat, errmsg)
      if (stat /= 0) return
      n = size(this%times)
      if (n > 1) then
        if (any(this%times(2:n) <= this%times(:n-1))) then
          stat = 1
          errmsg = 'times values not strictly increasing'
          return
        end if
      end if
    else
      allocate(this%times(0))
      n = 0
    end if

    call params%get('frequency', this%freq, stat, errmsg)
    if (stat /= 0) return
    if (any(this%freq <= 0.0_r8)) then
      stat = 1
      errmsg = 'frequency is <= 0.0'
      return
    else if (size(this%freq) /= n+1) then
      stat = 1
      errmsg = 'expect ' // i_to_c(n+1) // ' values for frequency'
      return
    endif

    call params%get('uniform-strength', this%const_src, stat, errmsg, default=spread(0.0_r8,1,n+1))
    if (stat /= 0) return
    if (size(this%const_src) /= n+1) then
      stat = 1
      errmsg = 'expect ' // i_to_c(n+1) // ' values for uniform-strength'
      return
    endif

    plist => params%sublist('coils')
    piter = parameter_list_iterator(plist, sublists_only=.true.)
    allocate(this%coil(piter%count()), this%current(piter%count(),n+1))
    do i = 1, size(this%coil)
      plist => piter%sublist()
      associate (coil => this%coil(i))
        call plist%get('center', array, stat, errmsg, default=spread(0.0_r8,1,3))
        if (stat /= 0) return
        if (size(array) /= 3) then
          stat = 1
          errmsg = piter%name() //': center requires a 3-vector value'
          return
        else
          coil%center = array
        end if
        call plist%get('radius', coil%radius, stat, errmsg)
        if (stat /= 0) return
        if (coil%radius <= 0.0_r8) then
          stat = 1
          errmsg = piter%name() // ': radius is <= 0.0'
          return
        end if
        call plist%get('num-loops', coil%nloop, stat, errmsg)
        if (stat /= 0) return
        if (coil%nloop <= 0) then
          stat = 1
          errmsg = piter%name() // ': num-loops is <= 0'
          return
        end if
        if (coil%nloop > 1) then
          call plist%get('length', coil%length, stat, errmsg)
          if (stat /= 0) return
          if (coil%length <= 0) then
            stat = 1
            errmsg = piter%name() // ': length is <= 0'
            return
          end if
        else
          coil%length = 0.0_r8
        end if
        call plist%get('current', array, stat, errmsg)
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

    call params%get('orientation', axis, stat, errmsg, default='Z')
    if (stat /= 0) return
    select case (axis)
    case ('x','X','y','Y','z','Z')
      this%axis = axis
    case default
      stat = 1
      errmsg = 'invalid value for orientation: ' // axis
      return
    end select

  end subroutine init

  !! Set the IH time to T. This effects the results returned by subsequent
  !! calls to H_FREQ, ALLOC_H_WAVEFORM_FUNC, and ALLOC_H_PROFILE_FUNC.

  subroutine set_time(this, t)
    class(ih_source_factory), intent(inout) :: this
    real(r8), intent(in) :: t
    this%t = t
  end subroutine

  !! Return the frequency of the H source at the current IH time.

  real(r8) function H_freq(this)
    class(ih_source_factory), intent(in) :: this
    H_freq = this%freq(this%data_index(this%t))
  end function

  !! Allocate a SCALAR_FUNC class object that is the time-periodic waveform
  !! source factor at the current IH time. NB: This is a 1-periodic function
  !! requiring time scaling / change in time units when setting up the Joule
  !! heat computation using the time domain EM solver.

  subroutine alloc_H_waveform_func(this, f)
    class(ih_source_factory), intent(in) :: this
    class(scalar_func), allocatable, intent(out) :: f
    allocate(waveform_func :: f)
  end subroutine

  !! Allocate a VECTOR_FUNC class object that is the space-dependent profile
  !! source factor at the current IH time. An optional multiplicative scale
  !! factor SCF may be specified.

  subroutine alloc_H_profile_func(this, f, scf)
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

  !! An auxiliary function that locates the IH time interval containing time T.
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

  !! Return an array [FREQ, CONST, CURRENTS] of the parameters defining the
  !! source at the given IH time T. These are the frequency, the strength
  !! of the uniform background field and the currents in the coils. Client
  !! code should use this only as a fingerprint of the source, for use as
  !! input to the SOURCE_DIFFERS and SOURCE_IS_SCALED methods.

  function source_data(this, t) result(src_data)
    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), allocatable :: src_data(:)
    integer :: n
    n = this%data_index(t)
    src_data = [this%freq(n), this%const_src(n), this%current(:,n)]
  end function

  !! Return true if the source is 0 at the given IH time T; otherwise false.

  logical function source_is_zero(this, t)
    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), allocatable :: d(:)
    d = this%source_data(t)
    source_is_zero = all(d(2:) == 0.0_r8)
  end function

  !! Return true if the source at IH times T1 and T2 differ; otherwise false.
  !! This is a specific function for the generic SOURCE_DIFFERS.

  logical function source_differs1(this, t1, t2) result(differs)
    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in) :: t1, t2
    differs = this%source_differs(this%source_data(t1), t2)
  end function

  !! Return true if the source at IH time T2 differs from the source with
  !! fingerprint D1, as returned by SOURCE_DATA. This is a specific function
  !! for the generic SOURCE_DIFFERS.

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

  !! Return true if the source at IH time T2 is a multiple of the source
  !! at IH time T1, and return the scale factor in the argument SCF. This
  !! is a specific function for the generic SOURCE_IS_SCALED.

  logical function source_is_scaled1(this, t1, t2, scf) result(scaled)
    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in)  :: t1, t2
    real(r8), intent(out) :: scf
    scaled = this%source_is_scaled(this%source_data(t1), t2, scf)
  end function source_is_scaled1

  !! Return true if the source strength at IH time T2 is a multiple of the
  !! strength of the source with fingerprint D1, as returned by SOURCE_DATA,
  !! and return the scale factor in the argument SCF. This is a specific
  !! function for the generic SOURCE_IS_SCALED. NB: A relatively arbitrary
  !! tolerance is used to decide between scaled and not-scaled.

  logical function source_is_scaled2(this, d1, t2, scf) result(scaled)

    class(ih_source_factory), intent(in) :: this
    real(r8), intent(in)  :: d1(:), t2
    real(r8), intent(out) :: scf

    real(r8) :: a, err
    real(r8), allocatable :: d2(:)

    d2 = this%source_data(t2)

    associate (f1 => d1(1), u1 => d1(2:), f2 => d2(1), u2 => d2(2:))
      ! Best scale factor in least-squares sense
      a = norm2(u1)
      if (a > 0.0_r8) then
        scf = dot_product(u1,u2) / a**2
      else
        scf = 0.0_r8
      end if
      err = norm2(u2-scf*u1) ! l2 error in best scaling
      scaled = (f1 == f2) .and. (err <= a*1.0e-6)
    end associate

  end function source_is_scaled2

  !! Returns the MD5 checksum of the source parameters that define the fixed
  !! geometry of the source which are invariant in an induction heating
  !! simulation. These are the symmetry axis and the parameters defining the
  !! geometry of the coils. This is useful for tagging joule heat data in
  !! checkpoints to ensure that the source geometry in restarts has not
  !! changed.

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

  !! The type-bound EVAL function for the WAVEFORM_FUNC extension of the
  !! SCALAR_FUNC class. This is a 1-periodic waveform function with initial
  !! fade-in to full strength. There are several potential options for the
  !! waveform, but it is currently hardwired to a simple sinusoid.

  function waveform(this, x) result(fx)
    class(waveform_func), intent(in) :: this
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
