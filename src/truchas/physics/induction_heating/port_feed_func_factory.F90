module port_feed_func_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use complex_scalar_func_class
  use complex_vector_func_class
  use physical_constants, only: vacuum_permittivity, vacuum_permeability
  implicit none
  private

  public :: alloc_te01_port_feed_func

  real(r8), parameter :: PI = 3.1415926535897932385_r8

contains

  subroutine alloc_te01_port_feed_func(omega, params, alpha, g, stat, errmsg)

    use const_complex_scalar_func_type
    use fptr_complex_vector_func_type

    real(r8), intent(in) :: omega
    type(parameter_list), intent(inout) :: params
    class(complex_scalar_func), allocatable, intent(out) :: alpha
    class(complex_vector_func), allocatable, intent(out) :: g
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: rvec(:)
    real(r8) :: x0(3), x_axis(3), y_axis(3), a, b, power
    real(r8) :: k0, Z0, h0, E0, p(10)

    call params%get('center', rvec, stat, errmsg)
    if (stat /= 0) return
    if (size(rvec) /= 3) then
      stat = 1
      errmsg = 'center is not a 3-vector'
      return
    end if
    x0 = rvec

    call params%get('x-axis', rvec, stat, errmsg)
    if (stat /= 0) return
    if (size(rvec) /= 3) then
      stat = 1
      errmsg = 'x-axis not a 3-vector'
      return
    end if
    if (norm2(rvec) == 0) then
      stat = 1
      errmsg = 'x-axis is 0'
      return
    end if
    x_axis = rvec / norm2(rvec)

    call params%get('y-axis', rvec, stat, errmsg)
    if (stat /= 0) return
    if (size(rvec) /= 3) then
      stat = 1
      errmsg = 'y-axis not a 3-vector'
      return
    end if
    rvec = rvec - dot_product(rvec, x_axis)*x_axis
    if (norm2(rvec) == 0) then
      stat = 1
      errmsg = 'y-axis not independent of x-axis'
      return
    end if
    y_axis = rvec / norm2(rvec)

    call params%get('x-width', a, stat, errmsg)
    if (stat /= 0) return
    if (a <= 0) then
      stat = 1
      errmsg = 'x-width is <= 0.0'
      return
    end if

    call params%get('y-width', b, stat, errmsg)
    if (stat /= 0) return
    if (b <= 0) then
      stat = 1
      errmsg = 'y-width is <= 0.0'
      return
    end if

    k0 = omega * sqrt(vacuum_permittivity*vacuum_permeability)
    h0 = sqrt(k0**2 - (PI/a)**2)
    call alloc_const_complex_scalar_func(alpha, cmplx(0.0_r8, h0, kind=r8))

    if (params%is_parameter('power')) then
      if (params%is_parameter('E-mag')) then
        stat = 1
        errmsg = 'parameters power and e-mag are mutually exclusive'
        return
      end if
      call params%get('power', power, stat, errmsg)
      if (stat /= 0) return
      if (power < 0) then
        stat = 1
        errmsg = 'power is < 0.0'
        return
      end if
      Z0 = sqrt(vacuum_permeability/vacuum_permittivity)
      E0 = sqrt((4*power/(a*b)) * (k0*Z0/h0))
    else
      call params%get('e-mag', E0, stat, errmsg, default=1.0_r8)
      if (stat /= 0) return
    end if

    if (E0 == 0) return ! no g function
    p = [x0, (PI/a)*x_axis, y_axis, 2*h0*E0] ! function parameter array
    call alloc_fptr_complex_vector_func(g, 3, te01_mode, p)

  end subroutine alloc_te01_port_feed_func

  function te01_mode(x, p, dim) result(fx)
    real(r8), intent(in) :: x(*), p(*)
    integer, value :: dim
    complex(r8) :: fx(dim)
    fx%re = 0.0_r8
    fx%im = (p(10) * cos(dot_product(x(1:3)-p(1:3),p(4:6)))) * p(7:9)
  end function

end module port_feed_func_factory
