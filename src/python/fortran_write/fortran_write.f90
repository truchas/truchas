! gfortran -fPIC -Wall -shared -o fortran_write.so fortran_write.f90

module fortran_write

  use,intrinsic :: iso_c_binding
  implicit none
  public

contains

  integer(c_int) function fopen(filename) bind(c)
    character(kind=c_char), intent(in) :: filename(*)
    open(newunit=fopen, file=trim(c_to_f_str(filename)), form="unformatted", status="replace")
  end function fopen

  subroutine fclose(unit) bind(c)
    integer(c_int), intent(in), value :: unit
    close(unit)
  end subroutine fclose

  subroutine fwrite_int(unit, x) bind(c)
    integer(c_int), intent(in), value :: unit
    integer(c_int), intent(in), value :: x
    write(unit) x
  end subroutine fwrite_int

  subroutine fwrite_i4x1(unit, x, n) bind(c)
    integer(c_int), intent(in), value :: unit
    integer(c_int), intent(in) :: x(*)
    integer(c_int), intent(in), value :: n
    write(unit) x(:n)
  end subroutine fwrite_i4x1

  subroutine fwrite_r8x1(unit, x, n) bind(c)
    integer(c_int), intent(in), value :: unit
    real(c_double), intent(in) :: x(*)
    integer(c_int), intent(in), value :: n
    write(unit) x(:n)
  end subroutine fwrite_r8x1

  subroutine fwrite_r8(unit, x) bind(c)
    integer(c_int), intent(in), value :: unit
    real(c_double), intent(in), value :: x
    write(unit) x
  end subroutine fwrite_r8

  subroutine fwrite_i8x1(unit, x, n) bind(c)
    integer(c_int), intent(in), value :: unit
    integer(c_int8_t), intent(in) :: x(*)
    integer(c_int), intent(in), value :: n
    write(unit) x(:n)
  end subroutine fwrite_i8x1

  subroutine fwrite_str(unit, x) bind(c)
    integer(c_int), intent(in), value :: unit
    character(kind=c_char), intent(in) :: x(*)
    write(unit) c_to_f_str(x)
  end subroutine fwrite_str


  ! convert a null-terminated C char* to a Fortran character scalar
  function c_to_f_str(cstr) result(fstr)

    character(kind=c_char), intent(in) :: cstr(*)
    character(:), allocatable :: fstr

    integer :: i

    i = 0
    do while (cstr(i+1) /= c_null_char)
      i = i + 1
    end do

    allocate(character(len=i) :: fstr)

    fstr = transfer(cstr(:i), fstr)
    ! do i = 1, len(fstr)
    !   fstr(i:i) = cstr(i)
    ! end do

  end function c_to_f_str

end module fortran_write
