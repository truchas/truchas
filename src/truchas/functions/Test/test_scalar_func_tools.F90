program test_scalar_func_tools

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_factories
  use scalar_func_tools
  implicit none

  integer :: status = 0

  call test_poly_antideriv

  if (status /= 0) error stop

contains

  subroutine test_poly_antideriv

    class(scalar_func), allocatable :: f, g
    integer :: stat
    character(:), allocatable :: errmsg

    call alloc_poly_scalar_func(f, [1.0_r8,3.0_r8], [-2,0], -1.0_r8)
    call alloc_scalar_func_antideriv(f, 1.0_r8, 6.0_r8, g, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test_alloc_poly_antideriv: allocation failed: ' // errmsg)
      return
    end if

    if (g%eval([1.0_r8]) /= 6.0_r8)  call write_fail('test_alloc_poly_antideriv: bad value at 1')
    if (g%eval([0.0_r8]) /= 2.5_r8)  call write_fail('test_alloc_poly_antideriv: bad value at 0')

  end subroutine

  subroutine write_fail(errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program
