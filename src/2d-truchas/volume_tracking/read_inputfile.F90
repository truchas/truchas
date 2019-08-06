module read_inputfile

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none

contains

  subroutine readfile(inputfile, xmin, xmax, nx, tsmax, dt, nmat, nvtrack, test_run)

  character(len=*), intent(in) :: inputfile
  real(r8), intent(inout) :: xmin(2), xmax(2), dt
  integer,  intent(inout) :: nx(2), tsmax, nmat, nvtrack, test_run

  integer :: n

  write(*,*) "reading input file: ", inputfile

  open(newunit=n, file=inputfile)
  read(n,*) nx(1), xmin(1), xmax(1)
  read(n,*) nx(2), xmin(2), xmax(2)
  read(n,*) tsmax, dt
  read(n,*) nmat
  read(n,*) nvtrack
  read(n,*) test_run
  close(n)

  end subroutine

end module read_inputfile
