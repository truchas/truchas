!!
!! Aditya K. Pandare <apandare@lanl.gov>, January 2020
!! SPDX-License-Identifier: BSD-3-Clause
!!
module read_inputfile

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none

contains

  subroutine readfile(inputfile, xmin, xmax, nx, dxeps, ptri, tsmax, dt, nmat, nvtrack)

  character(len=*), intent(in) :: inputfile
  real(r8), intent(inout) :: xmin(2), xmax(2), dxeps, ptri, dt
  integer,  intent(inout) :: nx(2), tsmax, nmat, nvtrack

  integer :: n

  write(*,*) "reading input file: ", inputfile

  open(newunit=n, file=inputfile)
  read(n,*) nx(1), xmin(1), xmax(1)
  read(n,*) nx(2), xmin(2), xmax(2)
  read(n,*) dxeps, ptri
  read(n,*) tsmax, dt
  read(n,*) nmat
  read(n,*) nvtrack
  close(n)

  end subroutine

end module read_inputfile
