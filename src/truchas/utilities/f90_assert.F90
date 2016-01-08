!!
!!  F90_ASSERT -- C-style assertions for Fortran.
!!
!!    Neil N. Carlson <nnc@newmexico.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  Usage: At the top of the source file, include the preprocessor file
!!  f90_assert.fpp which defines the Assert() preprocessor macro:
!!
!!    #include "f90_assert.fpp"
!!
!!  If the macro NDEBUG is undefined when the file is passed through the
!!  preprocessor, lines of the form
!!
!!    ASSERT( <scalar logical expression> )
!!
!!  will be expanded to Fortran code which tests whether the logical
!!  expression is true, and if not, calls the following routine which
!!  will print the file name and line number and then halt execution.
!!  If the macro NDEBUG is defined (e.g., -D NDEBUG), then the ASSERT()
!!  is expanded to a Fortran comment line.
!!
!!  This is intentionally not a module procedure.
!!
!!  NB: Use with Fortran-aware preprocessors like fpp is robust.  One
!!  can use the C preprocessor cpp, but if the expanded macro extends
!!  the line past 132 characters, a compiler error will probably result.
!!

subroutine f90_assert (file, line)

  use truchas_logging_services
  
  character(len=*), intent(in) :: file
  integer,          intent(in) :: line

  character(len=256) :: string
  
  write(string,fmt='(a,i4.4)') "Assertion failed at " // file // ":", line
  call TLS_panic (string)

end subroutine f90_assert
