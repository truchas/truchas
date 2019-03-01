!!
!! TRUCHAS_ENV
!!
!! This evolving module is intended to collect data and procedures that are
!! associated with the run-time environment of Truchas -- things like file
!! and directory paths.  Much of this is found in output_module, which is
!! slated for removal, and other places.  I expect this to eventually morph
!! into something reasonable.
!!
!! Neil Carlson <nnc@lanl.gov>
!! October 2012
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module truchas_env

  implicit none
  private
  
  public :: output_file_name, new_unit
  
  ! Stuff formerly in output_data_dir
  character(256), public :: prefix, input_dir, output_dir, title, input_file

contains

  !! A simplified replacement for make_file_name.  Given a string EXT, it
  !! returns the path of the output file in the output directory having
  !! EXT as the extention.  The base name of the file is determined by the
  !! input file name.  This does no more than concatenate strings -- I don't
  !! know how necessary this really is.

  function output_file_name (ext) result (ofile)
    character(*), intent(in) :: ext
    character(len_trim(prefix)+len_trim(ext)+1) :: ofile
    ofile = trim(prefix) // '.' // trim(ext)
  end function output_file_name

  !! Returns an unused unit number, or -1 if none was found.  This call does
  !! not reserve the returned unit number, so a file should be connected to
  !! the unit promptly to prevent another call to this procedure returning
  !! the same unit number.  Further, there is nothing to prevent other code,
  !! such as third-party libraries, attempting to use the same unit number for
  !! their own purposes (with likely disastrous affect).  When possible, code
  !! should instead use the NEWUNIT specifier on the OPEN statement (introduced
  !! in the Fortran 2008 standard) which is guaranteed to give a unit number
  !  that will not clash with any user-specified unit.

  subroutine new_unit (unit)

    use,intrinsic :: iso_fortran_env, only: input_unit, output_unit, error_unit

    integer, intent(out) :: unit

    integer, parameter :: PRECONNECTED_UNITS(3) = (/input_unit, output_unit, error_unit/)
    integer, parameter :: MAX_UNIT_NUMBER = 127 ! TBrook starts assigning at 128

    integer :: ios
    logical :: exists, opened

    do unit = 1, MAX_UNIT_NUMBER
      if (any(unit == PRECONNECTED_UNITS)) cycle
      inquire(unit, exist=exists, opened=opened, iostat=ios)
      if (exists .and. .not.opened .and. ios == 0) return
    end do

    unit = -1

  end subroutine new_unit

end module truchas_env
