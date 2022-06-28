!!
!! INPUT_UTILITIES
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This module provides some utilities that are useful in reading input files:
!!
!!    CALL SEEK_TO_NAMELIST (UNIT, NML, FOUND [, IOSTAT]) seeks to the next
!!      record in the file connected to UNIT that contains the start of a
!!      namelist data group for the namelist named NML.  This is a record whose
!!      first non white space character is an '&' followed (without intervening
!!      white space) by the name of the namelist.  Searching starts from the
!!      current position, and if such a record is found, the file is left
!!      positioned at that record and FOUND is set to true; otherwise FOUND is
!!      set to false.  The search is not case sensitive.  If an I/O error
!!      occurs, the nonzero I/O error code is returned in IOSTAT if it is
!!      present, otherwise an error message is written to stderr and the program
!!      is halted.
!!
!!      Note that the file connected to UNIT must be capable of being backspaced
!!      which is not the case, for example, when reading from the tty.
!!

module input_utilities

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  public :: seek_to_namelist

  integer,   parameter, public :: NULL_I = huge(1)
  real(r8),  parameter, public :: NULL_R = huge(1.0_r8)
  character, parameter, public :: NULL_C = char(0)

  !! White space characters: blank, tab
  character, parameter, private :: WHITE_SPACE(2) = [' ', achar(9)]

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SEEK_TO_NAMELIST
 !!
 !! This subroutine seeks to the next record in the file connected to UNIT that
 !! contains the start of a namelist data group for the namelist named NML.
 !! This is a record whose first non white space character is an '&' followed
 !! (without intervening white space) by the name of the namelist.  Searching
 !! starts from the current position, and if such a record is found, the file is
 !! left positioned at that record and FOUND is set to true; otherwise FOUND is
 !! set to false.  The search is not case sensitive.
 !!
 !! Implementation Note
 !!
 !! o If searching for the namelist 'foo', we need to be sure that we don't
 !!   mistakenly match '& foo' or '&foobar', thus we do a formatted read of a
 !!   string one character longer than the namelist name.  A matching input
 !!   string will be right delimited with a white space or EOR, and in the later
 !!   case the string read will be right padded with a blank.  For consistency,
 !!   we replace the rightmost character with a blank if it is a white space.
 !!   Thus a valid input will correctly match the namelist name during string
 !!   comparison.
 !!

  subroutine seek_to_namelist (unit, nml, found, iostat)

    use string_utilities, only: raise_case, i_to_c

    integer, intent(in) :: unit
    character(len=*), intent(in) :: nml
    logical, intent(out) :: found
    integer, intent(out), optional :: iostat

    integer :: ios
    integer :: n
    character(len=1) :: c
    character(len=1+len_trim(adjustl(nml))) :: name, string

    n = len(string)
    found = .false.
    if (present(iostat)) iostat = 0

    name = raise_case(adjustl(nml))
    !! Make sure we have a real name to search for.
    if (len_trim(name) == 0) return

    !!! Should we make sure that we are positioned at the start of a record? !!!

    do
      read(unit,fmt='(a)',advance='no',iostat=ios,END=999) c
      if (ios > 0) exit     ! unexpected read error; quit, returning false.
      if (ios < 0) cycle    ! EOR; file now positioned at the start of the next record.
      if (any(c == WHITE_SPACE)) cycle   ! leading white space; read the next character in this record.
      if (c == '&') then    ! possibly the start of the namelist; check the trailing string.
        read(unit,fmt='(a)',iostat=ios) string
        if (ios /= 0) cycle ! read error; not the namelist we're looking for.
        if (any(string(n:n) == WHITE_SPACE)) string(n:n) = ' '
        if (raise_case(string) /= name) cycle
        found = .true.
        backspace(unit)     ! leave file positioned ready to read the namelist.
        return
      else
        read(unit,*)
        cycle
      end if
    end do

    if (present(iostat)) then
      iostat = ios
    else if (ios /= 0) then
      write(0,'(4a)') 'SEEK_TO_NAMELIST: I/O error on unit ', i_to_c(unit), ': IOSTAT=', i_to_c(ios)
      stop
    end if

    999 continue

  end subroutine seek_to_namelist

end module input_utilities
