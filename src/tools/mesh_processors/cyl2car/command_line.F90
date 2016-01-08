!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module command_line

  implicit none
  private

  public  :: parse_command_line
  private :: usage_halt, usage_summary

  !! Public variables initialized from the command line.
  character(len=32), public :: prog

contains

  subroutine parse_command_line (infile, outfile)

    character(len=*), intent(out) :: infile, outfile

    integer :: i, n, num_arg
    character(len=128) :: arg

    call get_command_argument (0, arg)
    i = scan(arg, '/', back=.true.)
    prog = arg(i+1:)  ! Remove the leading path component, if any.

   !!!
   !!! PARSE COMMAND-LINE OPTIONS

    n = 0
    num_arg = command_argument_count()
    do while (n < num_arg)

      n = n + 1
      call get_command_argument (n, arg)

      select case (arg)

      case ('-h','--help')

        call usage_summary ()

      !!
      !! Example option with following numeric value
      !!
      !! case ('-n')
      !! 
      !!   n = n + 1
      !!   if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
      !!   call get_command_argument (n, arg)
      !!   read(unit=arg,fmt=*,iostat=ios) n_value
      !!   if (ios /= 0) call usage_halt ('invalid value for -n: ' // trim(arg))
      !!
      !! Continue to do sanity checking of the value.
      !!

      case default

        if (arg(1:1) == '-') call usage_halt ('invalid option: ' // trim(arg))
        n = n - 1 ! Back-up.
        exit

      end select

    end do

   !!!
   !!! INTERPRET REMAINING ARGUMENTS AS FILE NAMES

    if (num_arg - n /= 2) call usage_halt ('exactly two file arguments are required')
    call get_command_argument (n+1, infile)
    call get_command_argument (n+2, outfile)

  end subroutine parse_command_line

  subroutine usage_halt (message)

    character(len=*), intent(in) :: message

    write(unit=*, fmt='(a)') trim(prog) // ': ' // message
    write(unit=*, fmt='(a)') 'Try `' // trim(prog) // ' --help'' for more information.'
    stop

  end subroutine usage_halt

  subroutine usage_summary ()

    write(unit=*,fmt='(a)') 'Usage: ' // trim(prog) // ' [options] INFILE OUTFILE'
    write(unit=*,fmt='(a)') 'Map the Exodus-format mesh INFILE with (r,theta,z) cylindrical-coordinate'
    write(unit=*,fmt='(a)') 'data to an Exodus-format mesh OUTFILE with Cartesian-coordinate data,'
    write(unit=*,fmt='(a)') 'transforming element types as necessary.  The theta coordinate should be'
    write(unit=*,fmt='(a)') 'given in degrees.'
    write(unit=*,fmt='(a)') ' '
    write(unit=*,fmt='(a)') 'Options:'
    write(unit=*,fmt='(a)') '  -h, --help                  Display this help and exit.'
    stop

  end subroutine usage_summary

end module command_line
